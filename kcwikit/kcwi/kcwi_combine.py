"""
kcwi_combine.py -- align and co-add cubes that are *already stacked* by
kcwi.kcwi_stack (e.g. the same field from two runs, or two overlapping
pointings).

Inputs are the stacked products written by kcwi_stack:
    <name>_i<suffix>.fits   intensity   (surface brightness units)
    <name>_v<suffix>.fits   variance
    <name>_m<suffix>.fits   mask        (optional; rebuilt from the exposure map)
    <name>_e<suffix>.fits   exposure map
where <suffix> is e.g. 'cubes', 'cubed', 'cubes_l'. The '_wcs' intensity cube
written by kcwi_astrometry is also accepted; its WCS is then used for the
v/m/e cubes too (kcwi_astrometry does not write _wcs versions of them).

Steps:
    1. Collapse every cube into a white-light image over `wavebin`.
    2. Cross-correlate each image against the reference (the first cube) with a
       zero-normalized (Pearson) correlation over the pixels where both images
       overlap. The trial shifts are applied by resampling through the WCS, so
       different pixel scales or orientations are handled. There is a coarse
       integer-pixel pass, then a sub-pixel pass (1/upfactor).
    3. Correct CRVAL1/2 of each non-reference cube.
    4. Resample the wavelength axis to the reference grid if needed, then
       resample each cube spatially (NaN-aware bilinear) onto the output grid
       ('ref' = the reference cube's grid, 'union' = extended to cover all cubes).
    5. Co-add with exposure-map weights (as kcwi_stack does):
           I = sum(w I)/sum(w),  V = sum(w^2 V)/sum(w)^2,  E = sum(E)
       or with inverse-variance weights. Optionally (crr=True), cosmic rays
       and other positive outliers are rejected first by comparing each cube
       with the others on the common grid, using an absolute (noise) test AND a
       relative (source flux) test. See _flag_outliers.

Only numpy/scipy/astropy/matplotlib are required.
"""

import os
import re
import warnings
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs
from astropy import table
from scipy import ndimage
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _stack_filenames(icubefn):
    """Return the i/v/m/e file names of a kcwi_stack product."""
    base = os.path.basename(icubefn)
    m = re.search(r'_i(cube[sd](?:_[a-z])?)(_wcs)?\.fits$', base)
    if m is None:
        raise ValueError('Cannot parse stacked-cube name: {}. Expected '
                         '*_icubes.fits, *_icubed.fits, *_icubes_l.fits, ...'.format(icubefn))
    suffix = m.group(1)
    stem = icubefn[:len(icubefn) - len(m.group(0))]
    out = {'i': icubefn, 'suffix': suffix}
    for key in ['v', 'm', 'e']:
        fn = stem + '_' + key + suffix + '.fits'
        out[key] = fn if os.path.isfile(fn) else None
    return out


def _wave_axis(hdr):
    """Wavelength axis in header units (Angstrom for KCWI)."""
    cd = hdr['CD3_3'] if 'CD3_3' in hdr else hdr['CDELT3']
    return hdr['CRVAL3'] + (np.arange(hdr['NAXIS3']) + 1 - hdr['CRPIX3']) * cd


def _same_wave(h0, h1):
    if h0['NAXIS3'] != h1['NAXIS3']:
        return False
    w0, w1 = _wave_axis(h0), _wave_axis(h1)
    return np.allclose(w0, w1, rtol=0, atol=1e-3 * abs(w0[1] - w0[0]))


def _read_stack(icubefn):
    """Read i/v/e cubes (FITS order: [nz, ny, nx]) and build a validity mask."""
    fns = _stack_filenames(icubefn)
    hdu_i = fits.open(fns['i'])[0]
    hdr = hdu_i.header.copy()
    cube = hdu_i.data.astype(np.float64)

    if fns['v'] is None:
        raise FileNotFoundError('Variance cube not found for {}'.format(icubefn))
    var = fits.getdata(fns['v']).astype(np.float64)

    if fns['e'] is not None:
        exp = fits.getdata(fns['e']).astype(np.float64)
    else:
        print('  [warning] no exposure cube for {}; using unit exposure'.format(icubefn))
        exp = np.ones_like(cube)

    valid = np.isfinite(cube) & np.isfinite(var) & (var > 0) & (exp > 0)
    if fns['m'] is not None:
        valid &= (fits.getdata(fns['m']) == 0)

    cube[~valid] = np.nan
    var[~valid] = np.nan
    exp[~valid] = 0.
    return dict(fn=fns, hdr=hdr, i=cube, v=var, e=exp, valid=valid)


def _whitelight(st, wavebin):
    wave = _wave_axis(st['hdr'])
    q = (wave > wavebin[0]) & (wave < wavebin[1])
    if q.sum() == 0:
        raise ValueError('wavebin {} outside the wavelength range of {}'.format(
            wavebin, st['fn']['i']))
    sub = st['i'][q]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        img = np.nanmean(sub, axis=0)
    return img                                     # [ny, nx], NaN where no data


def _shift_map(w_out, shape_out, w_in, d=(0., 0.)):
    """Pixel coordinates in the input image of every output pixel, with the
    output grid shifted by d = (dx, dy) output pixels (content moves by +d)."""
    ny, nx = shape_out
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    ra, dec = w_out.wcs_pix2world(xx - d[0], yy - d[1], 0)
    xi, yi = w_in.wcs_world2pix(ra, dec, 0)
    return yi, xi


def _resample2d(img, yi, xi, order=1, minvalid=0.999):
    """NaN-aware resampling of a 2D image to the coordinates (yi, xi)."""
    good = np.isfinite(img)
    num = ndimage.map_coordinates(np.where(good, img, 0.), [yi, xi], order=order,
                                  mode='constant', cval=0., prefilter=False)
    den = ndimage.map_coordinates(good.astype(float), [yi, xi], order=order,
                                  mode='constant', cval=0., prefilter=False)
    out = np.full(num.shape, np.nan)
    ok = den >= minvalid
    out[ok] = num[ok] / den[ok]
    return out


def _ncc(a, b, box=None, min_frac=0.3):
    """Zero-normalized (Pearson) cross-correlation over the pixels where both are finite."""
    if box is not None:
        a = a[box[2]:box[3], box[0]:box[1]]
        b = b[box[2]:box[3], box[0]:box[1]]
    q = np.isfinite(a) & np.isfinite(b)
    if q.sum() < min_frac * np.isfinite(a).sum() or q.sum() < 10:
        return np.nan
    aa = a[q] - a[q].mean()
    bb = b[q] - b[q].mean()
    den = np.sqrt(np.sum(aa**2) * np.sum(bb**2))
    return np.sum(aa * bb) / den if den > 0 else np.nan


def _xcorr_map(img0, w0, img, w, dxs, dys, box, min_frac):
    crl = np.full((len(dxs), len(dys)), np.nan)
    for ii, dx in enumerate(dxs):
        for jj, dy in enumerate(dys):
            yi, xi = _shift_map(w0, img0.shape, w, (dx, dy))
            crl[ii, jj] = _ncc(img0, _resample2d(img, yi, xi), box, min_frac)
    return crl                                       # indexed [x, y]


def _resample_wave(st, hdr_ref, kind='cubic'):
    """Resample i/v/e onto the reference wavelength grid (NaN-aware)."""
    w_in, w_out = _wave_axis(st['hdr']), _wave_axis(hdr_ref)
    nz, ny, nx = st['i'].shape
    good = st['valid'].reshape(nz, -1).astype(float)

    def _interp(arr, k, fill=0.):
        a = np.nan_to_num(arr.reshape(nz, -1))
        f = interpolate.interp1d(w_in, a, kind=k, axis=0, bounds_error=False, fill_value=fill)
        return f(w_out)

    gm = _interp(good, 'linear', 0.)
    ok = gm >= 0.999
    i_new = np.where(ok, _interp(st['i'], kind), np.nan)
    v_new = np.where(ok, _interp(st['v'], 'linear'), np.nan)   # linear keeps V >= 0
    e_new = np.where(ok, _interp(st['e'], 'linear'), 0.)
    shp = (len(w_out), ny, nx)
    st['i'], st['v'], st['e'] = i_new.reshape(shp), v_new.reshape(shp), e_new.reshape(shp)
    st['valid'] = ok.reshape(shp)
    for key in ['NAXIS3', 'CRVAL3', 'CRPIX3', 'CD3_3']:
        if key in hdr_ref:
            st['hdr'][key] = hdr_ref[key]
    return st


def _resample_block(st, yi, xi, z0, z1, shape_out):
    """Spatially resample slices z0:z1 of i/v/e onto the output grid."""
    nzb = z1 - z0
    i_out = np.full((nzb,) + tuple(shape_out), np.nan)
    v_out = np.full_like(i_out, np.nan)
    e_out = np.zeros_like(i_out)
    for kk, k in enumerate(range(z0, z1)):
        if not st['valid'][k].any():
            continue
        i_out[kk] = _resample2d(st['i'][k], yi, xi)
        v_out[kk] = _resample2d(st['v'][k], yi, xi)
        ek = np.where(st['valid'][k], st['e'][k], np.nan)
        e_out[kk] = np.nan_to_num(_resample2d(ek, yi, xi))
    bad = ~np.isfinite(i_out) | ~np.isfinite(v_out)
    i_out[bad], v_out[bad], e_out[bad] = np.nan, np.nan, 0.
    return i_out, v_out, e_out


def _crr_zsmooth(a, size, step=1):
    """Running median along wavelength, evaluated every `step` channels.

    Used for the per-spaxel scale between two cubes, which is smooth along
    wavelength by construction, so there is no point evaluating it at every
    channel: scipy's median_filter costs O(size) per output voxel and is ~80%
    of the cosmic-ray rejection's run time on a long stack. Evaluating it on a
    decimated grid and interpolating is about eight times faster and agrees to
    ~1%, which then largely cancels in the ratio of two such medians.

    Args:
        a (ndarray): [nz, ny, nx].
        size (int): width of the running median [channels].
        step (int): evaluate every this many channels; 1 reproduces
            ndimage.median_filter(a, (size, 1, 1), mode='nearest') exactly.

    Returns:
        ndarray: same shape as `a`.
    """
    nz = a.shape[0]
    step = max(int(step), 1)
    if step == 1 or nz < 2 * step:
        return ndimage.median_filter(a, size=(int(size), 1, 1), mode='nearest')
    h = int(size) // 2
    zc = np.arange(0, nz, step)
    if zc[-1] != nz - 1:
        zc = np.append(zc, nz - 1)
    sm = np.empty((len(zc),) + a.shape[1:], dtype=a.dtype)
    for i, z in enumerate(zc):
        sm[i] = np.median(a[max(z - h, 0):min(z + h + 1, nz)], axis=0)
    w = np.interp(np.arange(nz), zc, np.arange(len(zc)))
    lo = np.floor(w).astype(int)
    hi = np.minimum(lo + 1, len(zc) - 1)
    fr = (w - lo).reshape((-1,) + (1,) * (a.ndim - 1))
    return sm[lo] * (1 - fr) + sm[hi] * fr


def _flag_outliers(I, V, good, nsig=5., rel=0.5, window=5, zwindow=3, zmed=15,
                   grow_xy=2, grow_z=2, grow_nsig=2., pad_xy=0, loo=None,
                   zratio=51, zstep=8, fsys=0.05, fgrad=0.5, grow_rel=0.1,
                   smax=5.):
    """
    Flag positive outliers (cosmic rays) in each of N aligned cubes by comparing
    each cube with the other ones on the same grid.

    For cube k, the comparison estimate `est` is the other cube (N=2) or the
    median of the other cubes (N>=3), with variance var_est. The excess is
        D = I_k - est.
    Differences caused by seeing, alignment or flux calibration follow the
    source's (smooth) spectrum, while cosmic rays are only a few pixels wide
    in wavelength. So a running median of D along wavelength (zmed pixels) is
    subtracted first:
        D_hp = D - median_z(D).
    A voxel is flagged only if D_hp passes BOTH tests:
        absolute:  D_hp > nsig * sqrt(V_k + var_est)   (significant compared with noise)
        relative:  D_hp > rel * L                      (large compared with the source)
    where L is the local maximum of |est| in a zwindow x window x window box
    (wavelength, y, x). The relative test protects bright sources and emission
    lines, where small mismatches leave residuals that are significant but
    fractionally small.

    Growth is anisotropic, because cosmic rays are clumped (tracks or blobs)
    spatially but narrow in wavelength (FWHM <~ the instrumental resolution):
      1. spectral: up to grow_z pixels along wavelength,
      2. spatial: up to grow_xy pixels in (x, y), 8-connected so diagonal
         tracks are followed.
    Each step is iterative and conditional: a neighbor is added only if its
    D_hp > grow_nsig * sigma, so the mask follows the connected clump rather
    than a fixed box. pad_xy then adds an unconditional spatial margin
    (e.g. 1 pixel) for the low-level wings spread by resampling.

    With N >= 5 cubes the median of ALL of them is used as the comparison
    estimate (one outlier does not move it), computed once instead of once per
    cube; with fewer, the median of the others is used. `loo` overrides that.


    Three guards keep a bright quasar from being clipped, where frame-to-frame
    seeing and transparency differences are large in absolute terms:
      - `zratio`: each cube is scaled to the comparison estimate per spaxel,
        smoothly along wavelength, before the difference is taken;
      - `fgrad`: a term proportional to the spatial gradient is added to the
        noise, so that a sub-pixel registration error is not significant on
        the steep flanks of a point source;
      - `fsys`: a fractional term is added as well, but see below -- at the
        default thresholds it never binds;
      - `grow_rel`: the relative test is applied while growing too.
    Set zratio=0, fsys=0, fgrad=0, grow_rel=0 for the previous behaviour.

    Args:
        I, V (ndarray [N, nz, ny, nx]): intensity and variance on a common grid.
        zratio (int): channels over which the per-spaxel scale between a cube
            and the comparison estimate is measured; 0 disables the matching.
        zstep (int): channels between evaluations of that scale, which is
            smooth along wavelength. 1 evaluates it everywhere and is the
            slowest part of this function; 8 is ~8x faster, agrees to ~1% and
            largely cancels in the ratio.
        fsys (float): fractional systematic added in quadrature to the noise.
            Inert at the default thresholds, and deliberately kept anyway.
            Because L is the local MAXIMUM of |scale*est|, L >= |scale*est|
            always, so the part of the absolute test this term controls,
                D > nsig * fsys * |scale*est|,
            is weaker than the relative test D > rel * L whenever
                nsig * fsys <= rel,
            which at the defaults is 5 * 0.05 = 0.25 against 0.5. Measured on
            a synthetic bright quasar it exceeded the relative threshold in
            0.000% of the bright voxels (median ratio 0.199, maximum 0.500,
            the bound). The relative test is itself a fractional floor, ten
            times stronger, so this one does nothing while rel is at 0.5.
            It matters when you LOWER rel: rel = 0.5 means a cosmic ray on
            the core of a quasar has to exceed half the quasar to be caught,
            so if you want those you reduce rel, and fsys is then the only
            brightness-aware guard left. Raising fsys above rel / nsig makes
            it active; do not do that without re-checking that real cosmic
            rays are still found.
        fgrad (float): the same for the spatial gradient of the estimate,
            which is what a sub-pixel registration error turns into: a frame
            offset by d pixels differs by about d * |grad I|, so fgrad is that
            offset in pixels. Unlike fsys this is not bounded by L -- on a
            steep flank the per-pixel gradient is a large fraction of the peak
            -- so it does bind, in about 28% of bright voxels, and removes
            roughly a quarter of the remaining false positives.
        grow_rel (float): relative threshold used while growing.
        smax (float): the per-spaxel scale is clipped to [1/smax, smax].
        good (bool ndarray [N, nz, ny, nx]): valid voxels.

    Returns:
        bool ndarray [N, nz, ny, nx]: outlier flags.
    """
    N = I.shape[0]
    Ig = np.where(good, I, np.nan)
    Vg = np.where(good, V, np.nan)
    flags = np.zeros(I.shape, dtype=bool)
    s_xy = np.zeros((3, 3, 3), dtype=bool)
    s_xy[1] = True                                    # 8-connected in (y, x), same slice
    s_z = np.zeros((3, 3, 3), dtype=bool)
    s_z[:, 1, 1] = True                               # along wavelength only
    if loo is None:
        loo = (N < 5)
    Esm_all = None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        if not loo:
            # median of all cubes: robust to one outlier, computed once
            n_all = np.sum(good, axis=0)
            est_all = np.nanmedian(Ig, axis=0)
            var_all = (np.where(n_all >= 3, np.pi / 2., 1.) *
                       np.nanmean(Vg, axis=0) / np.maximum(n_all, 1))
        for k in range(N):
            if not loo:
                est, var_est, n_oth = est_all, var_all, n_all - good[k]
            else:
                oth = [j for j in range(N) if j != k]
                n_oth = np.sum(good[oth], axis=0)
                if len(oth) == 1:
                    est, var_est = Ig[oth[0]], Vg[oth[0]]
                else:
                    est = np.nanmedian(Ig[oth], axis=0)
                    # variance of the median of n values ~ (pi/2) <V>/n; plain mean for n<=2
                    corr = np.where(n_oth >= 3, np.pi / 2., 1.)
                    var_est = corr * np.nanmean(Vg[oth], axis=0) / np.maximum(n_oth, 1)
            testable = good[k] & (n_oth > 0) & np.isfinite(est) & np.isfinite(var_est)
            e0 = np.nan_to_num(est)
            # Scale match before differencing. Seeing and transparency make one
            # frame a smoothly wavelength-dependent multiple of the others at
            # every spaxel, and high-passing that DIFFERENCE turns every narrow
            # spectral feature into a positive spike: where an absorption line
            # cuts the flux, the deficit of a fainter frame shrinks, and the
            # running median leaves a cosmic-ray-shaped residual. Dividing the
            # mismatch out first removes it at the source. The scale is the
            # ratio of the two cubes smoothed over zratio channels (and 3x3
            # spatially, so the faint halo is matched too), which cosmic rays,
            # being a few channels wide, cannot move.
            if zratio and zratio > 1:
                nsm = int(zratio) * 9
                Ism = ndimage.uniform_filter(
                    _crr_zsmooth(np.nan_to_num(Ig[k]), zratio, zstep),
                    size=(1, 3, 3), mode='nearest')
                if loo:
                    Esm = ndimage.uniform_filter(
                        _crr_zsmooth(e0, zratio, zstep), size=(1, 3, 3), mode='nearest')
                else:
                    # est is the same cube for every k, so smooth it once
                    if Esm_all is None:
                        Esm_all = ndimage.uniform_filter(
                            _crr_zsmooth(e0, zratio, zstep), size=(1, 3, 3), mode='nearest')
                    Esm = Esm_all
                ok_s = Esm > 3. * np.sqrt(np.nan_to_num(var_est) / nsm)
                scale = np.where(ok_s & (Esm != 0), Ism / np.where(Esm == 0, 1., Esm), 1.)
                scale = np.clip(np.nan_to_num(scale, nan=1.), 1. / smax, smax)
            else:
                scale = np.ones_like(e0)
            se = scale * e0
            D = np.where(testable, I[k] - se, 0.)
            if zmed and zmed > 1:
                D = D - ndimage.median_filter(D, size=(int(zmed), 1, 1), mode='nearest')
            # Noise budget. The photon term alone scales as sqrt(I), so on a
            # bright quasar any residual percent-level mismatch is formally
            # significant; fgrad and fsys add the systematic floor the
            # frame-to-frame scatter actually has. fgrad is the one that does
            # the work -- fsys is dominated by the relative test below unless
            # nsig * fsys > rel, see the docstring.
            var = V[k] + scale ** 2 * np.nan_to_num(var_est)
            if fsys:
                var = var + (fsys * np.abs(se)) ** 2
            if fgrad:
                var = var + (fgrad * np.maximum(np.abs(np.gradient(se, axis=1)),
                                                np.abs(np.gradient(se, axis=2)))) ** 2
            sig = np.sqrt(np.where(testable, var, np.inf))
            L = ndimage.maximum_filter(np.abs(se),
                                       size=(int(zwindow), int(window), int(window)),
                                       mode='nearest')
            f = testable & (D > nsig * sig) & (D > rel * L)
            if f.any():
                cand = testable & (D > grow_nsig * sig)
                if grow_rel:
                    # the relative test has to hold while growing as well, or a
                    # cosmic ray beside a bright source grows into it, where a
                    # 2-sigma condition is met by the mismatch alone
                    cand &= (D > grow_rel * L)
                for _ in range(int(grow_z)):
                    new = ndimage.binary_dilation(f, structure=s_z) & cand & ~f
                    if not new.any():
                        break
                    f |= new
                for _ in range(int(grow_xy)):
                    new = ndimage.binary_dilation(f, structure=s_xy) & cand & ~f
                    if not new.any():
                        break
                    f |= new
                if pad_xy > 0:
                    f = ndimage.binary_dilation(f, structure=s_xy, iterations=int(pad_xy)) & testable
            flags[k] = f
    return flags


# --------------------------------------------------------------------------
# main routine
# --------------------------------------------------------------------------

def kcwi_combine_stacks(cubefns, outname, wavebin=None, box=None, search_size=10,
                        upfactor=10, align=True, preshift=None, grid='ref',
                        weighting='exposure', fluxscale=None, wave_interp_method='cubic',
                        min_overlap=0.3, nsig_flux=3., crr=False, crr_nsig=5.,
                        crr_rel=0.5, crr_window=5, crr_zwindow=3, crr_zmed=15,
                        crr_grow_xy=2, crr_grow_z=2, crr_grow_nsig=2., crr_pad_xy=0,
                        crr_zratio=51, crr_zstep=8, crr_fsys=0.05, crr_fgrad=0.5,
                        crr_grow_rel=0.1,
                        zblock=200, display=True, overwrite=True):
    """
    Align and co-add stacked KCWI cubes (outputs of kcwi_stack).

    Args:
        cubefns (list of str): stacked intensity cubes, e.g.
            ['q0100_run1_icubes.fits', 'q0100_run2_icubes_wcs.fits'].
            The first one is the reference for both astrometry and flux.
        outname (str): output prefix. Writes <outname>_i<suffix>.fits and the
            matching _v, _m and _e cubes.
        wavebin ([float, float]): wavelength range (Angstrom) for the white-light
            images. Default: WAVGOOD0/1 of the reference header, or the central
            50% of the wavelength range.
        box ([x0, x1, y0, y1]): 0-based pixel box on the *reference* white-light
            image used for the cross-correlation. Default: the whole overlap. Use
            it to restrict the correlation to compact sources.
        search_size (int): coarse search half-width in reference pixels.
        upfactor (int): sub-pixel step = 1/upfactor pixel in the fine pass (+/-1 pix).
        align (bool): False skips the cross-correlation (only preshift is applied),
            e.g. when both cubes are already on absolute astrometry.
        preshift (list of [dRA, dDec] in arcsec, one per cube): an initial offset
            added to CRVAL1/2 before the search (same convention as the
            .preshift.list in kcwi_stack). Use it when offsets exceed search_size.
        grid (str): 'ref' (the reference cube's grid) or 'union' (extended to
            cover all cubes).
        weighting (str): 'exposure' (like kcwi_stack) or 'ivar'.
        fluxscale (None, 'auto', or list of float): multiplicative flux scale per
            cube. 'auto' scales each cube to the reference using the summed flux
            ratio of bright sources in the aligned white-light images. The ratio
            is always printed as a diagnostic.
        wave_interp_method (str): interpolation for the intensity when the
            wavelength grids differ (variance and exposure use 'linear').
        min_overlap (float): minimum fraction of the reference pixels that must
            overlap for a trial shift to count.
        nsig_flux (float): threshold for the pixels used in the flux-ratio estimate.
        crr (bool): reject outliers (cosmic rays) after alignment by comparing each
            cube with the others on the common grid (after flux scaling). A voxel
            is rejected only if its excess over the other cube(s) is both
            significant compared with the noise and large compared with the local
            source level. Only positive outliers are rejected. Voxels covered by
            only one cube cannot be tested.
            The excess has a running median along wavelength (crr_zmed pixels)
            subtracted first, so smooth differences from seeing, alignment or flux
            calibration are not flagged.
        crr_nsig (float): absolute test, excess > crr_nsig * sqrt(V_k + V_others).
        crr_rel (float): relative test, excess > crr_rel * (local max of |others|
            in a crr_zwindow x crr_window x crr_window box). Raise it if real
            differences are flagged on bright sources or lines; lower it to catch
            cosmic rays on top of sources.
        crr_window (int): spatial size (pixels) of the local-level box. It should
            be about the seeing FWHM or larger.
        crr_zwindow (int): spectral size (pixels) of the local-level box.
        crr_zmed (int): width (pixels) of the running median along wavelength that
            is subtracted from the excess. It should be several times wider than a
            cosmic ray; 0 disables it.
        crr_grow_xy (int): maximum spatial growth (pixels, 8-connected) of each
            rejection, iterative and restricted to connected neighbors with
            excess > crr_grow_nsig * sigma. It follows clumped cosmic-ray tracks.
        crr_grow_z (int): maximum spectral growth (pixels), same rule. Keep it
            small; cosmic rays are narrower than the instrumental resolution.
        crr_grow_nsig (float): threshold for growing.
        crr_pad_xy (int): unconditional spatial margin (pixels) added around the
            final mask, for the low-level wings spread by resampling. 0 = off.
        zblock (int): number of wavelength channels processed at a time
            (memory). Each block also reads `pad` channels on either side,
            pad = max(crr_zwindow, crr_zmed, crr_zratio)//2 + crr_grow_z + 1,
            i.e. 28 with the defaults, so a small block spends much of its
            time on overlap: 200 wastes 28%, 600 wastes 9%.
        crr_zstep (int): channels between evaluations of the per-spaxel scale,
            which is smooth along wavelength; 8 is ~8x faster than 1 and gives
            the same flags.
        display (bool): show the diagnostic plots (they are always saved).
        overwrite (bool): overwrite the outputs.

    Returns:
        astropy.table.Table: shifts applied to each cube (cross-correlation
            shift in reference pixels, and the total CRVAL change in arcsec,
            including any preshift), flux ratios, and the number of rejected
            voxels per cube. Also written to <outname>.combine.shift.list.
            With crr=True, a bit mask of rejected voxels (bit k = cube k) is
            written to kcwi_combine/<outname>.crmask.fits.
    """

    if len(cubefns) < 2:
        raise ValueError('Need at least two cubes.')
    if grid not in ['ref', 'union']:
        raise ValueError("grid must be 'ref' or 'union'")
    if weighting not in ['exposure', 'ivar']:
        raise ValueError("weighting must be 'exposure' or 'ivar'")

    if not display:
        oldbackend = matplotlib.get_backend()
        matplotlib.use('Agg')
    outdir = 'kcwi_combine'
    os.makedirs(outdir, exist_ok=True)

    # ---------------- read ----------------
    stacks = []
    for fn in cubefns:
        print('Reading {}'.format(fn))
        stacks.append(_read_stack(fn))
    ref = stacks[0]
    hdr_ref = ref['hdr']
    suffix = ref['fn']['suffix']

    bunits = [s['hdr'].get('BUNIT', '') for s in stacks]
    if len(set(bunits)) > 1:
        raise ValueError('Inconsistent BUNIT between cubes: {}'.format(bunits))
    for s in stacks:
        if s['hdr'].get('CTYPE3', 'WAVE') != hdr_ref.get('CTYPE3', 'WAVE'):
            raise ValueError('Inconsistent CTYPE3 (air vs vacuum?) in {}'.format(s['fn']['i']))

    if wavebin is None:
        if 'WAVGOOD0' in hdr_ref and 'WAVGOOD1' in hdr_ref:
            wavebin = [hdr_ref['WAVGOOD0'], hdr_ref['WAVGOOD1']]
        else:
            w = _wave_axis(hdr_ref)
            wavebin = [np.percentile(w, 25), np.percentile(w, 75)]
    print('White-light range: {:.1f} - {:.1f} A'.format(*wavebin))

    # ---------------- preshift ----------------
    if preshift is not None:
        for s, (dra, ddec) in zip(stacks, preshift):
            s['hdr']['CRVAL1'] += dra / 3600.
            s['hdr']['CRVAL2'] += ddec / 3600.

    w_ref = wcs.WCS(hdr_ref).celestial
    img_ref = _whitelight(ref, wavebin)
    pixscale = np.sqrt(np.abs(np.linalg.det(w_ref.pixel_scale_matrix))) * 3600.

    # ---------------- align ----------------
    shift_px = np.zeros((len(stacks), 2))
    shift_as = np.zeros((len(stacks), 2))
    flux_ratio = np.ones(len(stacks))
    thumbs = [img_ref]

    for k in range(1, len(stacks)):
        s = stacks[k]
        name = os.path.basename(s['fn']['i'])
        print('Aligning {} -> {}'.format(name, os.path.basename(ref['fn']['i'])))
        w_k = wcs.WCS(s['hdr']).celestial
        img_k = _whitelight(s, wavebin)

        d = np.array([0., 0.])
        if align:
            # coarse
            n = int(search_size)
            grid_c = np.arange(-n, n + 1).astype(float)
            crl_c = _xcorr_map(img_ref, w_ref, img_k, w_k, grid_c, grid_c, box, min_overlap)
            if not np.isfinite(crl_c).any():
                raise RuntimeError('No overlap for any trial shift; supply a preshift.')
            ic = np.unravel_index(np.nanargmax(crl_c), crl_c.shape)
            d = np.array([grid_c[ic[0]], grid_c[ic[1]]])
            dc = d.copy()
            if (min(ic) == 0) or (max(ic) == 2 * n):
                print('  [WARNING] correlation peak at the edge of the search area '
                      '({:+.0f}, {:+.0f}) pix. Increase search_size or add a preshift.'.format(*d))

            # fine
            grid_f = np.linspace(-1, 1, 2 * int(upfactor) + 1)
            crl_f = _xcorr_map(img_ref, w_ref, img_k, w_k, grid_f + d[0], grid_f + d[1],
                               box, min_overlap)
            jf = np.unravel_index(np.nanargmax(crl_f), crl_f.shape)
            d = d + np.array([grid_f[jf[0]], grid_f[jf[1]]])
            print('  shift = ({:+.2f}, {:+.2f}) ref pix, peak r = {:.3f}'.format(
                d[0], d[1], np.nanmax(crl_f)))

            # diagnostic plot
            fig, axs = plt.subplots(1, 3, figsize=(12, 4))
            ext_c = [grid_c[0] - 0.5, grid_c[-1] + 0.5, grid_c[0] - 0.5, grid_c[-1] + 0.5]
            axs[0].imshow(crl_c.T, origin='lower', extent=ext_c, cmap='Greys',
                          interpolation='nearest')
            st = grid_f[1] - grid_f[0]
            ext_f = [dc[0] + grid_f[0] - st / 2, dc[0] + grid_f[-1] + st / 2,
                     dc[1] + grid_f[0] - st / 2, dc[1] + grid_f[-1] + st / 2]
            axs[0].imshow(crl_f.T, origin='lower', extent=ext_f, cmap='plasma',
                          interpolation='nearest')
            axs[0].set_xlim(ext_c[:2]); axs[0].set_ylim(ext_c[2:])
            axs[0].plot(d[0], d[1], '+', color='r', ms=10)
            axs[0].set_title('NCC map (ref pix)')
            axs[0].set_xlabel(r'$\Delta x$'); axs[0].set_ylabel(r'$\Delta y$')

        # WCS correction: new CRVAL = w_ref(pix_ref(CRVAL) + d)
        crval = np.array([s['hdr']['CRVAL1'], s['hdr']['CRVAL2']])
        px, py = w_ref.wcs_world2pix(crval[0], crval[1], 0)
        ra_new, dec_new = w_ref.wcs_pix2world(px + d[0], py + d[1], 0)
        dra = (float(ra_new) - crval[0]) * 3600. * np.cos(np.deg2rad(crval[1]))
        ddec = (float(dec_new) - crval[1]) * 3600.
        s['hdr']['CRVAL1'] = float(ra_new)
        s['hdr']['CRVAL2'] = float(dec_new)
        shift_px[k] = d
        pre = np.array(preshift[k], dtype=float) if preshift is not None else np.zeros(2)
        pre[0] *= np.cos(np.deg2rad(crval[1]))
        shift_as[k] = [dra + pre[0], ddec + pre[1]]      # total, including preshift
        print('  CRVAL correction: dRA*cos(dec) = {:+.3f}", dDec = {:+.3f}" '
              '(total incl. preshift: {:+.3f}", {:+.3f}")'.format(dra, ddec, *shift_as[k]))

        # aligned thumbnail and flux ratio
        w_k = wcs.WCS(s['hdr']).celestial
        yi, xi = _shift_map(w_ref, img_ref.shape, w_k)
        thumb = _resample2d(img_k, yi, xi)
        thumbs.append(thumb)
        q = np.isfinite(img_ref) & np.isfinite(thumb)
        if q.sum() > 10:
            # aperture-sum ratio over bright sources (dilated, so it is robust to
            # seeing and resampling differences), with a local background subtracted
            thr0 = np.nanmedian(img_ref[q]) + nsig_flux * np.nanstd(img_ref[q])
            src = ndimage.binary_dilation(np.nan_to_num(img_ref, nan=-np.inf) > thr0, iterations=3)
            qq, qb = q & src, q & ~src
            if qq.sum() > 5 and qb.sum() > 5:
                b0, b1 = np.median(img_ref[qb]), np.median(thumb[qb])
                flux_ratio[k] = np.sum(thumb[qq] - b1) / np.sum(img_ref[qq] - b0)
        print('  flux ratio (this / ref) = {:.3f}'.format(flux_ratio[k]))

        if align:
            v0, v1 = np.nanpercentile(img_ref, [1, 99.5])
            axs[1].imshow(img_ref, origin='lower', cmap='Greys', vmin=v0, vmax=v1)
            lev = np.nanpercentile(thumb, [90, 97, 99.5])
            axs[1].contour(thumb, levels=np.unique(lev), colors='C1', linewidths=0.8)
            if box is not None:
                axs[1].add_patch(plt.Rectangle((box[0] - 0.5, box[2] - 0.5), box[1] - box[0],
                                               box[3] - box[2], fill=False, ec='C0'))
            axs[1].set_title('ref + aligned contours')
            resid = img_ref - thumb / flux_ratio[k]
            r = np.nanpercentile(np.abs(resid), 99)
            axs[2].imshow(resid, origin='lower', cmap='RdBu_r', vmin=-r, vmax=r)
            axs[2].set_title('ref - aligned/ratio')
            fig.suptitle(name)
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, os.path.basename(outname) + '_align_{}.png'.format(k)))
            if display:
                plt.show()
            plt.close(fig)

    # ---------------- flux scale ----------------
    if fluxscale is None:
        scale = np.ones(len(stacks))
    elif isinstance(fluxscale, str) and fluxscale == 'auto':
        scale = 1. / flux_ratio
    else:
        scale = np.asarray(fluxscale, dtype=float)
    print('Flux scale applied: {}'.format(np.round(scale, 4)))

    # ---------------- output grid ----------------
    hdr_out = hdr_ref.copy()
    ny0, nx0 = img_ref.shape
    if grid == 'union':
        xs, ys = [0, nx0 - 1], [0, ny0 - 1]
        for s in stacks[1:]:
            w_k = wcs.WCS(s['hdr']).celestial
            nyk, nxk = s['i'].shape[1:]
            cx = np.array([0, nxk - 1, 0, nxk - 1]); cy = np.array([0, 0, nyk - 1, nyk - 1])
            ra, dec = w_k.wcs_pix2world(cx, cy, 0)
            x, y = w_ref.wcs_world2pix(ra, dec, 0)
            xs += list(x); ys += list(y)
        x0, x1 = int(np.floor(min(xs))), int(np.ceil(max(xs)))
        y0, y1 = int(np.floor(min(ys))), int(np.ceil(max(ys)))
        hdr_out['NAXIS1'] = x1 - x0 + 1
        hdr_out['NAXIS2'] = y1 - y0 + 1
        hdr_out['CRPIX1'] = hdr_ref['CRPIX1'] - x0
        hdr_out['CRPIX2'] = hdr_ref['CRPIX2'] - y0
    w_out = wcs.WCS(hdr_out).celestial
    shape_out = (hdr_out['NAXIS2'], hdr_out['NAXIS1'])
    nz = hdr_ref['NAXIS3']
    print('Output grid: {} x {} x {} ({:.3f}"/pix)'.format(shape_out[1], shape_out[0], nz, pixscale))

    # ---------------- resample and co-add ----------------
    N = len(stacks)
    maps = []
    for k, s in enumerate(stacks):
        if not _same_wave(s['hdr'], hdr_ref):
            print('Resampling wavelength of {} ({})'.format(
                os.path.basename(s['fn']['i']), wave_interp_method))
            stacks[k] = _resample_wave(s, hdr_ref, kind=wave_interp_method)
        maps.append(_shift_map(w_out, shape_out, wcs.WCS(stacks[k]['hdr']).celestial))

    data = np.zeros((nz,) + shape_out)
    vdata = np.zeros_like(data)
    esum = np.zeros_like(data)
    crmask = np.zeros((nz,) + shape_out, dtype=np.int16) if crr else None
    nrej = np.zeros(N, dtype=np.int64)
    ntest = np.zeros(N, dtype=np.int64)
    # margin for the widest wavelength filter, which is the zratio scale
    # window. The truncation at a block edge mostly cancels in the ratio, so
    # this changes very few flags, but it makes the blocked result match the
    # unblocked one exactly -- see the note in kcwi.kcwi_stack.
    pad = (max(int(crr_zwindow), int(crr_zmed), int(crr_zratio)) // 2
           + int(crr_grow_z) + 1) if crr else 0

    print('Resampling and co-adding{}...'.format(' with outlier rejection' if crr else ''))
    for z0 in range(0, nz, zblock):
        z1 = min(z0 + zblock, nz)
        p0, p1 = max(z0 - pad, 0), min(z1 + pad, nz)
        I = np.empty((N, p1 - p0) + shape_out)
        V = np.empty_like(I)
        E = np.empty_like(I)
        for k in range(N):
            I[k], V[k], E[k] = _resample_block(stacks[k], maps[k][0], maps[k][1], p0, p1, shape_out)
            I[k] *= scale[k]
            V[k] *= scale[k]**2
        good = np.isfinite(I) & np.isfinite(V) & (V > 0) & (E > 0)

        if crr:
            fl = _flag_outliers(I, V, good, nsig=crr_nsig, rel=crr_rel, window=crr_window,
                                zwindow=crr_zwindow, zmed=crr_zmed, grow_xy=crr_grow_xy,
                                grow_z=crr_grow_z, grow_nsig=crr_grow_nsig, pad_xy=crr_pad_xy,
                                zratio=crr_zratio, zstep=crr_zstep, fsys=crr_fsys,
                                fgrad=crr_fgrad, grow_rel=crr_grow_rel)
            c = slice(z0 - p0, z0 - p0 + (z1 - z0))
            fl_c = fl[:, c]
            nrej += fl_c.reshape(N, -1).sum(axis=1)
            multi = good[:, c].sum(axis=0) > 1
            ntest += (good[:, c] & multi[None]).reshape(N, -1).sum(axis=1)
            for k in range(N):
                crmask[z0:z1] |= (fl_c[k].astype(np.int16) << k)
            good &= ~fl

        c = slice(z0 - p0, z0 - p0 + (z1 - z0))
        I, V, E, good = I[:, c], V[:, c], E[:, c], good[:, c]
        if weighting == 'exposure':
            wt = np.where(good, E / (scale**2)[:, None, None, None], 0.)   # as kcwi_stack
        else:
            wt = np.where(good, 1. / np.where(good, V, 1.), 0.)
        wsum = wt.sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            data[z0:z1] = np.where(wsum > 0, np.sum(wt * np.where(good, I, 0.), axis=0) / wsum, 0.)
            vdata[z0:z1] = np.where(wsum > 0, np.sum(wt**2 * np.where(good, V, 0.), axis=0) / wsum**2, 0.)
        esum[z0:z1] = np.sum(np.where(good, E, 0.), axis=0)
    del stacks
    mdata = (esum == 0).astype(np.int16)

    if crr:
        for k in range(N):
            frac = nrej[k] / max(ntest[k], 1)
            print('  {}: rejected {} voxels ({:.3g}% of overlapping voxels)'.format(
                os.path.basename(cubefns[k]), nrej[k], 100 * frac))
            if frac > 0.01:
                print('  [WARNING] >1% rejected. Check the alignment/flux scale, or raise '
                      'crr_rel / crr_nsig.')

    # ---------------- write ----------------
    hdr_out['HISTORY'] = 'kcwi_combine_stacks: co-added {} stacked cubes'.format(len(cubefns))
    for k, fn in enumerate(cubefns):
        hdr_out['COMB{}'.format(k)] = (os.path.basename(fn)[:60], 'input stack')
        hdr_out['COMBSX{}'.format(k)] = (round(shift_px[k, 0], 3), 'x shift [ref pix]')
        hdr_out['COMBSY{}'.format(k)] = (round(shift_px[k, 1], 3), 'y shift [ref pix]')
        hdr_out['COMBFS{}'.format(k)] = (round(float(scale[k]), 5), 'flux scale')
    hdr_out['COMBWT'] = (weighting, 'combination weighting')
    if crr:
        hdr_out['COMBCRR'] = (True, 'outlier rejection applied')
        hdr_out['CRRNSIG'] = (crr_nsig, 'absolute threshold [sigma]')
        hdr_out['CRRREL'] = (crr_rel, 'relative threshold [local level]')
        hdr_out['CRRWIN'] = (crr_window, 'local-level window, spatial [pix]')
        hdr_out['CRRZWIN'] = (crr_zwindow, 'local-level window, spectral [pix]')
        hdr_out['CRRZMED'] = (crr_zmed, 'spectral median high-pass [pix]')
        hdr_out['CRRGRXY'] = (crr_grow_xy, 'max spatial growth [pix]')
        hdr_out['CRRGRZ'] = (crr_grow_z, 'max spectral growth [pix]')
        hdr_out['CRRGRSIG'] = (crr_grow_nsig, 'growth threshold [sigma]')
        hdr_out['CRRPADXY'] = (crr_pad_xy, 'unconditional spatial pad [pix]')
        for k in range(len(cubefns)):
            hdr_out['CRRN{}'.format(k)] = (int(nrej[k]), 'rejected voxels in cube {}'.format(k))

    vhdr = hdr_out.copy()
    b = hdr_out.get('BUNIT', '')
    if b.startswith('1e-16 erg') or b.startswith('10^(-8)erg'):
        # '10^(-8)erg/s/cm3/arcsec2' is the pre-Nov-2026 kcwi_stack string
        vhdr['BUNIT'] = '1e-32 erg2/s2/cm4/Angstrom2/arcsec4'
    elif b.startswith('count'):
        vhdr['BUNIT'] = 'count2/s2/arcsec4'
    mhdr = hdr_out.copy()
    if 'BUNIT' in mhdr:
        del mhdr['BUNIT']
    ehdr = hdr_out.copy()
    ehdr['BUNIT'] = 's'

    fits.PrimaryHDU(data, header=hdr_out).writeto(outname + '_i' + suffix + '.fits', overwrite=overwrite)
    fits.PrimaryHDU(vdata, header=vhdr).writeto(outname + '_v' + suffix + '.fits', overwrite=overwrite)
    fits.PrimaryHDU(mdata, header=mhdr).writeto(outname + '_m' + suffix + '.fits', overwrite=overwrite)
    fits.PrimaryHDU(esum, header=ehdr).writeto(outname + '_e' + suffix + '.fits', overwrite=overwrite)

    if crr:
        chdr = hdr_out.copy()
        if 'BUNIT' in chdr:
            del chdr['BUNIT']
        chdr['COMMENT'] = 'bit k set = voxel rejected in input cube k (COMBk)'
        fits.PrimaryHDU(crmask, header=chdr).writeto(
            os.path.join(outdir, os.path.basename(outname) + '.crmask.fits'), overwrite=True)

    # aligned white-light thumbnails on the reference grid
    hdr2d = w_ref.to_header()
    fits.PrimaryHDU(np.array(thumbs), header=hdr2d).writeto(
        os.path.join(outdir, os.path.basename(outname) + '.thum.fits'), overwrite=True)

    tab = table.Table([[os.path.basename(f) for f in cubefns],
                       np.round(shift_px[:, 0], 3), np.round(shift_px[:, 1], 3),
                       np.round(shift_as[:, 0], 4), np.round(shift_as[:, 1], 4),
                       np.round(flux_ratio, 4), np.round(scale, 4), nrej],
                      names=['file', 'dx_pix', 'dy_pix', 'dRAcosd_arcsec', 'dDec_arcsec',
                             'flux_ratio', 'flux_scale', 'n_rejected'])
    ascii.write(tab, outname + '.combine.shift.list', overwrite=True, format='fixed_width')
    print('Wrote {}_[i,v,m,e]{}.fits'.format(outname, suffix))

    if not display:
        matplotlib.use(oldbackend)
    return tab
