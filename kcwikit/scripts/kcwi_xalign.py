#!/usr/bin/env python
"""
kcwi_xalign -- align every KCWI blue AND red frame to one common sky frame by
cross-correlation, without a point source, and write the .shift.list files
for kcwi_stack so that the blue and red stacks land on the same sky grid.

Why. kcwi_align registers each channel to its own first frame, and
kcwi_astrometry then needs a point source (ref_xy / ref_ad) in each stacked
channel. Without one, the blue and red stacks can end up offset from each
other by a fraction of a pixel -- which is enough to make a step at the
dichroic in a bright, extended source -- and the absolute WCS is whatever the
headers say.

What it does (step `shifts`)
  1. A white-light image of every frame, in a window next to the dichroic in
     each channel (where the continuum colours of blue and red are closest),
     with emission and sky lines cut out (--z), projected on one fine sky grid.
  2. Within each channel, every frame is registered to the leave-one-out
     template of the others (the template never contains the frame being
     measured, which would pull its shift towards zero), iterated to
     convergence. Mean shift of each channel = 0 (the mean header position).
  3. The red template is registered to the blue template: one blue-red offset.
  4. Optional (--ref image.fits): the combined template is registered to an
     external image (HST, Legacy Survey, Pan-STARRS, ...) for the absolute WCS.
     Without --ref the absolute position is the mean blue header position.
  The match is a least-squares fit  frame = a * template(x - d) + b + c x + e y
  over the overlap, so transparency, background levels and background
  gradients (e.g. sky-subtraction differences between the channels) are
  fitted out and do not move the shift.
  5. Writes <blue>.shift.list and <red>.shift.list (old ones backed up), a
     JSON solution and a diagnostic PNG in kcwi_xalign/.

Then stack both channels WITHOUT use_astrom (the absolute offset is already in
the shifts) and on the same output grid (same stack_ad / dimension / xpix /
ypix / orientation in both .par files; the script checks), and run step `wcs`:
it checks the two stacks against each other (and --ref), and writes
<stack>_icubes_wcs.fits.

    kcwi_xalign shifts IRAS03582+6012.list IRAS03582+6012-red.list --z 0.030 [--ref hst.fits]
    python -c "from kcwikit.kcwi.kcwi import kcwi_stack; kcwi_stack('IRAS03582+6012.list'); kcwi_stack('IRAS03582+6012-red.list')"
    kcwi_xalign wcs IRAS03582+6012.list IRAS03582+6012-red.list --z 0.030 [--ref hst.fits]

A source name instead of the two lists: `kcwi_xalign wcs IRAS03582+6012 --suffix-blue _skl`
reads IRAS03582+6012_skl.list and IRAS03582+6012-red.list. Many sources:
`kcwi_xalign batch SOURCES --steps shifts,stack,wcs` (see `kcwi_xalign batch -h`).

Run in the directory holding the .list / .par files (as for kcwi_align).
"""
import argparse
import json
import os
import shutil
import sys
import time
import warnings

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from astropy import wcs as awcs
from scipy import ndimage

try:
    from reproject import reproject_interp
except ImportError:                                    # pragma: no cover
    sys.exit('kcwi_xalign needs the reproject package (pip install reproject)')

# Rest-frame lines cut out of the white-light windows when --z is given [A]
LINES = [3727.4, 3868.8, 4101.7, 4340.5, 4861.3, 4958.9, 5006.8, 5875.6, 6300.3, 6363.8,
         6548.0, 6562.8, 6583.5, 6716.4, 6730.8, 7135.8, 9068.6, 9530.6]
SKY = [5577.34, 5889.95, 5895.92, 6300.30, 6363.78]   # observed [A]


# ===========================================================================
# frames
# ===========================================================================
def _readpar(parfn):
    try:
        from kcwikit.kcwi.kcwi import kcwi_stack_readpar
        return kcwi_stack_readpar(parfn)
    except Exception:
        par = {'stack_ad': np.array([-1., -1.]), 'stack_dimension': np.array([-1, -1]),
               'stack_xpix': -1., 'stack_ypix': -1., 'stack_orientation': -1000.}
        alias = {'ad': 'stack_ad', 'dimension': 'stack_dimension', 'xpix': 'stack_xpix',
                 'ypix': 'stack_ypix', 'orientation': 'stack_orientation'}
        if os.path.isfile(parfn):
            for line in open(parfn):
                s = line.split()
                if len(s) < 2:
                    continue
                k = alias.get(s[0], s[0])
                if k in ('stack_ad', 'stack_dimension'):
                    par[k] = np.array(s[1:3], dtype=float)
                elif k in par:
                    par[k] = float(s[1])
        return par


def read_list(fnlist, cubed=False):
    """Frames of a kcwi_stack .list: names (as in .shift.list), cube files, trims, preshifts."""
    suffix = 'cubed' if cubed else 'cubes'
    tab = ascii.read(fnlist, format='no_header')
    base = list(tab['col1'])
    fn = [b + '_i' + suffix + '.fits' for b in base]
    trim = np.array([tab['col2'], tab['col3']]).T.astype(int)
    pre = {}
    for cand in (fnlist.replace('.list', '.preshift.list'), fnlist.replace('.list', '.pre.list')):
        if os.path.isfile(cand):
            pt = ascii.read(cand, format='no_header')
            pre = {p + '_i' + suffix + '.fits': (float(r), float(d))
                   for p, r, d in zip(pt['col1'], pt['col2'], pt['col3'])}
            break
    pres = [pre.get(os.path.basename(f), (0., 0.)) for f in fn]
    return base, fn, trim, pres


def _wave(h):
    step = h['CD3_3'] if 'CD3_3' in h else h.get('CDELT3', 1.) * h.get('PC3_3', 1.)
    w = h['CRVAL3'] + (np.arange(h['NAXIS3']) + 1 - h['CRPIX3']) * step
    if str(h.get('CUNIT3', 'Angstrom')).strip().lower() in ('m', 'meter', 'metre'):
        w = w * 1e10
    return w


def _window(h, band, z, side):
    """Wavelength selection for the white-light image."""
    w = _wave(h)
    g0, g1 = h.get('WAVGOOD0', w[0]), h.get('WAVGOOD1', w[-1])
    if band is None:            # next to the dichroic
        band = (g1 - 420., g1 - 20.) if side == 'blue' else (g0 + 20., g0 + 420.)
    q = (w > max(band[0], g0)) & (w < min(band[1], g1))
    cut = [s for s in SKY]
    if z is not None:
        cut += [l * (1 + z) for l in LINES]
    for c in cut:
        q &= np.abs(w - c) > max(8., c * 1000. / 299792.458)      # +-1000 km/s
    return q, band


def frame_image(fn, band, z, side, trim=(0, 0), npix_trim=1, preshift=(0., 0.), use_regmask=True):
    """White-light image (ny, nx) of one cube, and its 2-D WCS with the preshift applied."""
    h = fits.getheader(fn)
    q, band = _window(h, band, z, side)
    if q.sum() < 20:
        raise ValueError('%s: only %d channels in the white-light window %s' % (fn, q.sum(), band))
    k = np.where(q)[0]
    k0, k1 = k[0], k[-1] + 1
    qq = q[k0:k1]
    # read only the planes of the window (stacked cubes are large)
    with fits.open(fn, memmap=False) as hl:
        cube = np.array(hl[0].section[k0:k1], dtype=float)[qq]
        flags = np.array(hl['FLAGS'].section[k0:k1])[qq] if (len(hl) >= 4 and 'FLAGS' in hl) else None
    bad = ~np.isfinite(cube) | (cube == 0)
    mfn = fn.replace('_icubes', '_mcubes').replace('_icubed', '_mcubed')
    if mfn != fn and os.path.isfile(mfn):
        with fits.open(mfn, memmap=False) as hm:
            bad |= np.array(hm[0].section[k0:k1])[qq] != 0
    c = np.where(bad, np.nan, cube)
    del cube, bad
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        img = np.nanmedian(c, axis=0)                     # robust to cosmic rays
    if flags is not None:                                  # pyDRP: whole-spaxel flags
        f = np.array(flags)
        n = npix_trim
        f[:, :n, :] = f[:, -n:, :] = f[:, :, :n] = f[:, :, -n:] = 128
        img[np.nanmedian(f, axis=0) == 128] = np.nan
    # kcwi_align-style trimming (rows = slices at the list's trims, plus a margin)
    good = np.isfinite(img)
    if good.any():
        yy, xx = np.where(good)
        img[:yy.min() + trim[0], :] = np.nan
        img[yy.max() - trim[1] + 1:, :] = np.nan
        img[:, :xx.min() + npix_trim] = np.nan
        img[:, xx.max() - npix_trim + 1:] = np.nan
    w2 = awcs.WCS(h).celestial
    if use_regmask:
        reg = fn.replace('.fits', '_rm.reg')
        if os.path.isfile(reg):
            # a region file that exists but cannot be applied is an error, not
            # a warning: silently skipping it changes the alignment result
            from kcwikit.kcwi.kcwi_regions import regions_to_mask
            hh = w2.to_header()
            hh['NAXIS'], hh['NAXIS1'], hh['NAXIS2'] = 2, img.shape[1], img.shape[0]
            m = regions_to_mask(reg, header=hh, shape=img.shape)
            img[m] = np.nan
            print('  region mask %s: %d pixel(s) masked' % (os.path.basename(reg), int(m.sum())))
    w2.wcs.crval = w2.wcs.crval + np.array(preshift) / 3600.
    return img, w2, band


# ===========================================================================
# common grid and image helpers
# ===========================================================================
def make_grid(wcss, shapes, center, pix_arcsec, margin_arcsec):
    w = awcs.WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = center
    w.wcs.cdelt = [-pix_arcsec / 3600., pix_arcsec / 3600.]
    w.wcs.crpix = [1., 1.]
    xs, ys = [], []
    for wf, sh in zip(wcss, shapes):
        ra, de = wf.all_pix2world([-0.5, sh[1] - 0.5, sh[1] - 0.5, -0.5], [-0.5, -0.5, sh[0] - 0.5, sh[0] - 0.5], 0)
        x, y = w.all_world2pix(ra, de, 0)
        xs += list(x); ys += list(y)
    m = margin_arcsec / pix_arcsec
    x0, x1, y0, y1 = min(xs) - m, max(xs) + m, min(ys) - m, max(ys) + m
    w.wcs.crpix = [1. - x0, 1. - y0]
    shape = (int(np.ceil(y1 - y0)) + 1, int(np.ceil(x1 - x0)) + 1)
    return w, shape


def project(img, wf, wg, shape):
    out, fp = reproject_interp((img, wf), wg, shape_out=shape, order='bilinear')
    out[~(fp > 0.999)] = np.nan
    return out


def nan_smooth(a, sigma_pix):
    if not sigma_pix:
        return a
    v = np.isfinite(a)
    num = ndimage.gaussian_filter(np.where(v, a, 0.), sigma_pix)
    den = ndimage.gaussian_filter(v.astype(float), sigma_pix)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = num / den
    out[~v] = np.nan
    return out


def highpass(a, sigma_pix):
    """Remove structure broader than sigma_pix (frame-to-frame sky-level patterns)."""
    if not sigma_pix:
        return a
    return a - nan_smooth(a, sigma_pix)


def nan_shift(a, d, order=3):
    """a shifted by d = (dy, dx): out(x) = a(x - d); invalid where any input was."""
    v = np.isfinite(a)
    out = ndimage.shift(np.where(v, a, 0.), d, order=order, mode='constant', cval=0.)
    vv = ndimage.shift(v.astype(float), d, order=1, mode='constant', cval=0.) > 0.999
    vv = ndimage.binary_erosion(vv, iterations=2 if order > 1 else 1)
    out[~vv] = np.nan
    return out


def _int_shift(a, dy, dx):
    out = np.full_like(a, np.nan)
    ny, nx = a.shape
    ys, yd = (slice(0, ny - dy), slice(dy, ny)) if dy >= 0 else (slice(-dy, ny), slice(0, ny + dy))
    xs, xd = (slice(0, nx - dx), slice(dx, nx)) if dx >= 0 else (slice(-dx, nx), slice(0, nx + dx))
    out[yd, xd] = a[ys, xs]
    return out


# ===========================================================================
# registration
# ===========================================================================
def _block(a, f):
    """Block-average by f, NaN-aware."""
    ny, nx = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    b = a[:ny, :nx].reshape(ny // f, f, nx // f, f)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(b, axis=(1, 3))


def _coarse(I, R, g0, r, minfrac):
    """
    Integer shift within r of g0 maximising the correlation coefficient.

    A candidate needs an overlap of at least minfrac x the overlap at g0 (and 50
    pixels), so the peak cannot run off to a shift where only a sliver of the
    images overlaps. Relative to the overlap at g0, not to the image's own size:
    two frames at different position angles overlap by much less than either.
    """
    vI = np.isfinite(I)
    n0 = int((vI & np.isfinite(_int_shift(R, int(g0[0]), int(g0[1])))).sum())
    nmin = max(50, minfrac * n0)
    best = (-np.inf, int(g0[0]), int(g0[1]))
    for dy in range(int(g0[0]) - r, int(g0[0]) + r + 1):
        for dx in range(int(g0[1]) - r, int(g0[1]) + r + 1):
            if (dy - g0[0]) ** 2 + (dx - g0[1]) ** 2 > r * r:
                continue
            Rs = _int_shift(R, dy, dx)
            m = vI & np.isfinite(Rs)
            if m.sum() < nmin:
                continue
            a, b = I[m], Rs[m]
            a = a - a.mean(); b = b - b.mean()
            den = np.sqrt((a * a).sum() * (b * b).sum())
            rho = (a * b).sum() / den if den > 0 else -np.inf
            if rho > best[0]:
                best = (rho, dy, dx)
    return best


def register(img, ref, guess=(0., 0.), radius=20, gradient=True, minfrac=0.6, err_scale=1.):
    """
    Displacement d = (dy, dx) [grid pixels] with img(x) ~ a ref(x - d) + b (+ c x + e y).

    Coarse: integer shifts within `radius` of `guess`, maximising the correlation
    coefficient over the overlap (>= minfrac of img's valid pixels). Fine: the
    residual of the linear model, minimised over continuous d on a FIXED set of
    pixels (so the objective is smooth). Error: from the curvature of the
    residual, times err_scale (pixels of the fine grid are correlated).
    """
    vi = np.isfinite(img)
    if vi.sum() < 50:
        raise ValueError('too few valid pixels in the image to register')
    yy, xx = np.where(vi)
    pad = int(radius) + int(np.ceil(np.abs(guess).max())) + 4
    sl = (slice(max(yy.min() - pad, 0), yy.max() + pad + 1), slice(max(xx.min() - pad, 0), xx.max() + pad + 1))
    I, R = img[sl], ref[sl]
    vI = np.isfinite(I)
    g0 = np.round(guess).astype(int)
    r = int(radius)
    if r > 15:
        # large search: first on a block-averaged copy, then +-2 blocks at full resolution
        f = int(np.ceil(r / 12.))
        gd = _coarse(_block(I, f), _block(R, f), np.round(g0 / float(f)).astype(int),
                     int(np.ceil(r / float(f))) + 1, minfrac)[1:]
        g0, r = np.array(gd) * f, 2 * f + 1
    rho_c, dy0, dx0 = _coarse(I, R, g0, r, minfrac)
    if not np.isfinite(rho_c):
        raise ValueError('no valid overlap within the search radius')
    if max(abs(dy0 - g0[0]), abs(dx0 - g0[1])) >= r - 1:
        print('  [WARNING] the correlation peak is at the edge of the search radius; raise --radius')

    # fixed pixel set for the fine fit
    m0 = vI & np.isfinite(nan_shift(R, (dy0, dx0)))
    m0 = ndimage.binary_erosion(m0, iterations=2)
    gy, gx = np.indices(I.shape)
    gy = (gy - gy.mean()) / I.shape[0]; gx = (gx - gx.mean()) / I.shape[1]
    y = I[m0]
    C = ndimage.spline_filter(np.where(np.isfinite(R), R, 0.), order=3)   # prefilter once

    def model_cols(d):
        Rs = ndimage.shift(C, d, order=3, mode='constant', cval=0., prefilter=False)[m0]
        cols = [Rs, np.ones_like(Rs)]
        if gradient:
            cols += [gx[m0], gy[m0]]
        return np.vstack(cols).T

    def S(d):
        X = model_cols(d)
        th, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ th
        return float(res @ res)

    def quad_step(c0, h):
        """S on a 3x3 (h small) or 5x5 grid around c0, 2-D quadratic fit: minimum and Hessian."""
        k = np.arange(-2, 3) if h >= 0.5 else np.arange(-1, 2)
        P = np.array([(c0[0] + i * h, c0[1] + j * h) for i in k for j in k])
        v = np.array([S(p) for p in P])
        u, w = P[:, 0] - c0[0], P[:, 1] - c0[1]
        A = np.vstack([np.ones_like(u), u, w, u * u, u * w, w * w]).T
        q, *_ = np.linalg.lstsq(A, v, rcond=None)
        H = np.array([[2 * q[3], q[4]], [q[4], 2 * q[5]]])
        try:
            st = -np.linalg.solve(H, q[1:3])
        except np.linalg.LinAlgError:
            st = np.zeros(2)
        if not np.all(np.isfinite(st)) or np.any(np.linalg.eigvalsh(H) <= 0):
            st = P[np.argmin(v)] - c0                       # not convex: best grid point
        lim = 2 * h * k.max()
        return c0 + np.clip(st, -lim, lim), H

    d, _ = quad_step(np.array([dy0, dx0], float), 0.5)      # +-1 px around the integer peak
    d, H = quad_step(d, 0.15)                                # refine
    d, H = quad_step(d, 0.05)
    X = model_cols(d)
    th, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ th
    n, p = len(y), X.shape[1] + 2
    s2 = float(res @ res) / max(n - p, 1)
    try:
        cov = 2. * s2 * np.linalg.inv(H) * err_scale ** 2
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.array([np.nan, np.nan])
    yv = y - y.mean()
    rho = 1. - float(res @ res) / float(yv @ yv) if (yv @ yv) > 0 else np.nan
    return dict(d=d, err=err, a=th[0], b=th[1], rho_coarse=float(rho_c), r2=float(rho), npix=int(n))


def align_channel(imgs, radius, gradient, err_scale, nit=8, tol=0.003, label=''):
    """Leave-one-out iterative registration of the frames of one channel."""
    N = len(imgs)
    d = np.zeros((N, 2)); a = np.ones(N); b = np.zeros(N)
    info = [None] * N
    rad_next = 3
    if N == 1:
        return d, np.zeros((N, 2)), [dict(d=np.zeros(2), err=np.zeros(2), a=1., b=0., r2=np.nan, npix=0)]
    for it in range(nit):
        aligned = [(nan_shift(imgs[j], -d[j]) - b[j]) / a[j] for j in range(N)]
        new = np.zeros_like(d)
        for k in range(N):
            stack = np.array([aligned[j] for j in range(N) if j != k])
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                tpl = np.nanmean(stack, axis=0)
                tpl[np.isfinite(stack).sum(axis=0) < max(1, (N - 1) // 2)] = np.nan
            r = register(imgs[k], tpl, guess=d[k], radius=radius if it == 0 else rad_next,
                         gradient=gradient, err_scale=err_scale)
            # The template of the others is itself off by -e_k/(N-1) (e = truth -
            # current estimate, sum e = 0), so the measured m_k = d_k + e_k/(N-1)
            # and the exact (linear) update is d_k = ((N-1) m_k + d_k_old) / N.
            new[k] = ((N - 1) * r['d'] + d[k]) / N
            a[k], b[k], info[k] = r['a'], r['b'], r
        new -= new.mean(axis=0)                         # gauge: mean header position
        rad_next = int(np.ceil(2 * np.abs(new - d).max())) + 3
        step = np.abs(new - d).max()
        d = new
        print('  %s iteration %d: largest change %.4f px' % (label, it + 1, step))
        if step < tol:
            break
    err = np.array([r['err'] for r in info])
    return d, err, info


def template(imgs, d, info):
    al = [(nan_shift(im, -dd) - r['b']) / r['a'] for im, dd, r in zip(imgs, d, info)]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(np.array(al), axis=0)


# ===========================================================================
# output
# ===========================================================================
def crpix_shift(wf_pre, D, wg):
    """CRPIX change for a frame whose image sits at +D (grid px) from the truth."""
    c = wf_pre.wcs.crval
    gx, gy = wg.all_world2pix(c[0], c[1], 0)
    ra0, de0 = wg.all_pix2world(gx, gy, 0)
    ra1, de1 = wg.all_pix2world(gx - D[1], gy - D[0], 0)
    dra, dde = float(ra1 - ra0), float(de1 - de0)     # correction to add to header coords
    x, y = wf_pre.all_world2pix(c[0] - dra, c[1] - dde, 1)
    return float(x - wf_pre.wcs.crpix[0]), float(y - wf_pre.wcs.crpix[1]), dra, dde


def _backup(fn, stamp):
    if os.path.isfile(fn):
        shutil.copy2(fn, fn + '.bak_' + stamp)


def check_par(blist, rlist):
    pb, pr = _readpar(blist.replace('.list', '.par')), _readpar(rlist.replace('.list', '.par'))
    keys = ['stack_ad', 'stack_dimension', 'stack_xpix', 'stack_ypix', 'stack_orientation']
    diff = [k for k in keys if not np.allclose(np.atleast_1d(pb[k]), np.atleast_1d(pr[k]))]
    if diff:
        print('[WARNING] the blue and red .par files define different stack grids (%s).' % ', '.join(diff))
        print('          GIST combines the two cubes spaxel by spaxel: give both the same values.')
    if np.atleast_1d(pb['stack_ad'])[0] == -1 or np.atleast_1d(pr['stack_ad'])[0] == -1:
        print('[NOTE] no stack_ad in the .par: each stack is centred on its own first frame. Add the')
        print('       same "stack_ad RA DEC" line to both .par files to pin the grid.')
    return pb, pr


def plot_summary(fn, res, TB, TR, tr_info, REF=None, ref_info=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ncol = 4 if REF is not None else 3
    fig, ax = plt.subplots(2, ncol, figsize=(4.2 * ncol, 8.2))
    v = np.isfinite(TB) | np.isfinite(TR)
    if v.any():
        yy, xx = np.where(v)
        crop = (slice(max(yy.min() - 15, 0), yy.max() + 16), slice(max(xx.min() - 15, 0), xx.max() + 16))
        TB, TR = TB[crop], TR[crop]
        REF = REF[crop] if REF is not None else None
    else:
        crop = (slice(None), slice(None))
    def show(a_, im, title):
        f = im[np.isfinite(im)]
        lo, hi = (np.percentile(f, [1, 99.5]) if f.size else (0, 1))
        a_.imshow(im, origin='lower', vmin=lo, vmax=hi, cmap='Greys_r')
        a_.set_title(title, fontsize=9); a_.set_xticks([]); a_.set_yticks([])
    show(ax[0, 0], TB, 'blue template')
    show(ax[0, 1], TR, 'red template (header position)')
    model = tr_info['a'] * nan_shift(TB, tr_info['d']) + tr_info['b']      # (gradient terms left out)
    resid = TR - model
    f = resid[np.isfinite(resid)]
    s = np.percentile(np.abs(f), 99) if f.size else 1
    ax[0, 2].imshow(resid, origin='lower', vmin=-s, vmax=s, cmap='RdBu_r')
    ax[0, 2].set_title('red - a*blue(aligned) - b', fontsize=9); ax[0, 2].set_xticks([]); ax[0, 2].set_yticks([])
    if REF is not None:
        show(ax[0, 3], REF, 'reference image (on the grid)')
    pix = res['grid_pix_arcsec']
    for k, (ch, col) in enumerate((('blue', 'C0'), ('red', 'C3'))):
        fr = res['frames'][ch]
        D = np.array([f_['D_grid'] for f_ in fr]) * pix; E = np.array([f_['err_grid'] for f_ in fr]) * pix
        a_ = ax[1, k]
        a_.errorbar(-D[:, 1], D[:, 0], xerr=E[:, 1], yerr=E[:, 0], fmt='o', color=col, ms=4)
        for i, f_ in enumerate(fr):
            a_.annotate(str(i), (-D[i, 1], D[i, 0]), fontsize=7, xytext=(3, 3), textcoords='offset points')
        a_.axhline(0, color='0.7', lw=0.5); a_.axvline(0, color='0.7', lw=0.5)
        a_.set_xlabel('header error, E [arcsec]'); a_.set_ylabel('header error, N [arcsec]')
        a_.set_title('%s frames (total, incl. blue-red and absolute)' % ch, fontsize=9)
        a_.invert_xaxis()
    a_ = ax[1, 2]; a_.axis('off')
    lines = ['grid %.3f"/px' % pix,
             'red vs blue: dE %+.3f"  dN %+.3f"  (+- %.3f, %.3f)' % (
                 -tr_info['d'][1] * pix, tr_info['d'][0] * pix, tr_info['err'][1] * pix, tr_info['err'][0] * pix),
             '   R^2 of the match %.4f' % tr_info['r2']]
    if ref_info is not None:
        lines += ['absolute (vs --ref): dE %+.3f"  dN %+.3f"' % (-ref_info['d'][1] * pix, ref_info['d'][0] * pix),
                  '   R^2 %.4f' % ref_info['r2']]
    else:
        lines += ['absolute: mean blue header position (no --ref)']
    lines += ['', 'windows: blue %.0f-%.0f A, red %.0f-%.0f A' % tuple(res['band_blue'] + res['band_red'])]
    a_.text(0, 1, '\n'.join(lines), va='top', family='monospace', fontsize=8)
    if REF is not None:
        ax[1, 3].axis('off')
    fig.tight_layout()
    fig.savefig(fn, dpi=110)
    plt.close(fig)


# ===========================================================================
# steps
# ===========================================================================
def _ref_image(fn, hdu, wg, shape, smooth_pix):
    with fits.open(fn) as hl:
        h = hl[hdu].header
        data = np.array(hl[hdu].data, dtype=float)
    w = awcs.WCS(h).celestial
    while data.ndim > 2:
        data = data[0]
    out, fp = reproject_interp((data, w), wg, shape_out=shape, order='bilinear')
    out[~(fp > 0.999)] = np.nan
    return nan_smooth(out, smooth_pix)


def step_shifts(a):
    stamp = time.strftime('%Y%m%d_%H%M%S')
    os.makedirs('kcwi_xalign', exist_ok=True)
    check_par(a.blue, a.red)
    chans = {}
    for ch, lst, band in (('blue', a.blue, a.wb), ('red', a.red, a.wr)):
        base, fn, trim, pres = read_list(lst, a.cubed)
        imgs, wcss, used = [], [], None
        for f, t, p in zip(fn, trim, pres):
            im, w2, used = frame_image(f, band, a.z, ch, trim=t, npix_trim=a.npix_trim,
                                       preshift=p, use_regmask=not a.no_regmask)
            imgs.append(im); wcss.append(w2)
        chans[ch] = dict(list=lst, base=base, fn=fn, imgs=imgs, wcs=wcss, band=[float(used[0]), float(used[1])])
        print('%s: %d frames, white light %.0f-%.0f A' % (ch, len(fn), used[0], used[1]))

    pb = _readpar(a.blue.replace('.list', '.par'))
    ad = np.atleast_1d(pb['stack_ad'])
    center = ad if ad[0] != -1 else chans['blue']['wcs'][0].wcs.crval
    allw = chans['blue']['wcs'] + chans['red']['wcs']
    alls = [im.shape for im in chans['blue']['imgs'] + chans['red']['imgs']]
    wg, shape = make_grid(allw, alls, center, a.gridpix, max(a.radius, a.ref_radius if a.ref else 0) + 1.)
    print('common grid: %d x %d at %.3f"/px' % (shape[1], shape[0], a.gridpix))

    # native pixel area / grid pixel area -> the grid pixels are correlated
    w0 = chans['blue']['wcs'][0]
    native = np.abs(np.linalg.det(w0.pixel_scale_matrix)) * 3600. ** 2
    err_scale = np.sqrt(max(native / a.gridpix ** 2, 1.))
    sm = a.smooth / a.gridpix / 2.3548 if a.smooth else 0
    rad = int(np.ceil(a.radius / a.gridpix))
    hp = a.highpass / a.gridpix / 2.3548 if a.highpass else 0

    for ch in ('blue', 'red'):
        c = chans[ch]
        c['g'] = [highpass(nan_smooth(project(im, w, wg, shape), sm), hp)
                  for im, w in zip(c['imgs'], c['wcs'])]
        print('aligning the %s frames' % ch)
        c['d'], c['err'], c['info'] = align_channel(c['g'], rad, not a.no_gradient, err_scale, label=ch)
        c['T'] = template(c['g'], c['d'], c['info'])

    print('red template vs blue template')
    tr = register(chans['red']['T'], chans['blue']['T'], radius=rad, gradient=not a.no_gradient, err_scale=err_scale)
    print('  red sits at dy=%+.3f dx=%+.3f grid px from blue (+- %.3f, %.3f); R^2 %.4f' % (
        tr['d'][0], tr['d'][1], tr['err'][0], tr['err'][1], tr['r2']))

    REF, rr = None, None
    D_abs = np.zeros(2)
    if a.ref:
        REF = highpass(_ref_image(a.ref, a.ref_hdu, wg, shape,
                                  a.ref_smooth / a.gridpix / 2.3548 if a.ref_smooth else 0), hp)
        # blue and red combined at the blue position
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            redb = (nan_shift(chans['red']['T'], -tr['d']) - tr['b']) / tr['a']
            comb = chans['blue']['T'] if a.ref_channel == 'blue' else \
                (redb if a.ref_channel == 'red' else np.nanmean([chans['blue']['T'], redb], axis=0))
        print('KCWI (%s) vs the reference image' % a.ref_channel)
        rr = register(comb, REF, radius=int(np.ceil(a.ref_radius / a.gridpix)), gradient=not a.no_gradient,
                      err_scale=err_scale)
        D_abs = rr['d']
        print('  KCWI header WCS is off by dE=%+.3f" dN=%+.3f" (+- %.3f"); R^2 %.4f' % (
            -D_abs[1] * a.gridpix, D_abs[0] * a.gridpix, np.hypot(*rr['err']) * a.gridpix, rr['r2']))
    elif a.astrom or a.abs_offset:
        # A known absolute correction, in the convention of kcwi_astrometry's .astrom.list
        # and kcwi_stack(use_astrom=True): CRVAL1 += dRA/3600 (RA coordinate, not on-sky
        # arcsec), CRVAL2 += dDec/3600. kcwi_astrometry measured it for a stack whose
        # first frame was unshifted, so it is applied to blue frame 0's header position.
        if a.astrom:
            t = ascii.read(a.astrom, format='no_header')
            dra, dde = float(t['col1'][0]), float(t['col2'][0])
        else:
            dra, dde = a.abs_offset
        c = np.array(wg.wcs.crpix) - 1
        ra0, de0 = wg.all_pix2world(c[0], c[1], 0)
        px, py = wg.all_world2pix(ra0 + dra / 3600., de0 + dde / 3600., 0)
        D_astrom = np.array([c[1] - py, c[0] - px])          # (dy, dx) image displacement
        D_abs = D_astrom - chans['blue']['d'][0]
        print('absolute: %s dRA %+.3f" (coordinate) dDec %+.3f", tied to blue frame 0 (%s)' % (
            a.astrom or '--abs-offset', dra, dde, os.path.basename(chans['blue']['fn'][0])))

    res = dict(created=stamp, blue=a.blue, red=a.red, z=a.z, grid_pix_arcsec=a.gridpix,
               band_blue=chans['blue']['band'], band_red=chans['red']['band'],
               red_vs_blue_grid=list(map(float, tr['d'])), red_vs_blue_err=list(map(float, tr['err'])),
               absolute_grid=list(map(float, D_abs)), ref=a.ref, astrom=a.astrom,
               abs_offset=a.abs_offset, frames={})
    for ch in ('blue', 'red'):
        c = chans[ch]
        rows, out = [], []
        for k in range(len(c['fn'])):
            D = c['d'][k] + (tr['d'] if ch == 'red' else 0.) + D_abs
            sx, sy, dra, dde = crpix_shift(c['wcs'][k], D, wg)
            out.append((c['fn'][k].replace('.fits', ''), sx, sy))
            rows.append(dict(frame=os.path.basename(c['fn'][k]), D_grid=list(map(float, D)),
                             err_grid=list(map(float, c['err'][k])), within_channel_grid=list(map(float, c['d'][k])),
                             crpix_shift=[sx, sy], dRA_coord_arcsec=dra * 3600., dDec_arcsec=dde * 3600.,
                             scale=float(c['info'][k]['a']), r2=float(c['info'][k]['r2'])))
            if np.isfinite(c['info'][k]['r2']) and c['info'][k]['r2'] < 0.8:
                print('  [WARNING] %s: poor match (R^2 %.2f) -- check the frame' % (rows[-1]['frame'], c['info'][k]['r2']))
        res['frames'][ch] = rows
        sfn = c['list'].replace('.list', '.shift.list')
        if not a.dry_run:
            _backup(sfn, stamp)
            ascii.write(Table(list(map(list, zip(*out)))), sfn, overwrite=True, format='no_header')
        print('%s %s' % ('would write' if a.dry_run else 'wrote', sfn))
        for o, r in zip(out, rows):
            print('   %-40s %+8.4f %+8.4f   (+-%.3f")' % (os.path.basename(o[0]), o[1], o[2],
                                                       np.hypot(*r['err_grid']) * a.gridpix))

    stem = os.path.join('kcwi_xalign', os.path.basename(a.blue).replace('.list', ''))
    with open(stem + '.xalign.json', 'w') as f:
        json.dump(res, f, indent=1)
    plot_summary(stem + '.xalign.png', res, chans['blue']['T'], chans['red']['T'], tr, REF, rr)
    print('wrote %s.xalign.json and .png' % stem)
    print('\nNext: kcwi_stack both lists with use_astrom=False (the absolute offset is in the shifts),')
    print('then: kcwi_xalign wcs %s %s%s' % (a.blue, a.red, (' --ref ' + a.ref) if a.ref else ''))


def step_wcs(a):
    stamp = time.strftime('%Y%m%d_%H%M%S')
    suffix = 'cubed' if a.cubed else 'cubes'
    sb = a.blue.replace('.list', '_i' + suffix + a.stack_suffix + '.fits')
    sr = a.red.replace('.list', '_i' + suffix + a.stack_suffix + '.fits')
    if a.stack_suffix and not a.check_only:
        sys.exit('--stack-suffix is for measuring only: add --check-only')
    hb, hr = fits.getheader(sb), fits.getheader(sr)
    wb, wr = awcs.WCS(hb).celestial, awcs.WCS(hr).celestial
    same_shape = (hb['NAXIS1'], hb['NAXIS2']) == (hr['NAXIS1'], hr['NAXIS2'])
    corners = np.array([[0, 0], [hb['NAXIS1'] - 1, 0], [0, hb['NAXIS2'] - 1], [hb['NAXIS1'] - 1, hb['NAXIS2'] - 1]], float)
    pr = wr.all_world2pix(wb.all_pix2world(corners, 0), 0) if same_shape else None
    off = np.abs(pr - corners).max() if same_shape else np.inf
    if off > 0.02:
        print('[WARNING] the blue and red stacks are NOT on the same pixel grid (max %.3f px apart%s).' % (
            off, '' if same_shape else ', different sizes'))
        print('          Same stack_ad / dimension / xpix / ypix / orientation in both .par files, then re-stack.')
    else:
        print('blue and red stacks share one pixel grid (to %.4f px)' % off)
    # measure the stacks against each other
    ib, w2b, bandb = frame_image(sb, a.wb, a.z, 'blue', npix_trim=1, use_regmask=False)
    ir, w2r, bandr = frame_image(sr, a.wr, a.z, 'red', npix_trim=1, use_regmask=False)
    pix = np.sqrt(np.abs(np.linalg.det(w2b.pixel_scale_matrix))) * 3600.
    sm = a.smooth / pix / 2.3548 if a.smooth else 0
    hp = a.highpass / pix / 2.3548 if a.highpass else 0
    ib, ir = highpass(nan_smooth(ib, sm), hp), highpass(nan_smooth(ir, sm), hp)
    if off <= 0.02:
        r = register(ir, ib, radius=int(np.ceil(a.radius / pix)), gradient=not a.no_gradient)
        print('red vs blue stack: dx=%+.3f dy=%+.3f px (+- %.3f, %.3f, formal), R^2 %.4f' % (
            r['d'][1], r['d'][0], r['err'][1], r['err'][0], r['r2']))
        if np.hypot(*r['d']) > 0.15:
            print('[WARNING] residual blue-red offset above 0.15 px: rerun `shifts` (or check the windows)')
    if a.ref:
        REF = highpass(_ref_image(a.ref, a.ref_hdu, w2b, ib.shape,
                                  a.ref_smooth / pix / 2.3548 if a.ref_smooth else 0), hp)
        r = register(ib, REF, radius=int(np.ceil(a.ref_radius / pix)), gradient=not a.no_gradient)
        print('blue stack vs reference: dx=%+.3f dy=%+.3f px (%.3f")' % (r['d'][1], r['d'][0], np.hypot(*r['d']) * pix))
    if a.check_only:
        return
    for s in (sb, sr):
        out = s.replace('.fits', '_wcs.fits')
        if os.path.isfile(out) and not a.overwrite:
            print('%s exists; --overwrite to replace (a backup is kept)' % out)
            continue
        _backup(out, stamp)
        with fits.open(s) as hl:
            hl[0].header['HISTORY'] = 'kcwi_xalign %s: WCS from cross-correlated .shift.list (no point source)' % stamp
            hl.writeto(out, overwrite=True)
        print('wrote', out)


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['shifts', 'wcs'])
    ap.add_argument('blue', help='blue .list (as for kcwi_align / kcwi_stack), or a source name NAME '
                                 '(lists NAME<suffix-blue>.list and NAME<suffix-red>-red.list)')
    ap.add_argument('red', nargs='?', help='red .list (omit when blue is a source name)')
    ap.add_argument('--suffix-blue', default='', help="with a source name: blue list NAME<suffix>.list, e.g. '_skl'")
    ap.add_argument('--suffix-red', default='', help='with a source name: red list NAME<suffix>-red.list')
    ap.add_argument('--z', type=float, help='redshift: cut emission lines out of the white-light windows')
    ap.add_argument('--wb', type=float, nargs=2, help='blue white-light window [A]; default: 400 A below WAVGOOD1')
    ap.add_argument('--wr', type=float, nargs=2, help='red white-light window [A]; default: 400 A above WAVGOOD0')
    ap.add_argument('--gridpix', type=float, default=0.1, help='common grid pixel [arcsec] (default 0.1)')
    ap.add_argument('--radius', type=float, default=2.0, help='search radius for frame/channel offsets [arcsec]')
    ap.add_argument('--smooth', type=float, default=0., help='Gaussian FWHM [arcsec] applied to all KCWI images')
    ap.add_argument('--no-gradient', action='store_true', help='fit only scale + offset, no background gradient')
    ap.add_argument('--highpass', type=float, default=0., metavar='FWHM',
                    help='subtract a Gaussian-smoothed copy of this FWHM [arcsec] from every image, to suppress '
                         'large-scale sky-subtraction differences between frames and channels (e.g. 3)')
    ap.add_argument('--ref', help='external image for the absolute WCS (any FITS with a celestial WCS)')
    ap.add_argument('--astrom', help='apply a known absolute correction: an .astrom.list from kcwi_astrometry '
                    '(e.g. the old point-source solution); used for BOTH channels; ignored with --ref')
    ap.add_argument('--abs-offset', type=float, nargs=2, metavar=('dRA', 'dDec'),
                    help='the same, given directly [arcsec; dRA in RA coordinate, as in .astrom.list]')
    ap.add_argument('--ref-hdu', type=int, default=0)
    ap.add_argument('--ref-radius', type=float, default=15., help='search radius against --ref [arcsec]')
    ap.add_argument('--ref-smooth', type=float, default=0., help='Gaussian FWHM [arcsec] for --ref (match the seeing)')
    ap.add_argument('--ref-channel', choices=['blue', 'red', 'both'], default='both')
    ap.add_argument('--cubed', action='store_true', help='use *_icubed instead of *_icubes')
    ap.add_argument('--npix-trim', type=int, default=1)
    ap.add_argument('--no-regmask', action='store_true', help='ignore *_rm.reg region masks')
    ap.add_argument('--dry-run', action='store_true', help='shifts: do not write the .shift.list files')
    ap.add_argument('--overwrite', action='store_true', help='wcs: replace existing *_wcs.fits')
    ap.add_argument('--check-only', action='store_true', help='wcs: only measure, write nothing')
    ap.add_argument('--stack-suffix', default='', help="wcs: measure e.g. the '_wcs' stacks instead")
    return ap


def source_lists(name, suffix_blue='', suffix_red=''):
    """NAME -> (NAME<sb>.list, NAME<sr>-red.list)."""
    b = name[:-len('.list')] if name.endswith('.list') else name
    return b + (suffix_blue or '') + '.list', b + (suffix_red or '') + '-red.list'


def parse(argv):
    ap = build_parser()
    a = ap.parse_args(argv)
    if a.red is None:
        a.blue, a.red = source_lists(a.blue, a.suffix_blue, a.suffix_red)
    elif not a.blue.endswith('.list'):
        ap.error('give two .list files, or one source name (with --suffix-blue / --suffix-red)')
    for l in (a.blue, a.red):
        if not os.path.isfile(l):
            sys.exit('no list %s' % l)
    return a


# ---------------------------------------------------------------------------
# batch: many sources from a file
# ---------------------------------------------------------------------------
BATCH_DOC = """kcwi_xalign batch SOURCES [--steps shifts,stack,wcs] [--suffix-blue S] [--suffix-red S]
                  [--channels blue|red|both] [--stack-kw KEY=VAL,...] [common kcwi_xalign options]

SOURCES: one source per line ('#' comments). The first word is a source name NAME
(lists NAME<suffix-blue>.list and NAME<suffix-red>-red.list) or two .list files;
the rest of the line is options for that source only (e.g. its --z, --ref),
overriding the common ones. The kcwi_skylevel batch file can be reused: options
this script does not know are an error, so keep per-source --z etc. in a
separate file if the two need different ones.

    # source            per-source options
    IRAS03582+6012      --z 0.030
    LEDA66688-final     --z 0.0216
    UGC2369N            --z 0.031 --ref ugc2369_f814w.fits

Steps (in this order, any subset):
  shifts  cross-correlation shift lists for both channels
  stack   kcwi_stack each channel in --channels (default both), with --stack-kw
          (default crr_split=True,crr_split_mode=frame,crr_save_mask=True,overwrite=True)
  wcs     check blue vs red and write the _icubes_wcs files (add --overwrite)

After kcwi_skylevel apply corrected one channel (e.g. blue -> NAME_skl.list, red
unchanged), do NOT re-run shifts: apply copied the shift list. Restack the
corrected channel and check:
    kcwi_xalign batch sources.txt --steps stack,wcs --suffix-blue _skl --channels blue --overwrite

Each source logs to <NAME><suffix-blue>.xalign.log; a summary is printed at the end
and a failing source does not stop the others.
"""


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for x in self.streams:
            x.write(s)

    def flush(self):
        for x in self.streams:
            x.flush()


def _kw(spec):
    out = {}
    for item in [x for x in (spec or '').split(',') if x.strip()]:
        k, v = item.split('=', 1)
        v = v.strip()
        for conv in (int, float):
            try:
                v = conv(v); break
            except ValueError:
                pass
        if v in ('True', 'False', 'None'):
            v = {'True': True, 'False': False, 'None': None}[v]
        out[k.strip()] = v
    return out


def batch(argv):
    import contextlib, io, re, traceback
    bp = argparse.ArgumentParser(prog='kcwi_xalign batch', description=BATCH_DOC,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    bp.add_argument('sources')
    bp.add_argument('--steps', default='shifts', help='any of shifts, stack, wcs (comma-separated; default shifts)')
    bp.add_argument('--suffix-blue', default='')
    bp.add_argument('--suffix-red', default='')
    bp.add_argument('--channels', choices=['blue', 'red', 'both'], default='both', help='stack: which channels')
    bp.add_argument('--stack-kw', default='crr_split=True,crr_split_mode=frame,crr_save_mask=True,overwrite=True',
                    help='keyword arguments for kcwi_stack')
    b, common = bp.parse_known_args(argv)
    steps = [x.strip() for x in b.steps.split(',') if x.strip()]
    if not steps or any(x not in ('shifts', 'stack', 'wcs') for x in steps):
        sys.exit('--steps: any of shifts, stack, wcs')
    steps = [x for x in ('shifts', 'stack', 'wcs') if x in steps]
    if 'shifts' in steps and (b.suffix_blue != b.suffix_red):
        print('NOTE: shifts with different blue/red suffixes rewrites the shift list of the unchanged channel too; '
              'after kcwi_skylevel apply, use --steps stack,wcs instead')
    skw = _kw(b.stack_kw)
    rows = []
    for line in open(b.sources):
        t = line.split('#')[0].split()
        if t:
            rows.append(t)
    res = []
    for row in rows:
        if len(row) >= 2 and row[0].endswith('.list') and row[1].endswith('.list'):
            blue, red, extra = row[0], row[1], row[2:]
        else:
            blue, red = source_lists(row[0], b.suffix_blue, b.suffix_red)
            extra = row[1:]
        name = os.path.basename(blue[:-len('.list')])
        buf = io.StringIO()
        r = dict(source=name, status='ok', grid='-', offset='-', warnings=0)
        with contextlib.redirect_stdout(_Tee(buf, sys.stdout)):
            try:
                print('=' * 100 + '\n%s  (%s, %s)  %s' % (name, blue, red, ' '.join(extra)))
                for l in (blue, red):
                    if not os.path.isfile(l):
                        sys.exit('no list %s' % l)
                for step in steps:
                    print('-- %s' % step)
                    if step == 'stack':
                        from kcwikit.kcwi.kcwi import kcwi_stack
                        for ch, l in (('blue', blue), ('red', red)):
                            if b.channels in (ch, 'both'):
                                kcwi_stack(l, **skw)
                        continue
                    a = parse([step, blue, red] + common + extra)
                    step_shifts(a) if step == 'shifts' else step_wcs(a)
            except SystemExit as e:
                r['status'] = 'FAILED: %s' % e
                print('FAILED:', e)
            except Exception as e:
                r['status'] = 'FAILED: %s: %s' % (type(e).__name__, e)
                traceback.print_exc(file=sys.stdout)
        log = buf.getvalue()
        r['warnings'] = sum(1 for l in log.splitlines() if 'WARNING' in l)
        m = re.findall(r'red vs blue stack: dx=([+-][\d.]+) dy=([+-][\d.]+) px', log)
        if m:
            r['offset'] = 'dx %s dy %s px' % m[-1]
        m = re.findall(r'share one pixel grid \(to ([\d.]+) px\)', log)
        if m:
            r['grid'] = 'same (%s px)' % m[-1]
        elif 'NOT on the same pixel grid' in log:
            r['grid'] = 'DIFFERENT'
        try:
            open(blue[:-len('.list')] + '.xalign.log', 'w').write(log)
        except OSError:
            pass
        res.append(r)
    print('\n%-30s %-16s %-24s %4s  %s' % ('source', 'grid (wcs)', 'red vs blue stack', 'warn', 'status'))
    for r in res:
        print('%-30s %-16s %-24s %4s  %s' % (r['source'], r['grid'], r['offset'], r['warnings'], r['status']))
    print('logs: <source>.xalign.log')
    if any(r['status'] != 'ok' for r in res):
        sys.exit(1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'batch':
        return batch(sys.argv[2:])
    a = parse(sys.argv[1:])
    step_shifts(a) if a.step == 'shifts' else step_wcs(a)


if __name__ == '__main__':
    main()
