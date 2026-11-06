#!/usr/bin/env python
"""
kcwi_skylevel -- match the sky level of the KCWI blue and red channels.

Each channel's sky subtraction can leave an additive error: an offset-sky frame
scaled by exposure time when the sky changed, a "sky" pointing that contains
galaxy light, ZAP's PCA taking up continuum of a galaxy that fills the field.
The error is the same in every spaxel, so it is invisible in bright regions and
dominates faint ones, and it differs between the channels, which makes a step
at the dichroic. This script measures it and puts it back, per channel:

    error(lambda) = beta_ch x T(lambda)

with T the continuum of the sky that was subtracted from each frame (smoothed
over sky lines), taken from the frame's own sky extension:
    SKYMODEL_ZAP (KSkyWizard/ZAP), SKY / SKYMODEL, or NOSKYSUB - data (KCWI_DRP);
or, with --template flat / counts, a constant in flux / in detector counts.
beta > 0: too much sky was subtracted (beta = 1: all of it); beta < 0: too little.

Step 1 -- measure (on the current stacks; needs the frames for the template):

    python kcwi_skylevel.py measure UGC2369S.list UGC2369S-red.list \
        --image ugc2369i.fits,7000,8950,hst  [--image legacy_g.fits,4100,5400,nmgy,4800] \
        [--cubewin ugc2369s_muse_windows.npz]

  Constraints, each one linear in (beta_blue, beta_red):
   * dichroic continuity (always): per spaxel F_blue_end = k F_red_start + c;
     c = k beta_red T_red - beta_blue T_blue must be 0.
   * absolute anchors (optional): KCWI = a x reference + b in a window, b = -beta T.
     --image FILE,LO,HI,UNIT[,PIVOT]  UNIT = hst (SCI e-/s x PHOTFLAM), nmgy
                                       (Legacy Survey / SDSS nanomaggies, needs PIVOT [A]),
                                       ps1 (Pan-STARRS stack cutout, needs PIVOT:
                                       g 4866, r 6215, i 7545), flam (erg/s/cm2/A per pixel)
     LO-HI = the KCWI window [A observed] compared with it (inside the band).
     --cubewin NPZ                    windows from cube_windows.py (e.g. MUSE)
  Which channel is trusted:
     --fixed red   (default) keep red as it is (beta_red = 0), solve for beta_blue
     --fixed blue  keep blue as it is (beta_blue = 0), solve for beta_red
     --fixed none  solve for both (needs at least one anchor, or --same)
     --same        one beta for both channels (both used the same sky pointing)
     --groups auto|FILE  one beta per group: a blue frame and the red frames taken at
                   the same pointing and PA (auto), or one group per line in FILE
                   (frame names). Measured on the frames themselves, for data where
                   the sky changed between exposures (e.g. twilight). apply then
                   corrects each frame with its group's beta (frame mode only).
  The dichroic constraint alone fixes only the difference between the channels,
  so without anchors one channel has to be fixed (or --same).
  Writes <blue base>.skylevel.json and prints every constraint with its residual.

Step 2 -- apply, then re-stack:

    python kcwi_skylevel.py apply UGC2369S.list UGC2369S-red.list [--suffix _skl] [--outdir DIR]

  For each corrected channel writes the frames as <frame>_skl_icubes.fits (all
  extensions kept; header SKLBETA, SKLTEMPL), the list UGC2369S_skl.list (or
  UGC2369S_skl-red.list) and copies of its .par / .shift.list / .preshift.list /
  .flx.list / .astrom.list. The fixed channel (beta = 0) is left alone: keep its
  original list and stack. Stack the new list as usual (kcwi_stack, kcwi_xalign wcs).
  `apply --stack` instead corrects the existing stack directly (quick check,
  no re-stack): <stack>_skl_icubes_wcs.fits etc.

Many sources at once:  kcwi_skylevel batch SOURCES [--steps measure,apply] [--jobs N] [options]
  (one source per line, optionally followed by its own options; `kcwi_skylevel batch -h`).

Check afterwards with kcwi_dichroic_step (c ~ 0; --red to name an unchanged red
stack) and gist_fitresid.py.
Run in the folder that holds the lists (the stack folder).
"""
import argparse
import json
import os
import re
import shutil
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter

warnings.simplefilter('ignore')
SKYLINES = [5577.34, 5889.95, 5895.92, 6300.30, 6363.78, 7316.3, 7340.9, 7750.6, 7794.1, 7821.5,
            7913.7, 7993.3, 8344.6, 8399.2, 8430.2, 8465.4, 8827.1, 8885.8, 8919.6]


# ---------------------------------------------------------------------------
# lists, frames, stacks
# ---------------------------------------------------------------------------
PATHMAP = []          # [(from, to)], set by --path-map


def local(p):
    for a, b in PATHMAP:
        if p.startswith(a):
            return b + p[len(a):]
    return p


def read_list(fn):
    rows = []
    for line in open(fn):
        s = line.split('#')[0].split()
        if s:
            rows.append(s)
    return rows


def frame_file(base):
    for suf in ('cubes', 'cubed'):
        f = base + '_i' + suf + '.fits'
        if os.path.isfile(f):
            return f
    sys.exit('no %s_icubes.fits' % base)


def wave(h):
    step = h.get('CD3_3', h.get('CDELT3'))
    return h['CRVAL3'] + (np.arange(h['NAXIS3']) + 1 - h['CRPIX3']) * step


def stack_files(lst):
    b = lst[:-len('.list')] if lst.endswith('.list') else lst
    for s in ('_icubes_wcs.fits', '_icubes.fits'):
        if os.path.isfile(b + s):
            return b, b + s, (b + '_mcubes.fits') if os.path.isfile(b + '_mcubes.fits') else None
    sys.exit('no stack %s_icubes(_wcs).fits for %s' % (b, lst))


def spaxel_area(h):
    w = WCS(h).celestial
    return abs(np.linalg.det(w.pixel_scale_matrix)) * 3600. ** 2


def smooth_continuum(w, s, width=60.):
    """Running median over `width` A after cutting sky lines: the sky continuum."""
    s = np.array(s, float)
    for l in SKYLINES:
        s[np.abs(w - l) < 6] = np.nan
    good = np.isfinite(s)
    if good.sum() < 10:
        return np.zeros_like(s)
    s = np.interp(w, w[good], s[good])
    n = max(3, int(width / np.median(np.diff(w))) | 1)
    return median_filter(s, size=n, mode='nearest')


def frame_template(fn, kind='sky'):
    """1-D template (frame units per spaxel) and a note on where it came from."""
    h = fits.open(fn, memmap=True)
    H = h[0].header
    w = wave(H)
    names = [e.name for e in h]
    nz = H['NAXIS3']
    if kind == 'flat':
        return w, np.ones(nz), 'flat'
    if kind == 'counts':
        cnt = fn.replace('_icubes', '_icube')
        if not os.path.isfile(cnt):
            sys.exit('--template counts needs %s' % cnt)
        c = fits.open(cnt, memmap=True)[0].data
        conv = np.zeros(nz)
        for a in range(0, nz, 300):
            r = np.array(h[0].data[a:a + 300], float) / np.array(c[a:a + 300], float)
            r[~np.isfinite(r)] = np.nan
            conv[a:a + 300] = np.nanmedian(r.reshape(r.shape[0], -1), axis=1)
        return w, smooth_continuum(w, conv, 30.), 'counts (flux per e-)'
    if 'SKYMODEL_ZAP' in names and 'NOSKYSUB' in names and 'UNZAPPED' in names and H.get('SKYCOR', False):
        # DRP subtracted a sky (SKYMAST), then ZAP ran on the result. Fine if ZAP only removed residuals
        # (its sky cube was DRP-subtracted too: ZAP continuum ~ 0); a double subtraction if ZAP's
        # continuum is comparable to the DRP sky (its sky cube was NOT sky-subtracted).
        m0, m1 = max(0, nz // 2 - 150), min(nz, nz // 2 + 150)
        zc = np.nanmedian(np.array(h['SKYMODEL_ZAP'].data[m0:m1], float))
        dc = np.nanmedian(np.array(h['NOSKYSUB'].data[m0:m1], float) - np.array(h['UNZAPPED'].data[m0:m1], float))
        frac = zc / dc if dc > 0 else np.nan
        if not np.isfinite(frac) or frac > 0.3:
            print('WARNING %s: sky subtracted twice -- by the DRP (SKYMAST %s) and again by ZAP, whose continuum '
                  'is %.0f%% of the DRP sky (its sky cube was not sky-subtracted). Re-run with "skip" for this '
                  'frame in kcwi.sky, or give ZAP a DRP-subtracted sky cube.' % (os.path.basename(fn), H.get('SKYMAST'), 100 * frac))
        src, get = ('ZAP + DRP sky' if not np.isfinite(frac) or frac > 0.3 else 'DRP sky + ZAP resid.'), \
            (lambda a, b, st: np.array(h['SKYMODEL_ZAP'].data[a:b:st], float)
                                     + np.array(h['NOSKYSUB'].data[a:b:st], float)
                                     - np.array(h['UNZAPPED'].data[a:b:st], float))
    elif 'SKYMODEL_ZAP' in names:
        src, get = 'SKYMODEL_ZAP', (lambda a, b, st: np.array(h['SKYMODEL_ZAP'].data[a:b:st], float))
    elif 'SKY' in names or 'SKYMODEL' in names:
        e = 'SKY' if 'SKY' in names else 'SKYMODEL'
        src, get = e, (lambda a, b, st: np.array(h[e].data[a:b:st], float))
    elif 'NOSKYSUB' in names:
        src, get = 'NOSKYSUB - data', (lambda a, b, st: np.array(h['NOSKYSUB'].data[a:b:st], float)
                                       - np.array(h[0].data[a:b:st], float))
    else:
        sys.exit('%s has no sky extension (SKYMODEL_ZAP, SKY, NOSKYSUB): use --template flat or counts' % fn)
    flags = h['FLAGS'].data if 'FLAGS' in names else None
    st = max(1, int(round(TEMPLATE_STEP / abs(H.get('CD3_3', 1.))))) if TEMPLATE_STEP > 0 else 1
    idx = np.arange(0, nz, st)
    prof_s = np.zeros(idx.size)
    blk = 300 * st                                    # ~300 sampled planes per read
    for a in range(0, nz, blk):
        s = get(a, a + blk, st)
        if flags is not None:
            s[np.array(flags[a:a + blk:st]) != 0] = np.nan
        s[s == 0] = np.nan
        i0 = a // st
        prof_s[i0:i0 + s.shape[0]] = np.nanmedian(s.reshape(s.shape[0], -1), axis=1)
    if st > 1:
        g = np.isfinite(prof_s)
        prof = np.interp(w, w[idx][g], prof_s[g]) if g.sum() > 1 else np.full(nz, np.nan)
    else:
        prof = prof_s
    # counts vs flux sanity check: counts would be ~1e4-1e6 x the flux values.
    # Ignore zeros (DAR padding / masked voxels are often >50% of the cube).
    dsl = np.abs(np.array(h[0].data[nz // 2 - 50:nz // 2 + 50], float))
    dsl = dsl[np.isfinite(dsl) & (dsl > 0)]
    dmed = np.median(dsl) if dsl.size else np.nan
    if np.isfinite(dmed) and np.nanmedian(np.abs(prof)) > 1e3 * dmed:
        print('WARNING %s: the sky extension looks like counts, not flux; template may be wrong' % os.path.basename(fn))
    return w, smooth_continuum(w, prof), src


CACHE = {}
TEMPLATE_STEP = 2.0   # [A] sample the sky extension every ~2 A (the template is smoothed over 60 A); 0 = every plane


def cached_template(fn, kind):
    key = '%s|%s|%d|s%g' % (os.path.abspath(fn), kind, int(os.path.getmtime(fn)), TEMPLATE_STEP)
    if key not in CACHE:
        CACHE[key] = frame_template(fn, kind)
    return CACHE[key]


def load_cache(fn):
    if fn and os.path.isfile(fn):
        try:
            z = json.load(open(fn))
            for k, v in z.items():
                CACHE[k] = (np.array(v[0]), np.array(v[1]), v[2])
        except Exception:
            pass


def save_cache(fn):
    if fn:
        json.dump({k: [list(map(float, v[0])), list(map(float, v[1])), v[2]] for k, v in CACHE.items()}, open(fn, 'w'))


def channel_template(lst, kind):
    """Per frame: (base, file, wave, template per spaxel, area, exptime, source)."""
    out = []
    for row in read_list(lst):
        fn = frame_file(local(row[0]))
        w, t, src = cached_template(fn, kind)
        H = fits.getheader(fn)
        out.append(dict(base=row[0], file=fn, wave=w, T=t, area=spaxel_area(H),
                        exptime=float(H.get('XPOSURE', H.get('EXPTIME', 1.))), src=src))
        print('  %-40s template %-15s median %.4g per spaxel (%.4g per arcsec2)' % (
            os.path.basename(fn), src, np.nanmedian(t), np.nanmedian(t) / out[-1]['area']))
    return out


def stack_template(frames, wst):
    """Exposure-weighted template in stack units (per arcsec2) on the stack grid."""
    num = np.zeros(wst.size); den = 0.
    for f in frames:
        num += f['exptime'] * np.interp(wst, f['wave'], f['T'] / f['area'], left=np.nan, right=np.nan)
        den += f['exptime']
    return num / den


def window_image(cube, mcube, lo, hi):
    h = fits.open(cube, memmap=True)
    H = h[0].header
    w = wave(H)
    k = (w >= lo) & (w <= hi)
    for l in SKYLINES:
        k &= np.abs(w - l) > 8
    i = np.where(k)[0]
    if i.size == 0:
        return None, None
    d = np.array(h[0].data[i[0]:i[-1] + 1], float)[k[i[0]:i[-1] + 1]]
    if mcube:
        m = np.array(fits.open(mcube, memmap=True)[0].data[i[0]:i[-1] + 1])[k[i[0]:i[-1] + 1]]
        d[m != 0] = np.nan
    d[d == 0] = np.nan
    return np.nanmedian(d, axis=0), k


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------
def robust_line(x, y, nboot=200):
    g = np.isfinite(x) & np.isfinite(y)
    ok = g.copy()
    for _ in range(10):
        A = np.vstack([x[ok], np.ones(ok.sum())]).T
        (a, b), *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        r = y - a * x - b
        s = 1.4826 * np.nanmedian(np.abs(r[ok]))
        new = g & (np.abs(r) < 3 * s)
        if (new == ok).all():
            break
        ok = new
    rng = np.random.default_rng(1)
    ii = np.where(ok)[0]
    bs = np.array([np.linalg.lstsq(np.vstack([x[j], np.ones(j.size)]).T, y[j], rcond=None)[0]
                   for j in (rng.choice(ii, ii.size) for _ in range(nboot))])
    return a, b, bs[:, 0].std(), bs[:, 1].std(), int(ok.sum())


def unit_scale(bunit):
    """Factor to 1e-16 erg/s/cm2/A from a BUNIT string like 10**(-20)erg... or 1e-20 erg..."""
    m = re.search(r'10\*\*\(?(-?\d+)\)?', bunit or '') or re.search(r'1e([+-]?\d+)', bunit or '')
    return 10. ** (int(m.group(1)) + 16) if m else 1.


def reference_on_grid(img, rw, kw, shape, fwhm, pix, register=True, kimg=None, maxshift=2.0):
    """Smooth a reference image to fwhm and sample it on the KCWI stack grid (optionally registered)."""
    foot = np.isfinite(img).astype(float)
    img = np.nan_to_num(img)
    sig = fwhm / 2.3548 / pix
    sm = gaussian_filter(img, sig)
    fs = gaussian_filter(foot, sig)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    ra, de = kw.all_pix2world(xx, yy, 0)
    px, py = rw.all_world2pix(ra, de, 0)
    best = (0, 0, -2)
    if register and kimg is not None:
        n = int(maxshift / pix)
        stp = max(1, n // 12)
        for rnd in range(2):
            cx, cy = best[0], best[1]
            rng_ = range(-n, n + 1, stp) if rnd == 0 else range(-stp, stp + 1)
            for dx in rng_:
                for dy in rng_:
                    sx, sy = (dx, dy) if rnd == 0 else (cx + dx, cy + dy)
                    mm = map_coordinates(sm, [py + sy, px + sx], order=1, cval=np.nan)
                    ff = map_coordinates(fs, [py + sy, px + sx], order=1, cval=0)
                    g = np.isfinite(mm) & np.isfinite(kimg) & (ff > 0.99)
                    if g.sum() < 200:
                        continue
                    cc = np.corrcoef(mm[g], kimg[g])[0, 1]
                    if cc > best[2]:
                        best = (sx, sy, cc)
    mm = map_coordinates(sm, [py + best[1], px + best[0]], order=1, cval=np.nan)
    ff = map_coordinates(fs, [py + best[1], px + best[0]], order=1, cval=0)
    mm[ff < 0.99] = np.nan
    return mm, (best[0] * pix, best[1] * pix, best[2])


def load_image(spec):
    p = spec.split(',')
    fn, lo, hi, unit = p[0], float(p[1]), float(p[2]), p[3].lower()
    h = fits.open(fn)
    e = 'SCI' if 'SCI' in h else next(i for i, x in enumerate(h) if x.data is not None and x.data.ndim == 2)
    H = h[e].header
    d = np.array(h[e].data, float)
    w = WCS(H).celestial
    pix = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600.
    if unit == 'hst':
        flam = H.get('PHOTFLAM') or h[0].header.get('PHOTFLAM')
        d = d * flam
    elif unit == 'nmgy':
        lam = float(p[4])
        d = d * 3.631e-29 * 2.99792458e18 / lam ** 2
    elif unit == 'ps1':
        # Pan-STARRS stack cutout: asinh-compressed counts (BSOFTEN, BOFFSET), ZP = 25 + 2.5 log10(EXPTIME)
        lam = float(p[4])
        if 'BSOFTEN' in H and 'BOFFSET' in H:
            d = H['BOFFSET'] + H['BSOFTEN'] * 2. * np.sinh(d * np.log(10.) / 2.5)
        zp = 25. + 2.5 * np.log10(H.get('EXPTIME', h[0].header.get('EXPTIME', 1.)))
        d = d * 10. ** (-0.4 * (zp + 48.6)) * 2.99792458e18 / lam ** 2
    elif unit != 'flam':
        sys.exit('unknown unit %s (hst, nmgy, ps1, flam)' % unit)
    return dict(name=os.path.basename(fn), lo=lo, hi=hi, img=d / pix ** 2 / 1e-16, wcs=w, pix=pix)


# ---------------------------------------------------------------------------
# groups: one blue frame + the red frames taken at the same pointing
# ---------------------------------------------------------------------------
def frame_pointing(fn):
    H = fits.getheader(fn)
    return (float(H.get('CRVAL1', np.nan)), float(H.get('CRVAL2', np.nan)), float(H.get('ROTPOSN', np.nan)),
            str(H.get('DATE-BEG', '')))


def make_groups(a, chans):
    """[(name, [blue frames], [red frames])]; from --groups FILE or by pointing (--groups auto)."""
    fr = {side: chans[side]['fr'] for side in ('blue', 'red')}
    groups = []
    if a.groups != 'auto':
        # one group per line: any frame names (base names, with or without path / suffix)
        for i, row in enumerate(read_list(a.groups)):
            g = {'blue': [], 'red': []}
            for tok in row:
                key = os.path.basename(tok).replace('_icubes.fits', '').replace('.fits', '')
                for side in ('blue', 'red'):
                    g[side] += [f for f in fr[side] if os.path.basename(f['base']) == key
                                or os.path.basename(f['base']).startswith(key + '_')]
            groups.append(('group %d' % (i + 1), g['blue'], g['red']))
    else:
        tol, rtol = 1. / 3600., 0.5
        pts = []
        for side in ('blue', 'red'):
            for f in fr[side]:
                f['pt'] = frame_pointing(f['file'])
                pts.append((side, f))
        for side, f in sorted(pts, key=lambda x: x[1]['pt'][3]):
            ra, de, rot, _ = f['pt']
            for g in groups:
                gra, gde, grot = g[3]
                if (abs((ra - gra) * np.cos(np.radians(de))) < tol and abs(de - gde) < tol
                        and abs(((rot - grot) + 180.) % 360. - 180.) < rtol):
                    g[1 if side == 'blue' else 2].append(f)
                    break
            else:
                groups.append(['', [f] if side == 'blue' else [], [f] if side == 'red' else [], (ra, de, rot)])
        groups = [('PA %.1f, %s' % (g[3][2], (g[1] or g[2])[0]['pt'][3][11:19]), g[1], g[2]) for g in groups]
    used = {side: set() for side in ('blue', 'red')}
    out = []
    for name, gb, gr in groups:
        if not gb or not gr:
            print('WARNING %s: %s has no %s frames; not measured' % (
                name, ', '.join(os.path.basename(f['base']) for f in gb + gr), 'red' if gb else 'blue'))
            continue
        for side, gl in (('blue', gb), ('red', gr)):
            for f in gl:
                if f['base'] in used[side]:
                    sys.exit('%s is in more than one group' % f['base'])
                used[side].add(f['base'])
        out.append((name, gb, gr))
    for side in ('blue', 'red'):
        left = [os.path.basename(f['base']) for f in fr[side] if f['base'] not in used[side]]
        if left:
            print('WARNING %s frames in no group (beta = 0): %s' % (side, ', '.join(left)))
    return out


def frame_window(f, lo, hi):
    """Median image of one frame in a window, per arcsec2 (flagged and zero voxels left out)."""
    h = fits.open(f['file'], memmap=True)
    w = wave(h[0].header)
    k = (w >= lo) & (w <= hi)
    for l in SKYLINES:
        k &= np.abs(w - l) > 8
    i = np.where(k)[0]
    d = np.array(h[0].data[i[0]:i[-1] + 1], float)[k[i[0]:i[-1] + 1]]
    for e in ('FLAGS', 'MASK'):
        if e in h:
            m = np.array(h[e].data[i[0]:i[-1] + 1])[k[i[0]:i[-1] + 1]]
            d[m != 0] = np.nan
    d[d == 0] = np.nan
    return np.nanmedian(d, axis=0) / f['area'], WCS(h[0].header).celestial


def same_grid(w1, w2, shape1, shape2):
    if shape1 != shape2:
        return False
    return (np.allclose(w1.wcs.crval, w2.wcs.crval, rtol=0, atol=1e-6 / 3600.)      # 1 micro-arcsec
            and np.allclose(w1.wcs.crpix, w2.wcs.crpix, rtol=0, atol=1e-3)
            and np.allclose(w1.pixel_scale_matrix, w2.pixel_scale_matrix, rtol=1e-6, atol=1e-12))


def template_window(f, lo, hi):
    k = (f['wave'] >= lo) & (f['wave'] <= hi)
    return float(np.nanmedian(f['T'][k])) / f['area']


def measure_groups(a, chans, wb, wr):
    mode = 'same' if a.same else a.fixed
    if mode == 'none':
        sys.exit('--groups solves one beta per group from its dichroic step: use --fixed red|blue or --same')
    if a.image or a.cubewin:
        print('WARNING --groups: anchors (--image / --cubewin) are not used per group; ignored')
    groups = make_groups(a, chans)
    if not groups:
        sys.exit('no group has both blue and red frames')
    print('\ngroups (dichroic %.0f-%.0f / %.0f-%.0f A; c = k beta_red T_red - beta_blue T_blue):' % (wb + wr))
    beta = {'blue': {}, 'red': {}}
    res = []
    for name, gb, gr in groups:
        # blue: exposure-weighted mean of the group's blue frames (normally one)
        imgs = []
        for side, gl, win in (('blue', gb, wb), ('red', gr, wr)):
            ref_w, stack_ = None, []
            for f in gl:
                im, wcs_ = frame_window(f, *win)
                if ref_w is None:
                    ref_w, shape = wcs_, im.shape
                elif not same_grid(ref_w, wcs_, shape, im.shape):
                    pix = np.sqrt(abs(np.linalg.det(wcs_.pixel_scale_matrix))) * 3600.
                    im, _ = reference_on_grid(im, wcs_, ref_w, shape, 0.01, pix, register=False)
                stack_.append((f['exptime'], im))
            wsum = sum(t for t, _ in stack_)
            with np.errstate(invalid='ignore'):
                img = np.nansum([t * np.nan_to_num(im) for t, im in stack_], axis=0) / \
                    np.nansum([t * np.isfinite(im) for t, im in stack_], axis=0)
            imgs.append((img, ref_w, shape))
        (ib, wcb, shb), (ir, wcr, shr) = imgs
        if not same_grid(wcb, wcr, shb, shr):     # red onto the blue frame grid
            pix = np.sqrt(abs(np.linalg.det(wcr.pixel_scale_matrix))) * 3600.
            ir, _ = reference_on_grid(ir, wcr, wcb, shb, 0.01, pix, register=False)
        k, c, ek, ec, n = robust_line(ir.ravel(), ib.ravel())
        tb = np.average([template_window(f, *wb) for f in gb], weights=[f['exptime'] for f in gb])
        tr = np.average([template_window(f, *wr) for f in gr], weights=[f['exptime'] for f in gr])
        err = max(ec, a.floor)
        if mode == 'red':
            bb, br, ebb, ebr = -c / tb, 0., err / tb, 0.
        elif mode == 'blue':
            bb, br, ebb, ebr = 0., c / (k * tr), 0., err / (k * tr)
        else:
            bb = br = c / (k * tr - tb); ebb = ebr = err / abs(k * tr - tb)
        for f in gb:
            beta['blue'][f['base']] = float(bb)
        for f in gr:
            beta['red'][f['base']] = float(br)
        names = '%s + %s' % (', '.join(os.path.basename(f['base']) for f in gb),
                             ', '.join(os.path.basename(f['base']) for f in gr))
        print('  %-20s %s' % (name, names))
        print('  %-20s k = %.3f  c = %+.4f +- %.4f  (%d spaxels)  T_blue %.4f  T_red %.4f  ->  '
              'beta_blue %+.3f +- %.3f  beta_red %+.3f +- %.3f' % ('', k, c, ec, n, tb, tr, bb, ebb, br, ebr))
        res.append(dict(name=name, blue=[f['base'] for f in gb], red=[f['base'] for f in gr], k=float(k), c=float(c),
                        err=float(err), T_blue=float(tb), T_red=float(tr), beta_blue=float(bb), beta_red=float(br),
                        err_blue=float(ebb), err_red=float(ebr), n=int(n)))
    for side in ('blue', 'red'):
        for f in chans[side]['fr']:
            beta[side].setdefault(f['base'], 0.)
    mean = {side: float(np.average([beta[side][f['base']] for f in chans[side]['fr']],
                                   weights=[f['exptime'] for f in chans[side]['fr']])) for side in ('blue', 'red')}
    if len(res) > 1 and mode != 'same':
        side = 'blue' if mode == 'red' else 'red'
        b_ = np.array([r['beta_' + side] for r in res]); e_ = np.array([r['err_' + side] for r in res])
        chi2 = np.sum(((b_ - np.average(b_, weights=1 / e_ ** 2)) / e_) ** 2)
        print('\n  %s beta differs between groups: chi2 = %.1f for %d groups (%s)' % (
            side, chi2, len(res), 'significant: per-group correction needed' if chi2 > 2 * (len(res) - 1)
            else 'consistent with one beta'))
    print('  exposure-weighted mean: beta_blue %+.3f, beta_red %+.3f' % (mean['blue'], mean['red']))
    out = dict(blue=dict(list=a.blue, beta=mean['blue'], frames=beta['blue']),
               red=dict(list=a.red, beta=mean['red'], frames=beta['red']),
               template=a.template, fixed=mode, grouped=True, groups=res)
    fn = a.json or (chans['blue']['base'] + '.skylevel.json')
    json.dump(out, open(fn, 'w'), indent=1)
    print('wrote', fn)


# ---------------------------------------------------------------------------
# measure
# ---------------------------------------------------------------------------
def measure(a):
    chans = {}
    for side, lst in (('blue', a.blue), ('red', a.red)):
        print('%s frames (%s):' % (side, lst))
        frames = channel_template(lst, a.template)
        base, cube, mcube = stack_files(lst)
        newer = [os.path.basename(f['file']) for f in frames if os.path.getmtime(f['file']) > os.path.getmtime(cube)]
        if newer:
            print('WARNING %s is older than %s: the stack was built from earlier versions of the frames '
                  '(re-reduced since?). Restack before measuring.' % (os.path.basename(cube), ', '.join(newer)))
        H = fits.getheader(cube)
        wst = wave(H)
        chans[side] = dict(list=lst, base=base, cube=cube, mcube=mcube, wave=wst, T=stack_template(frames, wst),
                           wg=(H.get('WAVGOOD0', wst[0]), H.get('WAVGOOD1', wst[-1])), header=H,
                           frames=[dict(base=f['base'], src=f['src']) for f in frames], fr=frames)
    B, R = chans['blue'], chans['red']
    Tw = lambda ch, lo, hi: float(np.nanmedian(ch['T'][(ch['wave'] >= lo) & (ch['wave'] <= hi)]))
    rows = []      # (label, coef_blue, coef_red, value, error)

    # dichroic continuity
    wb = (B['wg'][1] - 20 - a.width, B['wg'][1] - 20)
    wr = (R['wg'][0] + 20, R['wg'][0] + 20 + a.width)
    if a.groups:
        return measure_groups(a, chans, wb, wr)
    fb, _ = window_image(B['cube'], B['mcube'], *wb)
    fr, _ = window_image(R['cube'], R['mcube'], *wr)
    if fb.shape != fr.shape:
        sys.exit('blue and red stacks are on different grids')
    k, c, ek, ec, n = robust_line(fr.ravel(), fb.ravel())
    tb, tr = Tw(B, *wb), Tw(R, *wr)
    rows.append(('dichroic %.0f-%.0f / %.0f-%.0f A (k=%.3f, %d spaxels)' % (wb + wr + (k, n)), -tb, k * tr, c, max(ec, a.floor)))

    # image anchors
    kw = WCS(B['header']).celestial
    for spec in a.image or []:
        ref = load_image(spec)
        side = 'blue' if ref['hi'] <= B['wg'][1] + 1 else 'red'
        ch = chans[side]
        lo, hi = max(ref['lo'], ch['wg'][0]), min(ref['hi'], ch['wg'][1])
        kimg, _ = window_image(ch['cube'], ch['mcube'], lo, hi)
        rimg, reg = reference_on_grid(ref['img'], ref['wcs'], kw, kimg.shape, a.fwhm, ref['pix'], kimg=kimg)
        s, b, es, eb, n = robust_line(rimg.ravel(), kimg.ravel())
        t = Tw(ch, lo, hi)
        rows.append(('%s %s %.0f-%.0f A (slope %.3f, shift %+.2f" %+.2f"%s, %d px)' % (
            ref['name'], side, lo, hi, s, reg[0], reg[1], ', r=%.2f' % reg[2] if reg[2] > -2 else '', n),
            -t if side == 'blue' else 0., -t if side == 'red' else 0., b, max(eb, a.floor)))

    # cube-window anchors (cube_windows.py)
    for npz in a.cubewin or []:
        z = np.load(npz)
        rw = WCS(fits.Header.fromstring(str(z['wcs'])))
        pix = np.sqrt(float(z['pixarea']))
        scl = unit_scale(str(z['bunit']))
        for key in [k_ for k_ in z.files if re.match(r'^w\d+_\d+$', k_)]:
            lo, hi = [float(v) for v in key[1:].split('_')]
            side = 'blue' if hi <= B['wg'][1] + 1 else ('red' if lo >= R['wg'][0] - 1 else None)
            if side is None:
                continue
            ch = chans[side]
            kimg, _ = window_image(ch['cube'], ch['mcube'], lo, hi)
            rimg, reg = reference_on_grid(z[key] * scl, rw, kw, kimg.shape, a.fwhm, pix, kimg=kimg)
            s, b, es, eb, n = robust_line(rimg.ravel(), kimg.ravel())
            t = Tw(ch, lo, hi)
            rows.append(('%s %s %.0f-%.0f A (slope %.3f, shift %+.2f" %+.2f"%s, %d px)' % (
                os.path.basename(npz), side, lo, hi, s, reg[0], reg[1], ', r=%.2f' % reg[2] if reg[2] > -2 else '', n),
                -t if side == 'blue' else 0., -t if side == 'red' else 0., b, max(eb, a.floor)))

    # solve
    A = np.array([[r[1], r[2]] for r in rows])
    y = np.array([r[3] for r in rows])
    e = np.array([r[4] for r in rows])
    mode = 'same' if a.same else a.fixed      # 'red', 'blue', 'none' or 'same'
    if mode == 'same':
        M = A.sum(1, keepdims=True)
    elif mode == 'red':         # red kept: beta_red = 0, solve beta_blue
        M = A[:, 0:1]
    elif mode == 'blue':        # blue kept: beta_blue = 0, solve beta_red
        M = A[:, 1:2]
    else:
        M = A
    if M.shape[1] > len(rows):
        sys.exit('two unknowns but only %d constraint(s): add an anchor (--image / --cubewin), '
                 'or use --fixed red|blue or --same' % len(rows))
    W = 1. / e
    p, *_ = np.linalg.lstsq(M * W[:, None], y * W, rcond=None)
    cov = np.linalg.pinv((M * W[:, None]).T @ (M * W[:, None]))
    if mode == 'same':
        bb = br = p[0]; ebb = ebr = np.sqrt(cov[0, 0])
    elif mode == 'red':
        bb, br, ebb, ebr = p[0], 0., np.sqrt(cov[0, 0]), 0.
    elif mode == 'blue':
        bb, br, ebb, ebr = 0., p[0], 0., np.sqrt(cov[0, 0])
    else:
        bb, br = p; ebb, ebr = np.sqrt(np.diag(cov))
    model = A @ np.array([bb, br])
    print('\nconstraints (value = coef_blue x beta_blue + coef_red x beta_red):')
    for r, m in zip(rows, model):
        print('  %-90s %+.4f +- %.4f   model %+.4f   (%+.1f sigma)' % (r[0], r[3], r[4], m, (r[3] - m) / r[4]))
    label = {'same': 'same beta in both channels', 'none': 'both channels free'}.get(mode, 'fixed %s' % mode)
    print('\nsolution (%s): beta_blue = %+.3f +- %.3f, beta_red = %+.3f +- %.3f' % (label, bb, ebb, br, ebr))
    for side, beta in (('blue', bb), ('red', br)):
        ch = chans[side]
        if mode == side:
            print('  %s: fixed, unchanged' % side)
            continue
        print('  %s: adds %+.4f (median, per arcsec2) = %+.0f%% of the subtracted sky continuum' % (
            side, beta * np.nanmedian(ch['T']), 100 * beta))
    if len(rows) > M.shape[1]:
        chi2 = np.sum(((y - model) / e) ** 2)
        print('  chi2 = %.1f for %d constraints: %s' % (chi2, len(rows), 'consistent' if chi2 < 2 * len(rows) else
              'the constraints disagree; check the anchors (sky level of the reference, PSF, field coverage)'))
    out = dict(blue=dict(list=a.blue, beta=float(bb), err=float(ebb)), red=dict(list=a.red, beta=float(br), err=float(ebr)),
               template=a.template, fixed=mode,
               constraints=[dict(label=r[0], coef_blue=float(r[1]), coef_red=float(r[2]), value=float(r[3]),
                                 error=float(r[4]), model=float(m)) for r, m in zip(rows, model)])
    fn = a.json or (B['base'] + '.skylevel.json')
    json.dump(out, open(fn, 'w'), indent=1)
    print('wrote', fn)


# ---------------------------------------------------------------------------
# apply
# ---------------------------------------------------------------------------
def rename_frames(txt, pairs):
    """Rename frame entries (full path or bare name, followed by _icube... / space / end of line)."""
    for old, new in pairs:
        for o, n in ((old, new), (os.path.basename(old), os.path.basename(new))):
            txt = re.sub(r'(?<![\w/.+-])' + re.escape(o) + r'(?=_i|\s|$)', n.replace('\\', r'\\'), txt, flags=re.M)
    return txt


def renamed_list(lst, suffix):
    if lst.endswith('-red.list'):
        return lst[:-len('-red.list')] + suffix + '-red.list'
    return lst[:-len('.list')] + suffix + '.list'


def apply(a):
    sol = json.load(open(a.json or (stack_files(a.blue)[0] + '.skylevel.json')))
    tmpl = sol.get('template', 'sky')
    fixed = sol.get('fixed', {'red': 'blue', 'blue': 'red', 'free': 'none'}.get(sol.get('tie'), sol.get('tie')))
    changed = []
    perframe = {side: sol[side].get('frames') for side in ('blue', 'red')}
    if a.stack and any(perframe.values()):
        sys.exit('per-group betas (measure --groups) differ between frames: they cannot be applied to a stack. '
                 'Run apply without --stack and restack.')
    for side, lst in (('blue', a.blue), ('red', a.red)):
        beta = sol[side]['beta']
        if perframe[side] and any(v != 0 for v in perframe[side].values()):
            beta = None
        if fixed == side or beta == 0:
            print('%s: fixed (beta = 0), nothing written; keep %s and its stack' % (side, lst))
            continue
        changed.append(side)
        print('%s: beta = %s (%s)' % (side, 'per group' if beta is None else '%+.3f' % beta, lst))
        if a.stack:
            base, cube, mcube = stack_files(lst)
            frames = channel_template(lst, tmpl)
            h = fits.open(cube)
            T = stack_template(frames, wave(h[0].header))
            d = h[0].data
            cov = np.isfinite(d) & (d != 0)      # leave no-coverage voxels (0) at 0
            h[0].data = np.where(cov, d + beta * np.nan_to_num(T)[:, None, None], d).astype(d.dtype)
            h[0].header['SKLBETA'] = (beta, 'kcwi_skylevel: + beta x sky continuum')
            nbase = base[:-4] + a.suffix + '-red' if base.endswith('-red') else base + a.suffix
            out = cube.replace(base, nbase, 1)
            h.writeto(out, overwrite=True)
            print('  wrote', out)
            for kind in ('v', 'm', 'e'):
                src = base + '_%scubes.fits' % kind
                dst = nbase + '_%scubes.fits' % kind
                if os.path.isfile(src) and not os.path.exists(dst):
                    try:
                        os.symlink(os.path.abspath(src), dst)
                    except OSError:
                        shutil.copy(src, dst)
            continue
        rows = read_list(lst)
        new = []
        for row in rows:
            lb = local(row[0])
            fn = frame_file(lb)
            w, T, src = cached_template(fn, tmpl)
            beta = perframe[side].get(row[0], 0.) if perframe[side] else sol[side]['beta']
            outdir = a.outdir or os.path.dirname(fn)
            out = os.path.join(outdir, os.path.basename(lb) + a.suffix + fn[len(lb):])
            nb = (os.path.join(a.outdir, os.path.basename(row[0])) if a.outdir else row[0]) + a.suffix
            if not a.dry_run:
                os.makedirs(outdir, exist_ok=True)
                h = fits.open(fn)
                d = h[0].data
                cov = np.isfinite(d) & (d != 0)  # leave DAR padding / no-data voxels (0) at 0
                h[0].data = np.where(cov, d + beta * T[:, None, None], d).astype(d.dtype)
                h[0].header['SKLBETA'] = (beta, 'kcwi_skylevel: + beta x sky continuum')
                h[0].header['SKLTEMPL'] = (src[:60], 'kcwi_skylevel template')
                h[0].header['HISTORY'] = 'kcwi_skylevel: added %+.4f x sky continuum (%s)' % (beta, src)
                h.writeto(out, overwrite=True)
                reg = fn.replace('.fits', '.thum.reg')
                if os.path.isfile(reg):
                    shutil.copy(reg, out.replace('.fits', '.thum.reg'))
            print('  %s -> %s  (beta %+.3f, adds median %+.4g per spaxel)' % (
                os.path.basename(fn), os.path.basename(out), beta, beta * np.nanmedian(T)))
            new.append([nb] + row[1:])
        nl = renamed_list(lst, a.suffix)
        if not a.dry_run:
            with open(nl, 'w') as f:
                for r in new:
                    f.write(' '.join(r) + '\n')
            pairs = [(row[0], r[0]) for row, r in zip(rows, new)]
            for ext in ('.par', '.shift.list', '.preshift.list', '.pre.list', '.flx.list', '.astrom.list'):
                src_ = lst[:-len('.list')] + ext
                if os.path.isfile(src_):
                    txt = open(src_).read()
                    if ext != '.par':
                        txt = rename_frames(txt, pairs)
                    open(nl[:-len('.list')] + ext, 'w').write(txt)
        print('  list %s (+ .par / .shift.list ... copies, frame names updated)' % nl)
    if changed:
        bbase = stack_files(a.blue)[0]
        rbase = stack_files(a.red)[0]
        nb = bbase + a.suffix if 'blue' in changed else bbase
        nr = (rbase[:-4] + a.suffix + '-red' if rbase.endswith('-red') else rbase + a.suffix) if 'red' in changed else rbase
        step = 'kcwi_dichroic_step %s' % nb + ('' if nr == nb + '-red' else ' --red %s' % nr)
        if a.stack:
            print('\ncheck: %s   (c should be ~0)' % step)
        else:
            print('\nnext: stack the new %s (same settings as before), kcwi_xalign wcs, then\n  %s   (c should be ~0)'
                  % (' and '.join(renamed_list(l, a.suffix) for s_, l in (('blue', a.blue), ('red', a.red)) if s_ in changed), step))


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['measure', 'apply'])
    ap.add_argument('blue'); ap.add_argument('red')
    ap.add_argument('--template', default='sky', choices=['sky', 'flat', 'counts'])
    ap.add_argument('--template-step', type=float, default=2.0,
                    help='sample the sky extension every this many A (default 2; 0 = every plane)')
    ap.add_argument('--image', action='append', help='FILE,LO,HI,UNIT[,PIVOT]  (UNIT: hst, nmgy, ps1, flam)')
    ap.add_argument('--cubewin', action='append', help='.npz from cube_windows.py (e.g. MUSE)')
    ap.add_argument('--fixed', choices=['red', 'blue', 'none'], default='red',
                    help='channel kept as it is (default red); none = solve for both (needs an anchor)')
    ap.add_argument('--same', action='store_true', help='one beta for both channels (overrides --fixed)')
    ap.add_argument('--groups', metavar='auto|FILE',
                    help='measure one beta per group of frames (a blue frame + the red frames at the same pointing) '
                         'from the frames themselves: auto = group by pointing (RA, Dec, PA); FILE = one group per line')
    ap.add_argument('--tie', choices=['free', 'same', 'red', 'blue'], help=argparse.SUPPRESS)  # old name
    ap.add_argument('--width', type=float, default=150., help='dichroic windows [A] (default 150)')
    ap.add_argument('--fwhm', type=float, default=1.0, help='reference images smoothed to this [arcsec]')
    ap.add_argument('--floor', type=float, default=0.003, help='minimum error per constraint [stack units]')
    ap.add_argument('--json', help='solution file (default <blue stack>.skylevel.json)')
    ap.add_argument('--suffix', default='_skl')
    ap.add_argument('--outdir', help='apply: write the corrected frames here (default: next to the originals)')
    ap.add_argument('--stack', action='store_true', help='apply: correct the existing stacks instead of the frames')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--path-map', action='append', default=[], metavar='FROM=TO',
                    help='rewrite frame paths of the lists (e.g. /Volumes/X=/mnt/X); the new lists keep FROM')
    return ap


def run(argv):
    """One source, one step (measure or apply). Returns the parsed options."""
    global TEMPLATE_STEP
    a = build_parser().parse_args(argv)
    if a.tie:    # old option: --tie red meant "correct red" (= --fixed blue)
        a.fixed, a.same = {'red': ('blue', False), 'blue': ('red', False), 'free': ('none', False),
                           'same': (a.fixed, True)}[a.tie]
        print('note: --tie is renamed; --tie %s = %s' % (a.tie, '--same' if a.same else '--fixed ' + a.fixed))
    TEMPLATE_STEP = a.template_step
    PATHMAP[:] = []
    for pm in a.path_map:
        f, t = pm.split('=', 1)
        PATHMAP.append((f, t))
    CACHE.clear()
    cache = (stack_files(a.blue)[0] + '.skylevel_templates.json')
    load_cache(cache)
    try:
        measure(a) if a.step == 'measure' else apply(a)
    finally:
        save_cache(cache)
    return a


# ---------------------------------------------------------------------------
# batch: many sources from a file
# ---------------------------------------------------------------------------
BATCH_DOC = """kcwi_skylevel batch SOURCES [--steps measure|apply|measure,apply] [--jobs N] [common options]

SOURCES: one source per line; '#' starts a comment. The first word is the source:
a stack base name (blue list = NAME.list, red list = NAME-red.list), or two list
files. Anything after it is per-source options, which override the common ones
given on the command line (options that can be repeated, like --image, add up):

    # source              per-source options
    IRAS03582+6012
    LEDA66688-final       --groups auto
    UGC2369S              --same
    TMASXJ1735+2047       --fixed blue
    blue/X.list red/X-red.list

Each source writes <NAME>.skylevel.json as usual, and its full output to
<NAME>.skylevel.log; a summary table is printed at the end. A source that fails
is reported and the batch goes on. --jobs N runs N sources at once (separate
processes; worth it when the frames are on different disks or on fast storage).
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


def _source_argvs(row):
    if len(row) >= 2 and row[0].endswith('.list') and row[1].endswith('.list'):
        blue, red, extra = row[0], row[1], row[2:]
    else:
        b = row[0][:-len('.list')] if row[0].endswith('.list') else row[0]
        blue, red, extra = b + '.list', b + '-red.list', row[1:]
    return blue, red, extra


def _summary(a, blue, log):
    fn = a.json or (stack_files(blue)[0] + '.skylevel.json')
    sol = json.load(open(fn))
    mode = sol.get('fixed', '?') + (' groups' if sol.get('grouped') else '')
    if sol.get('grouped'):
        side = 'blue' if sol.get('fixed') == 'red' else 'red'
        bstr = '/'.join('%+.2f' % g['beta_' + side] for g in sol['groups'])
        cstr = '/'.join('%+.3f' % g['c'] for g in sol['groups'])
        bb = bstr if side == 'blue' else '%+.3f' % sol['blue']['beta']
        br = bstr if side == 'red' else '%+.3f' % sol['red']['beta']
    else:
        bb, br = '%+.3f' % sol['blue']['beta'], '%+.3f' % sol['red']['beta']
        d = [c for c in sol.get('constraints', []) if c['label'].startswith('dichroic')]
        cstr = '%+.4f' % d[0]['value'] if d else '-'
    warn = sum(1 for l in log.splitlines() if l.startswith('WARNING'))
    return dict(mode=mode, beta_blue=bb, beta_red=br, c=cstr, warnings=warn)


def _run_source(job):
    """Run the steps for one source with its output going to <name>.skylevel.log (and the terminal)."""
    import contextlib, io, traceback
    name, blue, red, argv_common, extra, steps, echo = job
    buf = io.StringIO()
    out = _Tee(buf, sys.stdout) if echo else buf
    res = dict(source=name, status='ok')
    with contextlib.redirect_stdout(out):
        try:
            print('=' * 100 + '\n%s  (%s, %s)  %s' % (name, blue, red, ' '.join(extra)))
            a = None
            for step in steps:
                a = run([step, blue, red] + argv_common + extra)
            if 'measure' in steps or os.path.isfile(stack_files(blue)[0] + '.skylevel.json'):
                res.update(_summary(a, blue, buf.getvalue()))
        except SystemExit as e:
            res['status'] = 'FAILED: %s' % e
            print('FAILED:', e)
        except Exception as e:
            res['status'] = 'FAILED: %s: %s' % (type(e).__name__, e)
            traceback.print_exc(file=sys.stdout)
    try:
        logname = (blue[:-len('.list')] if blue.endswith('.list') else blue) + '.skylevel.log'
        open(logname, 'w').write(buf.getvalue())
        res['log'] = logname
    except OSError:
        pass
    return res


def batch(argv):
    bp = argparse.ArgumentParser(prog='kcwi_skylevel batch', description=BATCH_DOC,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    bp.add_argument('sources')
    bp.add_argument('--steps', default='measure', help='measure, apply, or measure,apply (default measure)')
    bp.add_argument('--jobs', type=int, default=1, help='sources run at once (default 1)')
    b, common = bp.parse_known_args(argv)
    steps = [x.strip() for x in b.steps.split(',') if x.strip()]
    if not steps or any(x not in ('measure', 'apply') for x in steps):
        sys.exit('--steps: measure, apply or measure,apply')
    if '--json' in common:
        sys.exit('--json would make every source write the same file; give it per source in SOURCES if needed')
    rows = read_list(b.sources)
    jobs = []
    for row in rows:
        blue, red, extra = _source_argvs(row)
        name = os.path.basename(blue[:-len('.list')])
        jobs.append((name, blue, red, common, extra, steps, b.jobs <= 1))
    print('%d sources, steps: %s%s' % (len(jobs), ', '.join(steps), ', %d at once' % b.jobs if b.jobs > 1 else ''))
    if b.jobs <= 1:
        res = [_run_source(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=b.jobs) as ex:
            res = []
            for r in ex.map(_run_source, jobs):
                print('  done: %-28s %s' % (r['source'], r['status']))
                res.append(r)
    print('\n%-28s %-14s %-24s %-24s %-24s %4s  %s' % ('source', 'mode', 'beta_blue', 'beta_red', 'c (dichroic)', 'warn', 'status'))
    for r in res:
        print('%-28s %-14s %-24s %-24s %-24s %4s  %s' % (r['source'], r.get('mode', '-'), r.get('beta_blue', '-'),
              r.get('beta_red', '-'), r.get('c', '-'), r.get('warnings', '-'), r['status']))
    print('logs: <source>.skylevel.log')
    if any(r['status'] != 'ok' for r in res):
        sys.exit(1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'batch':
        return batch(sys.argv[2:])
    run(sys.argv[1:])


if __name__ == '__main__':
    main()
