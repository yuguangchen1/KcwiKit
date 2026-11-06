#!/usr/bin/env python
"""
kcwi_dichroic_step.py -- is there an additive (sky-level) mismatch between the blue
and red KCWI stacks at the dichroic? Model-free: no fit to the spectra.

A galaxy's continuum is smooth across the dichroic gap, so spaxel by spaxel the
flux just blue of the gap is a fixed multiple of the flux just red of it:
    F_blue = k F_red            (k = continuum colour across the gap)
An additive error a_b in blue and a_r in red (sky over/under-subtraction,
scattered light) adds a constant:
    F_blue = k F_red + c,       c = a_b - k a_r
c is per spaxel, so it is invisible where the galaxy is bright and dominates
where it is faint. The script measures F_blue and F_red per spaxel in a window
on each side of the gap (median over wavelength, sky lines and masked voxels
left out), fits the line robustly, and prints F_blue/F_red by brightness.

    python kcwi_dichroic_step.py UGC2369S              # <base>_icubes_wcs.fits and <base>-red_...
    python kcwi_dichroic_step.py TMASXJ1735+2047 --width 150 --plot step.png
    python kcwi_dichroic_step.py IRAS03582+6012_skl --red IRAS03582+6012-red   # red stack unchanged

c = 0 and a flat F_blue/F_red: no additive mismatch. c > 0: blue high / red low
in faint regions. Caveat: a colour gradient that correlates with brightness also
tilts F_blue/F_red, but over a few hundred A it is a few per cent at most.
"""
import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits

SKY = [5577.34, 5889.95, 5895.92, 6300.30, 6363.78]


def wave(h):
    return h['CRVAL3'] + (np.arange(h['NAXIS3']) + 1 - h['CRPIX3']) * h['CD3_3']


def cubefile(base, kind, suffix):
    for s in (suffix, ''):
        fn = '%s_%scubes%s.fits' % (base, kind, s)
        if os.path.isfile(fn):
            return fn
    return None


def window_flux(base, side, width, gapoff, suffix, skyhw=8.):
    fi = cubefile(base, 'i', suffix)
    if fi is None:
        sys.exit('no %s_icubes%s.fits' % (base, suffix))
    hi = fits.open(fi, memmap=True)
    h = hi[0].header
    w = wave(h)
    g0, g1 = h.get('WAVGOOD0', w[0]), h.get('WAVGOOD1', w[-1])
    if side == 'blue':
        lo, hi_ = g1 - gapoff - width, g1 - gapoff
    else:
        lo, hi_ = g0 + gapoff, g0 + gapoff + width
    k = (w >= lo) & (w <= hi_)
    for s in SKY:
        k &= np.abs(w - s) > skyhw
    idx = np.where(k)[0]
    d = np.array(hi[0].data[idx[0]:idx[-1] + 1], float)[k[idx[0]:idx[-1] + 1]]
    fm = cubefile(base, 'm', '')
    if fm is not None:
        m = np.array(fits.open(fm, memmap=True)[0].data[idx[0]:idx[-1] + 1])[k[idx[0]:idx[-1] + 1]]
        d[m != 0] = np.nan
    d[d == 0] = np.nan
    fv = cubefile(base, 'v', '')
    var = None
    if fv is not None:
        var = np.array(fits.open(fv, memmap=True)[0].data[idx[0]:idx[-1] + 1], float)[k[idx[0]:idx[-1] + 1]]
        var[~np.isfinite(d)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        f = np.nanmedian(d, axis=0)
        n = np.isfinite(d).sum(0)
        e = np.sqrt(np.nanmean(var, axis=0) / np.maximum(n, 1)) * 1.2533 if var is not None else None
    return f, e, (lo, hi_), h


def robust_line(x, y, it=10, clip=3.):
    ok = np.isfinite(x) & np.isfinite(y)
    for _ in range(it):
        A = np.vstack([x[ok], np.ones(ok.sum())]).T
        (k, c), *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        r = y - (k * x + c)
        s = 1.4826 * np.nanmedian(np.abs(r[ok]))
        new = np.isfinite(x) & np.isfinite(y) & (np.abs(r) < clip * s)
        if (new == ok).all():
            break
        ok = new
    # bootstrap errors
    rng = np.random.default_rng(1)
    ii = np.where(ok)[0]
    bs = []
    for _ in range(200):
        j = rng.choice(ii, ii.size)
        A = np.vstack([x[j], np.ones(j.size)]).T
        bs.append(np.linalg.lstsq(A, y[j], rcond=None)[0])
    bs = np.array(bs)
    return k, c, bs[:, 0].std(), bs[:, 1].std(), ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('base', help='blue stack base name')
    ap.add_argument('--red', help='red stack base name (default <base>-red)')
    ap.add_argument('--width', type=float, default=150., help='window width each side [A observed] (default 150)')
    ap.add_argument('--gapoff', type=float, default=20., help='stay this far inside WAVGOOD [A] (default 20)')
    ap.add_argument('--suffix', default='_wcs', help="flux cube suffix (default '_wcs', falls back to none)")
    ap.add_argument('--plot')
    a = ap.parse_args()
    fb, eb, wb, hb = window_flux(a.base, 'blue', a.width, a.gapoff, a.suffix)
    fr, er, wr, hr = window_flux(a.red or a.base + '-red', 'red', a.width, a.gapoff, a.suffix)
    if fb.shape != fr.shape:
        sys.exit('blue and red stacks are on different grids %s vs %s' % (fb.shape, fr.shape))
    x, y = fr.ravel(), fb.ravel()
    good = np.isfinite(x) & np.isfinite(y)
    k, c, ek, ec, ok = robust_line(x, y)
    print('%s + %s: %d spaxels with data in both windows; blue %.0f-%.0f A, red %.0f-%.0f A observed' % (
        a.base, a.red or a.base + '-red', good.sum(), wb[0], wb[1], wr[0], wr[1]))
    print('fit  F_blue = k F_red + c :  k = %.3f +- %.3f   c = %+.4g +- %.2g  (%s per spaxel)' % (
        k, ek, c, ec, hb.get('BUNIT', 'cube units')))
    q = np.nanpercentile(x[good], [0, 10, 25, 50, 75, 90, 100])
    print('\n  F_red range            n     median F_red   median F_blue/F_red   c / (k F_red)')
    for lo, hi in zip(q[:-1], q[1:]):
        s = good & (x >= lo) & (x <= hi)
        r = np.nanmedian(y[s] / x[s]) if s.any() else np.nan
        print('  %8.3g - %-8.3g %6d   %10.4g       %8.3f            %+.1f%%' % (
            lo, hi, s.sum(), np.nanmedian(x[s]), r, 100 * c / (k * np.nanmedian(x[s]))))
    sig = abs(c) / ec if ec > 0 else np.inf
    print('\nadditive mismatch c is %.1f sigma from zero; in the faintest tenth it is %+.0f%% of the continuum' % (
        sig, 100 * c / (k * np.nanmedian(x[good & (x <= q[1])]))))
    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        ax[0].plot(x[good], y[good], ',', color='0.6')
        ax[0].plot(x[ok], y[ok], ',', color='C0')
        xx = np.linspace(0, np.nanpercentile(x[good], 99), 10)
        ax[0].plot(xx, k * xx + c, 'C3-', label='k=%.3f, c=%+.3g' % (k, c))
        ax[0].plot(xx, k * xx, 'k--', lw=0.7, label='c = 0')
        ax[0].set_xlim(np.nanpercentile(x[good], [0.5, 99])); ax[0].set_ylim(np.nanpercentile(y[good], [0.5, 99]))
        ax[0].set_xlabel('F red start'); ax[0].set_ylabel('F blue end'); ax[0].legend(fontsize=8)
        ax[1].semilogx(x[good], y[good] / x[good], ',', color='C0')
        ax[1].semilogx(xx[1:], k + c / xx[1:], 'C3-')
        ax[1].axhline(k, color='k', ls='--', lw=0.7)
        ax[1].set_ylim(0, 3 * k); ax[1].set_xlabel('F red start'); ax[1].set_ylabel('F blue / F red')
        fig.suptitle(os.path.basename(a.base)); fig.tight_layout(); fig.savefig(a.plot, dpi=110)
        print('wrote', a.plot)


if __name__ == '__main__':
    main()
