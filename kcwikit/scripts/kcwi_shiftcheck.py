#!/usr/bin/env python
"""
kcwi_shiftcheck -- how wrong is a .shift.list? Compare two of them.

For every shift list given, each frame is placed on the sky exactly as
kcwi_stack would place it (preshift on CRVAL, then CRPIX += shift), and then
checked against the average of the OTHER frames placed the same way:

  * residual offset of each frame [mas, E/N]: what is left to correct. It is
    converted from the measured offset r_k against the others' average with
    e_k = r_k (N-1)/N (the average itself carries -e_k/(N-1)), so it is the
    frame's own error, up to a common offset of the whole channel.
  * R^2 of the frame against that average (1 = identical up to scale/offset/
    gradient): higher means better aligned.
  * with a blue AND a red list: the offset between the blue and red
    averages, i.e. how far apart the two stacks would come out.
  * the difference between the two lists, frame by frame, on the sky.

Nothing is fitted to either list: both are judged by the same yardstick.
The yardstick is the same least-squares image match kcwi_xalign uses, so for
a list made by kcwi_xalign this is a consistency check; for one made by
kcwi_align (or by hand) it is an independent one.

    kcwi_shiftcheck Mrk248.list Mrk248-red.list --a '.shift copy.list' --b .shift.list --z 0.0355
    kcwi_shiftcheck Mrk248.list --a .shift.list.bak_20261023_190500 --b .shift.list

--a / --b replace '.list' in each list name (default --a .shift.list; omit --b
to check one set only). A full path to a file also works with a single list.
"""
import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import ascii

from kcwikit.scripts import kcwi_xalign as X


def shift_file(lst, spec):
    if os.path.isfile(spec):
        return spec
    return lst[:-len('.list')] + spec if lst.endswith('.list') else lst + spec


def load_shifts(fn, nframes):
    t = ascii.read(fn, format='no_header')
    s = np.array([t['col2'], t['col3']], float).T
    if len(s) != nframes:
        raise ValueError('%s has %d rows, the list has %d frames' % (fn, len(s), nframes))
    return s, [os.path.basename(str(n)) for n in t['col1']]


def placed(wf, s):
    w = wf.deepcopy()
    w.wcs.crpix = w.wcs.crpix + np.asarray(s)
    return w


def check_set(g, radius_pix):
    """Per-frame residual (dy, dx) [grid px], its error, R^2, and the channel average."""
    N = len(g)
    r = np.zeros((N, 2)); e = np.zeros((N, 2)); r2 = np.zeros(N)
    infos = []
    for k in range(N):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            tpl = np.nanmean([g[j] for j in range(N) if j != k], axis=0)
        res = X.register(g[k], tpl, radius=radius_pix)
        r[k], e[k], r2[k] = res['d'], res['err'], res['r2']
        infos.append(res)
    if N > 1:
        r = r * (N - 1) / N
        # The channel's common offset cannot be tested here; refer the residuals to
        # the MEDIAN frame, so that one or two bad frames show up as such instead of
        # dragging the reference and making every frame look off.
        r -= np.median(r, axis=0)
    # the channel average as the stack would have it (frames where the shifts put them)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        T = np.nanmean([(g[k] - i['b']) / i['a'] for k, i in enumerate(infos)], axis=0) if N > 1 else g[0]
    return r, e, r2, T


def to_EN(d, pix):
    """grid (dy, dx) -> (E, N) on the sky, mas (grid x points West)."""
    d = np.asarray(d)
    return np.stack([-d[..., 1], d[..., 0]], axis=-1) * pix * 1000.


def sky_corr(wf, s):
    """Sky correction [mas E, N] that adding s to CRPIX applies to the frame."""
    c = wf.wcs.crpix - 1
    a = wf.all_pix2world([[c[0], c[1]]], 0)[0]
    b = wf.all_pix2world([[c[0] - s[0], c[1] - s[1]]], 0)[0]
    return np.array([(b[0] - a[0]) * np.cos(np.deg2rad(a[1])), b[1] - a[1]]) * 3.6e6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('lists', nargs='+', help='one .list, or a blue and a red .list')
    ap.add_argument('--a', default='.shift.list', help="first shift set (suffix replacing '.list', or a file)")
    ap.add_argument('--b', help='second shift set, to compare with the first')
    ap.add_argument('--z', type=float, help='redshift: cut emission lines out of the white-light windows')
    ap.add_argument('--wb', type=float, nargs=2); ap.add_argument('--wr', type=float, nargs=2)
    ap.add_argument('--gridpix', type=float, default=0.1)
    ap.add_argument('--radius', type=float, default=1.5, help='search radius [arcsec] (default 1.5)')
    ap.add_argument('--highpass', type=float, default=0., help='as in kcwi_xalign [arcsec FWHM]')
    ap.add_argument('--tol', type=float, default=50., help='flag residuals above this [mas] (default 50)')
    ap.add_argument('--cubed', action='store_true')
    ap.add_argument('--plot', help='write a summary PNG here')
    a = ap.parse_args()
    if len(a.lists) > 2:
        sys.exit('give one list, or a blue and a red list')

    sets = [('a', a.a)] + ([('b', a.b)] if a.b else [])
    chans = []
    for i, lst in enumerate(a.lists):
        side = 'blue' if i == 0 else 'red'
        base, fn, trim, pres = X.read_list(lst, a.cubed)
        imgs, wcss = [], []
        for f, t, p in zip(fn, trim, pres):
            im, w2, band = X.frame_image(f, a.wb if side == 'blue' else a.wr, a.z, side, trim=t, preshift=p)
            imgs.append(im); wcss.append(w2)
        ch = dict(list=lst, side=side, fn=fn, imgs=imgs, wcs=wcss, band=band, shifts={})
        for key, spec in sets:
            sf = shift_file(lst, spec)
            ch['shifts'][key] = (sf,) + load_shifts(sf, len(fn))
        chans.append(ch)
        print('%s: %s, %d frames, white light %.0f-%.0f A' % (side, lst, len(fn), band[0], band[1]))

    # one grid for everything (all frames, all sets)
    allw = [placed(w, ch['shifts'][k][1][i]) for ch in chans for k, _ in sets for i, w in enumerate(ch['wcs'])]
    alls = [im.shape for ch in chans for k, _ in sets for im in ch['imgs']]
    wg, shape = X.make_grid(allw, alls, chans[0]['wcs'][0].wcs.crval, a.gridpix, a.radius + 1.)
    rad = int(np.ceil(a.radius / a.gridpix))
    hp = a.highpass / a.gridpix / 2.3548 if a.highpass else 0
    pix = a.gridpix

    results = {}
    for ch in chans:
        for key, _ in sets:
            sf, s, names = ch['shifts'][key]
            g = [X.highpass(X.project(im, placed(w, s[i]), wg, shape), hp)
                 for i, (im, w) in enumerate(zip(ch['imgs'], ch['wcs']))]
            r, e, r2, T = check_set(g, rad)
            results[(ch['side'], key)] = dict(r=r, e=e, r2=r2, T=T, file=sf)

    for ch in chans:
        print('\n== %s (%s)' % (ch['side'], ch['list']))
        hdr = '   %-28s' % 'frame'
        for key, _ in sets:
            hdr += ' | set %s: residual E,N [mas]   R^2   ' % key
        if a.b:
            hdr += ' | b - a on the sky [mas]'
        print(hdr)
        for k, f in enumerate(ch['fn']):
            line = '   %-28s' % os.path.basename(f)[:28]
            for key, _ in sets:
                R = results[(ch['side'], key)]
                en = to_EN(R['r'][k], pix)
                flag = ' !' if np.hypot(*en) > a.tol else '  '
                line += ' | %+7.0f %+7.0f%s  %.4f   ' % (en[0], en[1], flag, R['r2'][k])
            if a.b:
                sa, sb = ch['shifts']['a'][1][k], ch['shifts']['b'][1][k]
                dd = sky_corr(ch['wcs'][k], sb) - sky_corr(ch['wcs'][k], sa)
                line += ' | %+7.0f %+7.0f' % tuple(dd)
            print(line)
        for key, _ in sets:
            R = results[(ch['side'], key)]
            en = to_EN(R['r'], pix)
            print('   set %s (%s): rms residual %.0f mas, worst %.0f mas, mean R^2 %.4f' % (
                key, os.path.basename(R['file']), np.sqrt((en ** 2).sum(1).mean()),
                np.hypot(*en.T).max(), R['r2'].mean()))
    if len(chans) == 2:
        print('\n== blue vs red stacks (average of each channel, as placed by the shifts)')
        # stack pixel scale from the blue .par file, when it has one
        stack_pix = X._readpar(a.lists[0].replace('.list', '.par')).get('stack_xpix', -1.)
        for key, _ in sets:
            tr = X.register(results[('red', key)]['T'], results[('blue', key)]['T'], radius=rad)
            en = to_EN(tr['d'], pix)
            d = np.hypot(*en)
            in_px = ' = %.2f stack px of %.3f"' % (d / (stack_pix * 1000.), stack_pix) if stack_pix > 0 else ''
            print('   set %s: red sits %+.0f mas E, %+.0f mas N of blue (|d| = %.0f mas%s)%s' % (
                key, en[0], en[1], d, in_px, '  !' if d > a.tol else ''))
    print('\nResiduals are relative to each channel\'s median frame; "!" marks > %.0f mas.' % a.tol)
    print('A residual of a few tens of mas is the noise floor for faint frames; compare the sets, not')
    print('just the numbers: the better list has the smaller residuals. R^2 changes little at high S/N.')

    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, len(chans), figsize=(5.2 * len(chans), 4.8), squeeze=False)
        for c_, ch in enumerate(chans):
            for key, mk in zip([k for k, _ in sets], ['o', 's']):
                en = to_EN(results[(ch['side'], key)]['r'], pix)
                ax[0, c_].plot(en[:, 0], en[:, 1], mk, mfc='none' if key == 'a' else None,
                               label='set %s: %s' % (key, os.path.basename(results[(ch['side'], key)]['file'])))
                for i in range(len(en)):
                    ax[0, c_].annotate(str(i), en[i], fontsize=7, xytext=(3, 3), textcoords='offset points')
            t = np.linspace(0, 2 * np.pi, 100)
            ax[0, c_].plot(a.tol * np.cos(t), a.tol * np.sin(t), color='0.6', lw=0.7)
            ax[0, c_].axhline(0, color='0.85', lw=0.5); ax[0, c_].axvline(0, color='0.85', lw=0.5)
            ax[0, c_].invert_xaxis(); ax[0, c_].set_aspect('equal', 'datalim')
            ax[0, c_].set_xlabel('residual E [mas]'); ax[0, c_].set_ylabel('residual N [mas]')
            ax[0, c_].set_title('%s frames' % ch['side'], fontsize=10); ax[0, c_].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(a.plot, dpi=120)
        print('wrote', a.plot)


if __name__ == '__main__':
    main()
