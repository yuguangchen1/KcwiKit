#!/usr/bin/env python
"""
kcwi_cube_windows -- collapse a (big) IFU cube into median images in a few
wavelength windows, so the images can be compared elsewhere without moving the cube.

Run it where the cube lives (e.g. on bifrost, not over a network mount):

    kcwi_cube_windows UGC_2369_S_DATACUBE_FINAL.fits ugc2369s_muse_windows.npz \
        --win 4800 4950 --win 5408 5558 --win 5757 5907 --win 7020 7170

Writes an .npz with one image per window (median over the window; the sky lines
5577, 5890/5896, 6300, 6364, 7316, 7341 A +-8 A left out), in the cube's flux unit
per arcsec^2, plus the celestial WCS header as a string. Reads the cube in chunks of
planes, so memory stays at a few hundred MB even for a 3.5 GB MUSE cube.
"""
import argparse
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

SKY = [5577.34, 5889.95, 5895.92, 6300.30, 6363.78, 7316.3, 7340.9]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cube')
    ap.add_argument('out')
    ap.add_argument('--win', type=float, nargs=2, action='append', required=True, metavar=('LO', 'HI'))
    ap.add_argument('--ext', default=None, help='data extension (default: DATA if present, else 0)')
    a = ap.parse_args()
    warnings.simplefilter('ignore')
    h = fits.open(a.cube, memmap=True)
    ext = a.ext if a.ext is not None else ('DATA' if 'DATA' in h else 0)
    H = h[ext].header
    step = H.get('CD3_3', H.get('CDELT3'))
    w = H['CRVAL3'] + (np.arange(H['NAXIS3']) + 1 - H['CRPIX3']) * step
    wcs = WCS(H).celestial
    pixarea = abs(np.linalg.det(wcs.pixel_scale_matrix)) * 3600. ** 2
    out = {'wcs': wcs.to_header_string(), 'bunit': str(H.get('BUNIT', '')) + ' per arcsec2', 'pixarea': pixarea}
    for lo, hi in a.win:
        k = (w >= lo) & (w <= hi)
        for s in SKY:
            k &= np.abs(w - s) > 8
        idx = np.where(k)[0]
        if idx.size == 0:
            print('window %g-%g: outside the cube, skipped' % (lo, hi))
            continue
        d = np.array(h[ext].data[idx[0]:idx[-1] + 1], np.float32)[k[idx[0]:idx[-1] + 1]]
        img = np.nanmedian(d, axis=0) / pixarea
        out['w%d_%d' % (lo, hi)] = img
        print('window %g-%g: %d planes, median %.4g' % (lo, hi, idx.size, np.nanmedian(img)))
    np.savez_compressed(a.out, **out)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
