"""
DS9 region files -> boolean pixel masks, with the semantics the old
``pyregion`` code gave us:

* region files in sky coordinates (fk5, icrs, galactic, ...) are converted to
  pixels with the WCS of the supplied header;
* ``physical`` coordinates are treated as ``image`` (KCWI frames are never
  binned relative to the DS9 display);
* shapes with a leading ``-`` are *excluded* (subtracted after the included
  shapes are OR-ed together);
* a shape that lies entirely outside the image contributes nothing instead of
  raising.

The astropy ``regions`` package is used when available; ``pyregion`` is the
fallback so that either dependency works.

    >>> from kcwikit.kcwi.kcwi_regions import regions_to_mask
    >>> mask = regions_to_mask('frame.reg', header2d, shape=(ny, nx))
    >>> img[mask] = np.nan

``shape`` is (ny, nx), i.e. the numpy shape of the 2-D image, and the returned
mask has the same shape.
"""

import os
import warnings

import numpy as np
from astropy.wcs import WCS

__all__ = ['read_region_string', 'regions_to_mask', 'header_to_2d']


def header_to_2d(header):
    """
    Return a copy of ``header`` reduced to its two spatial axes, which is what
    both ``regions`` and ``pyregion`` need to interpret a region file.
    A 2-D header is returned unchanged (as a copy).
    """
    hdr = header.copy()
    if hdr.get('NAXIS', 2) <= 2:
        return hdr
    for key in list(hdr.keys()):
        if key.endswith('3') and key[:-1] in ('NAXIS', 'CRVAL', 'CRPIX', 'CDELT', 'CTYPE',
                                             'CUNIT', 'CNAME', 'CROTA'):
            del hdr[key]
        elif key in ('CD3_3', 'CD1_3', 'CD2_3', 'CD3_1', 'CD3_2',
                     'PC3_3', 'PC1_3', 'PC2_3', 'PC3_1', 'PC3_2'):
            del hdr[key]
    hdr['NAXIS'] = 2
    hdr['WCSAXES'] = 2
    return hdr


def read_region_string(regfn):
    """Read a region file and normalise it (physical -> image)."""
    with open(regfn, 'r') as f:
        regstr = f.read()
    if 'physical' in regstr:
        regstr = regstr.replace('physical', 'image')
    return regstr


def _mask_with_regions(regstr, header, shape):
    from regions import Regions, SkyRegion

    regs = Regions.parse(regstr, format='ds9')
    include = np.zeros(shape, dtype=bool)
    exclude = np.zeros(shape, dtype=bool)
    wcs2d = None
    for reg in regs:
        if isinstance(reg, SkyRegion):
            if header is None:
                raise ValueError('Region file is in sky coordinates but no header was given.')
            if wcs2d is None:
                wcs2d = WCS(header_to_2d(header)).celestial
            reg = reg.to_pixel(wcs2d)
        try:
            img = reg.to_mask(mode='center').to_image(shape)
        except Exception as e:  # e.g. an unsupported shape (text, ruler, ...)
            warnings.warn('Skipping region %s: %s' % (reg.__class__.__name__, e))
            continue
        if img is None:        # shape entirely outside the image
            continue
        m = img.astype(bool)
        if reg.meta.get('include', 1) in (0, False):
            exclude |= m
        else:
            include |= m
    return include & ~exclude


def _mask_with_pyregion(regstr, header, shape):
    import pyregion
    regs = pyregion.parse(regstr)
    if header is not None:
        regs = regs.as_imagecoord(header_to_2d(header))
    return regs.get_mask(shape=shape).astype(bool)


def regions_to_mask(regfn, header=None, shape=None, prefer='regions'):
    """
    Rasterise a DS9 region file into a boolean mask.

    Args:
        regfn (str): region file name.
        header (astropy.io.fits.Header): header of the image the regions refer
            to (2-D or 3-D; only the spatial axes are used). Required for sky
            coordinates; used for ``shape`` when ``shape`` is None.
        shape (tuple): (ny, nx) of the output mask.
        prefer (str): 'regions' (default) or 'pyregion'; the other one is used
            as a fallback when the preferred package is not installed.

    Returns:
        numpy.ndarray(bool): True inside the (included minus excluded) shapes.
    """
    if shape is None:
        if header is None:
            raise ValueError('Give either shape or a header with NAXIS1/NAXIS2.')
        shape = (int(header['NAXIS2']), int(header['NAXIS1']))
    shape = tuple(int(s) for s in shape)

    if not os.path.isfile(regfn):
        raise FileNotFoundError(regfn)
    regstr = read_region_string(regfn)

    backends = [_mask_with_regions, _mask_with_pyregion]
    if prefer == 'pyregion':
        backends = backends[::-1]
    last_err = None
    for fn in backends:
        try:
            return fn(regstr, header, shape)
        except ImportError as e:
            last_err = e
            continue
    raise ImportError('Neither the "regions" nor the "pyregion" package is installed '
                      '(%s)' % last_err)
