#!/usr/bin/env python

"""
kcwi_automask.py -- make KCWI/KCRM masks automatically instead of drawing them
by hand in DS9. Many modes:

INTF MODE (sky mask for the DRP, <base>_intf.fits -> <base>.reg -> <base>_smsk.fits)
    In an _intf frame every slice is a vertical strip with wavelength along y, so
    a continuum source is a narrow column that stands above the rest of its slice
    at (almost) every wavelength, while sky lines run along x and cosmic rays are
    compact. Slice by slice:
        1. slices are found from the blank gaps between the strips;
        2. the sky spectrum (median over the slice's columns at each row) is
           subtracted, removing the sky lines;
        3. the median along wavelength of that residual gives one number per
           column, immune to cosmic rays and sky-line residuals;
        4. a running median across columns gives the slice background, and the
           scatter of what is left sets the threshold: a column above nsig is a
           detection, and the flag then follows its neighbours down to grow_nsig.
    Slices are trimmed by a couple of columns first, because the partly
    illuminated columns beside each gap carry offsets that mimic sources.
    The mask is then followed as a function of wavelength rather than being a
    plain rectangle: a bright object scatters light along the slice in proportion
    to its count rate, so the masked width grows where the object is bright and
    shrinks where it is faint (--rect restores the old rectangles). Light that
    leaks into neighbouring slices is picked up as its own trace there.
    Only continuum traces are found. A source seen solely in an emission line
    leaves the column median unchanged, so keep drawing those by hand.

THUM MODE (mask for median filtering, <base>_icube.thum.fits ->
           <base>_icube.thum.reg -> <base>_icube.mask.fits)
    The thumbnail is a white-light image with one column per slice. Each slice's
    own level is removed (the striping), sources are detected above nsig * MAD,
    grown by a small amount only, and written as ellipses. Two constraints keep
    the median filter well behaved: detections touching the edge of the field are
    dropped (--edge), and no slice is allowed to be masked by more than --maxfrac
    of its length, with the weakest sources dropped first when it would be.

BCRMSK MODE (astroscrappy protection mask for the DRP, <base>_intf.fits ->
             <base>_bcrmsk.reg -> <base>_intf_bcrmsk.fits)
    KCWI_DRP's RemoveCosmicRays hands <base>_intf_bcrmsk.fits to astroscrappy as
    `inmask`, the pixels detect_cosmics must ignore. A bright quasar is exactly
    what its Laplacian filter is built to find -- narrow and steep -- so parts
    of the object get "cleaned" away. Measured against the DRP's own flags on
    eight 600 s frames, astroscrappy's rate rises from 0.08% of pixels at sky
    level to 1.4% above 2500 counts, a rate that should not depend on
    brightness at all: about 150 object pixels destroyed per frame.
    The sky mask has to over-cover, because source flux left in it biases the
    sky model. This one has to under-cover: inside it a real cosmic ray now
    survives, to be caught later by crr_split in kcwi_stack or by a crmsk. So
    it masks the bright spots themselves rather than whole traces -- the sky
    spectrum is subtracted per slice, a short running median along wavelength
    (--bcr-zmed) removes anything narrower than the resolution, so a cosmic ray
    the DRP missed on the first pass cannot be protected on the second, and
    what is left above --bcr-thresh counts is masked. On those eight frames
    that covered 80% of the false flags with 0.83% of the frame, against 76%
    and 2.51% for full-height strips (--bcr-mode strips), and gave up three
    times fewer real cosmic rays -- better on both counts. Needs an _intf
    frame, so this is a second-pass product: reduce once, build the masks, set
    bcrmsk in kcwi_blue.cfg, reduce again.

WL MODE (source mask for KSkyWizard/ZAP, <base>_wlimg.fits -> <base>.reg
         -> <base>_zap_mask.fits)
    The white-light image KSkyWizard writes beside each frame, detected the same
    way as a thumbnail but with three differences. The region file goes out in
    physical coordinates as ellipses, which is what the GUI reads back. There is
    no per-slice budget (--maxfrac defaults to 1): ZAP estimates the background
    per pixel with a PCA, so unlike the median filter it does not need part of
    each slice left unmasked. And the usable area comes from KSkyWizard's own
    flag criterion when that is informative -- on a frame whose flags are raised
    at every wavelength it covers the whole image, and detection then falls back
    to all finite pixels (--min-usable) -- after which the illumination ramp at
    the two ends of the slice is cut by level (--trim-nsig), a row at a time,
    while its median runs away from the interior. That strip reaches thousands of
    counts against an interior of ~100 and would otherwise register as a long
    source lying along the boundary. Because the ramp is removed by level, --edge
    defaults to 0 here rather than 2; a blind erosion on top of the trim was
    throwing away real sources a row or two inside the field, at the cost of the
    occasional two-pixel blip in a corner, which --edge 1 suppresses.

All four modes write a DS9 region file, following the naming
that kcwi_makemask / kcwi_makemask_medfilter expect, so you can check and edit
in DS9 and then convert as usual:

    kcwi_automask.py kb261007_*_intf.fits
    ds9 kb261007_00074_intf.fits -regions kb261007_00074.reg
    kcwi_makemask.py .

    kcwi_automask.py --json kcwi_2026oct07_red.json --thum
    kcwi_makemask_medfilter.py . --red -f kcwi_2026oct07_red.json

Options:
    --json FILE   take the frames from a grouping JSON ({group: [files]}), the
                  first file of each group, as the open_ds9 scripts do; --all
                  processes every file instead.
    --wl          white-light mode: write the KSkyWizard region file <base>.reg
                  from <base>_wlimg.fits (ellipses, physical coordinates).
    --bcrmsk      astroscrappy protection mask: <base>_intf_bcrmsk.fits (and
                  <base>_bcrmsk.reg) from <base>_intf.fits. The FITS is written
                  by default here, because it is what the DRP reads. Tune the
                  cosmic-ray cost with --bcr-thresh: higher masks less.
    --convert     do not detect anything: build the mask FITS from the region
                  file that is already there. This is the mode to use after
                  hand-editing a .reg in DS9 -- the region file is only read.
                  With --json --all the group's first frame supplies the regions
                  for every frame of that group, as kcwi_makemask_medfilter does.
    --no-clobber  in detection mode, skip any frame that already has a region
                  file, so hand edits are never overwritten.
    --plot        write a quick-look PNG per frame showing what was detected.
    --smsk        also write the mask FITS directly (_smsk.fits / _icube.mask.fits).

Examples:
    kcwi_automask.py kb261007_00086_intf.fits --convert     # after editing the .reg
    kcwi_automask.py --json red.json --thum --all --convert # one .reg per group
    kcwi_automask.py kb*_intf.fits --plot
    kcwi_automask.py --json red.json --bcrmsk --plot    # then re-run the DRP
    kcwi_automask.py kb261007_00086_intf.fits --bcrmsk --convert
    kcwi_automask.py --json blue.json --plot --smsk
    kcwi_automask.py --json red.json --thum --grow 1 --maxfrac 0.25 --plot
"""

import os
import re
import json
import argparse
import warnings

import numpy as np
from astropy.io import fits
from scipy import ndimage


# --------------------------------------------------------------------------
# generic helpers
# --------------------------------------------------------------------------

def _mad(x):
    """Median absolute deviation, scaled to a Gaussian sigma."""
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def _robust_sigma(r, outer=False):
    """Robust scatter of an array.

    Takes the larger of the MAD and half the 16-84 interpercentile range.

    With outer=True the 90th percentile of |r|, scaled to a Gaussian, joins
    them. That is for the _intf column profiles only: there the data are
    quantized and the sky spectrum is a median over the slice's own columns, so
    a typical column cancels exactly and more than half the values can be
    exactly zero, which collapses any centre-based estimator (one slice was
    reporting sigma seven times too small and masking 20 columns of noise).
    The 2D images have no such pile-up, and there the outer percentile is
    inflated by the wings of real sources instead -- on a white-light image it
    came out 3x the MAD and hid everything faint -- so they use the default.
    Falls back to the scatter of the first differences and to a clipped
    standard deviation if both are degenerate.
    """
    r = np.asarray(r)
    r = r[np.isfinite(r)]
    if r.size < 3:
        return np.nan
    p16, p84 = np.percentile(r, [16, 84])
    cands = [_mad(r),                                   # centre-based
             0.5 * (p84 - p16)]                         # interpercentile
    if outer:
        cands.append(np.percentile(np.abs(r), 90) / 1.6449)
    s = max([c for c in cands if np.isfinite(c)] or [0.])
    if s > 0:
        return s
    s = _mad(np.diff(r)) / np.sqrt(2.)
    if np.isfinite(s) and s > 0:
        return s
    sd = np.std(r)
    keep = np.abs(r - np.median(r)) < 3 * sd if sd > 0 else np.ones(r.size, bool)
    s = np.std(r[keep])
    return s if np.isfinite(s) and s > 0 else np.nan


def _runs(flag):
    """Contiguous True runs of a 1D boolean array, as (start, stop) pairs."""
    out = []
    i, n = 0, len(flag)
    while i < n:
        if flag[i]:
            j = i
            while j < n and flag[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def _running_median(y, size):
    return ndimage.median_filter(y, size=max(3, int(size) | 1), mode='nearest')


def _write_regions(regfn, lines, srcname, coord='image'):
    with open(regfn, 'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write('# generated by kcwi_automask.py from {}\n'.format(srcname))
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" '
                'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write(coord + '\n')
        for ln in lines:
            f.write(ln + '\n')


def _trim_edge_ramp(data, cover, nsig=5., maxtrim=5):
    """Shrink `cover` at the ends of each axis while the boundary level runs away.

    The illumination correction breaks down in the outermost rows of the field,
    so a white-light image often ends in a strip of very large positive or
    negative values. Those are not sources, but they are bright enough to
    detect: on one frame the last three rows along the slice ran to 1300-3700
    against an interior level of 90. KSkyWizard's flag criterion normally keeps
    them out of the usable area, and a fixed erosion covers them only if it
    happens to be deep enough, so this cuts them by level instead. Starting at
    each end, a row is dropped while its median sits more than `nsig` interior
    scatters off the interior median, and the walk stops at the first row that
    does not.

    Args:
        data (ndarray): the image.
        cover (bool ndarray): the usable area.
        nsig (float): how far off the interior level a boundary row has to be.
        maxtrim (int): most rows dropped at one end of one axis.

    Returns:
        (cover, trimmed) with `trimmed` the number of rows (or columns) cut at
        each end of each axis, as ((lo, hi), (lo, hi)) for axes 0 and 1.
    """
    out = cover.copy()
    trimmed = []
    for axis in (0, 1):
        idx = np.where(out.any(axis=1 - axis))[0]
        if idx.size < 8:
            trimmed.append((0, 0))
            continue
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            masked = np.where(out, data, np.nan)
            prof = np.nanmedian(masked, axis=1 - axis)[idx]
        # the interior half sets the level and the scatter to compare against,
        # so a source sitting at one end cannot define its own reference
        q = idx.size // 4
        ref_p = prof[q:idx.size - q]
        ref_p = ref_p[np.isfinite(ref_p)]
        if ref_p.size < 4:
            trimmed.append((0, 0))
            continue
        ref = np.median(ref_p)
        sig = _robust_sigma(ref_p - ref)
        if not np.isfinite(sig) or sig <= 0:
            trimmed.append((0, 0))
            continue
        cut = [0, 0]
        for end in (0, 1):
            order = range(idx.size) if end == 0 else range(idx.size - 1, -1, -1)
            for k in order:
                if cut[end] >= maxtrim:
                    break
                if np.isfinite(prof[k]) and abs(prof[k] - ref) <= nsig * sig:
                    break
                if axis == 0:
                    out[idx[k], :] = False
                else:
                    out[:, idx[k]] = False
                cut[end] += 1
        trimmed.append((cut[0], cut[1]))
    return out, tuple(trimmed)



# --------------------------------------------------------------------------
# intf mode
# --------------------------------------------------------------------------

def find_slices(data, min_width=20, min_level=0.5, edge_trim=2):
    """Locate the slices of an _intf frame from the illumination level.

    A gap between slices has a column median of zero, and the one or two columns
    on either side of it are only partly illuminated. Those ramp columns carry
    large offsets that look like sources, so slices are cut where the column
    median falls below `min_level` times the typical slice level, and then
    trimmed by `edge_trim` columns on each side.

    Args:
        data (ndarray): the _intf image, [ny, nx].
        min_width (int): ignore runs narrower than this.
        min_level (float): threshold on the column median, as a fraction of the
            typical slice level.
        edge_trim (int): columns dropped at each end of a slice.

    Returns:
        list of (x0, x1): slice column ranges, x1 exclusive, already trimmed.
    """
    lev = np.median(np.where(np.isfinite(data), data, 0.), axis=0)
    pos = lev[lev > 0]
    if pos.size == 0:
        return []
    ref = np.median(pos)
    runs = [(a, b) for a, b in _runs(lev > min_level * ref) if b - a >= min_width]
    if not runs:
        return []
    # a gap that is not perfectly blank can merge neighbouring slices: split any
    # run that is a multiple of the typical width
    typical = np.median([b - a for a, b in runs])
    out = []
    for a, b in runs:
        n = int(round((b - a) / typical)) if typical > 0 else 1
        if n > 1:
            edges = np.linspace(a, b, n + 1).astype(int)
            out += list(zip(edges[:-1], edges[1:]))
        else:
            out.append((a, b))
    return [(a + edge_trim, b - edge_trim) for a, b in out if b - a > 2 * edge_trim + 4]


def detect_slice(sub, nsig=8., grow_nsig=3., grow=1, bkg_window=11., min_frac_valid=0.2):
    """Flag the columns of one slice that carry an object trace.

    A column is a detection when it rises above nsig, and the flag then follows
    the neighbouring columns down to grow_nsig, so the wings of the trace are
    covered without widening isolated noise spikes.

    Args:
        sub (ndarray): the slice's pixels, [ny, nw].
        nsig (float): threshold for the core of a detection.
        grow_nsig (float): the flag extends over connected columns above this.
        grow (int): extra columns added on each side afterwards.
        bkg_window (float): running-median window for the background, in columns.
            It should be a few times the source width; 0 uses half the slice width.
        min_frac_valid (float): columns with less valid data than this are skipped.

    Returns:
        (flag, snr): boolean array [nw], and the profile significance [nw].
    """
    nw = sub.shape[1]
    good = np.isfinite(sub) & (sub != 0)
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)              # all-NaN rows
        sky = np.nanmedian(np.where(good, sub, np.nan), axis=1)      # sky spectrum
        resid = np.where(good, sub - sky[:, None], np.nan)
        prof = np.nanmedian(resid, axis=0)                           # one value per column

    usable = good.mean(axis=0) > min_frac_valid
    prof = np.where(usable, prof, np.nan)
    if np.isfinite(prof).sum() < 5:
        return np.zeros(nw, dtype=bool), np.zeros(nw)

    window = bkg_window if bkg_window > 0 else max(5, nw // 2)
    filled = np.where(np.isfinite(prof), prof, np.nanmedian(prof))
    bkg = _running_median(filled, window)
    sig = _robust_sigma(prof - bkg, outer=True)
    if not np.isfinite(sig) or sig <= 0:
        return np.zeros(nw, dtype=bool), np.zeros(nw)

    snr = np.where(np.isfinite(prof), (prof - bkg) / sig, 0.)
    core = snr > nsig
    flag = np.zeros(nw, dtype=bool)
    if core.any():
        for a, b in _runs(snr > grow_nsig):
            if core[a:b].any():
                flag[a:b] = True
        if grow > 0:
            flag = ndimage.binary_dilation(flag, np.ones(2 * int(grow) + 1, dtype=bool))
    return flag & usable, snr



def trace_boxes(sub, cores, x0, ny_full, ybin=100, wing_nsig=3., seed_nsig=5.,
                maxfrac=0.5, far_pad=10):
    """Follow a trace's width as a function of wavelength.

    A bright QSO scatters light along the slice, and the amount scattered tracks
    the count rate, which varies strongly with wavelength (in one slice of a real
    frame the residual at the trace runs from 14 to 1244 counts between the blue
    and red ends). A single full-height rectangle is therefore too narrow where
    the object is bright and too wide where it is faint.

    The slice is divided into bins of `ybin` rows. In each bin the trace's wings
    are measured against the slice's own baseline -- taken from the columns more
    than `far_pad` away from any core, since a running median would simply follow
    a broad wing and subtract it -- and the mask extends outwards while the
    profile stays above wing_nsig. Bins that end up with the same extent are
    merged back together, so a faint trace still produces a single box.

    Args:
        sub (ndarray): the slice's pixels, [ny, nw].
        cores (list of (a, b, is_seed)): column ranges to start from. A seed
            (propagated from a bright neighbouring slice) only produces a mask in
            the bins where it carries signal above seed_nsig; a detection always does.
        x0 (int): the slice's first column in the frame.
        ny_full (int): rows in the frame.
        ybin (int): rows per wavelength bin.
        wing_nsig (float): threshold for extending the mask sideways.
        seed_nsig (float): a seeded position must reach this in a bin to be masked
            there; higher than wing_nsig because a seed has no detection of its own.
        maxfrac (float): most of the slice that may be masked in any one bin.

    Returns:
        list of (bx0, bx1, by0, by1): boxes in frame coordinates, stop-exclusive.
    """
    ny, nw = sub.shape
    if not cores:
        return []
    # merge overlapping entries: several bright slices can seed the same place,
    # and a seed touching a detection belongs to that detection
    merged = []
    for (a, b, seed) in sorted(cores):
        if merged and a <= merged[-1][1]:
            pa, pb, pseed = merged[-1]
            merged[-1] = (pa, max(pb, b), pseed and seed)
        else:
            merged.append((a, b, seed))
    cores = merged
    good = np.isfinite(sub) & (sub != 0)
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        sky = np.nanmedian(np.where(good, sub, np.nan), axis=1)
        resid = np.where(good, sub - sky[:, None], np.nan)

    far = np.ones(nw, dtype=bool)
    for (ca, cb, _) in cores:
        far[max(ca - far_pad, 0):cb + far_pad] = False
    if far.sum() < 5:
        far = np.ones(nw, dtype=bool)

    nb = max(1, int(np.ceil(ny / float(ybin))))
    edges = [(int(i * ny / nb), int((i + 1) * ny / nb)) for i in range(nb)]
    colmax = max(1, int(np.floor(maxfrac * nw)))

    spans = []                                   # [bin][core] = (lo, hi) or None
    for (y0, y1) in edges:
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            prof = np.nanmedian(resid[y0:y1], axis=0)
        ok = np.isfinite(prof)
        row = []
        if ok.sum() < 5:
            spans.append([(a, b) if not seed else None for (a, b, seed) in cores])
            continue
        base = np.median(prof[far & ok])
        sig = _robust_sigma(prof[far & ok] - base, outer=True)
        snr = np.where(ok, (prof - base) / sig, 0.) if (np.isfinite(sig) and sig > 0) \
            else np.zeros(nw)
        for (ca, cb, seed) in cores:
            if seed and np.nanmax(snr[ca:cb]) <= seed_nsig:
                row.append(None)                 # no signal here: mask nothing
                continue
            lo, hi = ca, cb
            while lo > 0 and snr[lo - 1] > wing_nsig and (hi - lo) < colmax:
                lo -= 1
            while hi < nw and snr[hi] > wing_nsig and (hi - lo) < colmax:
                hi += 1
            row.append((lo, hi))
        spans.append(row)

    # a one-column flicker between neighbouring bins would split the mask into
    # many boxes, so widen each bin to the envelope of itself and its neighbours
    smoothed = []
    for ib in range(nb):
        row = []
        for ic in range(len(cores)):
            near = [spans[j][ic] for j in (ib - 1, ib, ib + 1)
                    if 0 <= j < nb and spans[j][ic] is not None]
            if spans[ib][ic] is None and not near:
                row.append(None)
            elif spans[ib][ic] is None:
                row.append(None)
            else:
                row.append((min(p[0] for p in near), max(p[1] for p in near)))
        smoothed.append(row)
    spans = smoothed

    # merge consecutive bins with the same extent
    boxes = []
    for ic in range(len(cores)):
        b_start = 0
        for ib in range(1, nb + 1):
            same = ib < nb and spans[ib][ic] == spans[b_start][ic]
            if not same:
                span = spans[b_start][ic]
                if span is not None:
                    lo, hi = span
                    y0 = edges[b_start][0]
                    y1 = edges[ib - 1][1]
                    boxes.append((x0 + lo, x0 + hi, y0, y1))
                b_start = ib
    return boxes


def automask_intf(imfname, nsig=8., grow_nsig=3., grow=1, maxfrac=0.5, bkg_window=11.,
                  edge_trim=2, adaptive=True, ybin=200, wing_nsig=3., seed_nsig=5., bright_nsig=50.,
                  bright_reach=3, regfn=None, write_mask=False, mask_out=None, plot=False,
                  verbose=True):
    """Detect object traces in an _intf frame and write the sky mask.

    Args:
        imfname (str): the _intf.fits frame.
        nsig (float): threshold for the core of a detection, in units of the
            robust scatter of the slice's column profile.
        grow_nsig (float): the flag follows connected columns down to this level.
        grow (int): extra columns added on each side afterwards.
        maxfrac (float): never mask more than this fraction of a slice's columns;
            the threshold is raised for that slice until the limit is met, so the
            sky model always keeps enough sky.
        bkg_window (float): running-median window for the slice background [columns].
        edge_trim (int): columns dropped at each end of a slice, where the
            illumination ramps down and the profile is unreliable.
        adaptive (bool): let the masked width follow the brightness as a function
            of wavelength (see trace_boxes). False gives full-height rectangles.
        ybin (int): rows per wavelength bin when adaptive.
        wing_nsig (float): threshold for extending the mask sideways [sigma].
        seed_nsig (float): threshold a seeded position must reach to be masked.
        bright_nsig (float): a trace above this is bright enough to scatter light
            into neighbouring slices, so its position is seeded there as well.
        bright_reach (int): how many slices away a seed is carried. A seed only
            masks where it carries signal above wing_nsig.
        regfn (str): output region file. Default <base>.reg, which is what
            kcwi_makemask expects next to <base>_intf.fits.
        write_mask (bool): also write the mask FITS (uint8, 1 = masked).
        mask_out (str): where that mask goes. Default <base>_smsk.fits.
        plot (bool): write a quick-look PNG.

    Returns:
        (regfn, nboxes, frac)
    """
    data = fits.getdata(imfname)
    if data.ndim != 2:
        raise ValueError('{} is not a 2D _intf image'.format(imfname))
    ny, nx = data.shape
    slices = find_slices(data, edge_trim=edge_trim)
    if verbose:
        print('{}: {} x {}, {} slices'.format(os.path.basename(imfname), nx, ny, len(slices)))

    mask = np.zeros((ny, nx), dtype=bool)
    boxes = []                                     # (x0, x1, y0, y1), stop-exclusive
    snr_all = np.zeros(nx)

    # pass 1: detect the traces slice by slice
    flags, peaks = [], []
    for isl, (x0, x1) in enumerate(slices):
        sub = data[:, x0:x1]
        thresh = nsig
        flag = np.zeros(x1 - x0, dtype=bool)
        for _ in range(12):                        # raise the threshold if too much is flagged
            flag, snr = detect_slice(sub, nsig=thresh, grow_nsig=grow_nsig, grow=grow,
                                     bkg_window=bkg_window)
            if flag.mean() <= maxfrac:
                break
            thresh *= 1.4
        else:
            flag = np.zeros_like(flag)
        snr_all[x0:x1] = snr
        if thresh > nsig and verbose:
            print('  slice {:2d}: threshold raised to {:.1f} sigma to keep the masked '
                  'fraction under {:.0%}'.format(isl, thresh, maxfrac))
        flags.append(flag)
        peaks.append(float(np.nanmax(snr)) if np.isfinite(snr).any() else 0.)

    # pass 2: a bright trace scatters light into its neighbours at the same
    # position along the slice, so seed it there (it only masks where it shows)
    seeds = [[] for _ in slices]
    nbright = 0
    if adaptive and bright_reach > 0:
        for isl, (x0, x1) in enumerate(slices):
            if peaks[isl] < bright_nsig:
                continue
            nbright += 1
            nw = x1 - x0
            for (a, b) in _runs(flags[isl]):
                for d in range(-int(bright_reach), int(bright_reach) + 1):
                    j = isl + d
                    if d == 0 or not (0 <= j < len(slices)):
                        continue
                    mw = slices[j][1] - slices[j][0]
                    ja = int(round(a * mw / float(nw)))
                    jb = int(round(b * mw / float(nw)))
                    seeds[j].append((max(ja, 0), min(max(jb, ja + 1), mw)))
    if nbright and verbose:
        print('  {} bright trace(s) above {:.0f} sigma: seeded into +/-{} slices'.format(
            nbright, bright_nsig, bright_reach))

    # pass 3: build the boxes
    for isl, (x0, x1) in enumerate(slices):
        sub = data[:, x0:x1]
        if adaptive:
            cores = [(a, b, False) for a, b in _runs(flags[isl])]
            for (a, b) in seeds[isl]:
                if not flags[isl][a:b].any():      # not already detected here
                    cores.append((a, b, True))
            for bx in trace_boxes(sub, sorted(cores), x0, ny, ybin=ybin,
                                  wing_nsig=wing_nsig, seed_nsig=seed_nsig,
                                  maxfrac=maxfrac):
                boxes.append(bx)
                mask[bx[2]:bx[3], bx[0]:bx[1]] = True
        else:
            for w0, w1 in _runs(flags[isl]):
                boxes.append((x0 + w0, x0 + w1, 0, ny))
                mask[:, x0 + w0:x0 + w1] = True

    if regfn is None:
        regfn = imfname.replace('_intf.fits', '.reg')
    if mask_out is None:
        mask_out = imfname.replace('_intf', '_smsk')
    frac = mask.mean()

    if plot:
        _plot_intf(imfname, data, slices, boxes, snr_all, nsig, verbose=verbose,
                   pngfn=(imfname.replace('_intf.fits', '_bcrmsk.png')
                          if mask_out and mask_out.endswith('_bcrmsk.fits') else None))

    if len(boxes) == 0:
        # an empty region file makes kcwi_masksky_ds9.py fail (it builds no mask),
        # so write nothing and leave any existing .reg alone
        if verbose:
            print('  nothing detected: no region file written '
                  '(lower --nsig if a source is expected)')
        if write_mask:
            fits.PrimaryHDU(np.uint8(mask)).writeto(mask_out, overwrite=True)
        return None, 0, 0.

    lines = ['box({:.1f},{:.1f},{:d},{:d},0)'.format(
        0.5 * (bx0 + bx1) + 0.5, 0.5 * (by0 + by1) + 0.5, bx1 - bx0, by1 - by0)
        for (bx0, bx1, by0, by1) in boxes]
    _write_regions(regfn, lines, os.path.basename(imfname))
    if verbose:
        print('  {} boxes, {:.2%} of the frame masked -> {}'.format(
            len(boxes), frac, os.path.basename(regfn)))

    if write_mask:
        outfile = mask_out
        fits.PrimaryHDU(np.uint8(mask)).writeto(outfile, overwrite=True)
        if verbose:
            print('  wrote {}'.format(os.path.basename(outfile)))
    return regfn, len(boxes), frac


# --------------------------------------------------------------------------
# thum mode
# --------------------------------------------------------------------------


def _cover_with_ellipses(grown, snr, max_radius=6.):
    """Describe a detection footprint with one ellipse, or several if it is long.

    The ellipse written to the region file is what the DRP finally masks, so it
    has to enclose the detected pixels. A single ellipse with its semi-axes
    capped at max_radius silently truncates an extended source -- on a real
    thumbnail a feature spanning rows 31-49 came out as an ellipse reaching only
    row 43, leaving the upper half unmasked. So when the footprint is longer than
    the cap allows, it is tiled along its major axis into overlapping ellipses,
    each enclosing its own share of the pixels.

    Args:
        grown (bool ndarray): the footprint, [ny, nx].
        snr (ndarray): significance image, used to weight the centroid.
        max_radius (float): largest semi-axis of any one ellipse [pixels].

    Returns:
        list of (xc, yc, a, b, angle_deg), 0-based centres.
    """
    ys, xs = np.nonzero(grown)
    if len(xs) == 0:
        return []
    w = np.clip(snr[grown], 0, None)
    if w.sum() <= 0:
        w = np.ones(len(xs), dtype=float)
    xc, yc = np.average(xs, weights=w), np.average(ys, weights=w)
    cxx = np.average((xs - xc)**2, weights=w)
    cyy = np.average((ys - yc)**2, weights=w)
    cxy = np.average((xs - xc) * (ys - yc), weights=w)
    th = 0.5 * np.arctan2(2 * cxy, cxx - cyy)
    ct, st = np.cos(th), np.sin(th)
    u = (xs - xc) * ct + (ys - yc) * st                  # along the major axis
    v = -(xs - xc) * st + (ys - yc) * ct

    span = u.max() - u.min()

    def _fit(ntile):
        """Ellipses for `ntile` slabs along the major axis, each scaled to
        enclose its own pixels exactly (an ellipse does not cover the corners
        of the slab, so a uniform inflation is applied)."""
        out = []
        if ntile == 1:
            groups = [np.ones(len(u), dtype=bool)]
        else:
            edges = np.linspace(u.min() - 0.5, u.max() + 0.5, ntile + 1)
            groups = [(u >= edges[i] - 0.5) & (u <= edges[i + 1] + 0.5)
                      for i in range(ntile)]
        for sel in groups:
            if sel.sum() == 0:
                continue
            um = 0.5 * (u[sel].min() + u[sel].max())
            vm = 0.5 * (v[sel].min() + v[sel].max())
            aa = max(0.5 * (u[sel].max() - u[sel].min()), 0.5) + 0.2
            bb = max(0.5 * (v[sel].max() - v[sel].min()), 0.5) + 0.2
            t = ((u[sel] - um) / aa)**2 + ((v[sel] - vm) / bb)**2
            scale = max(np.sqrt(t.max()), 1.0)
            aa, bb = aa * scale + 0.2, bb * scale + 0.2
            cx = xc + um * ct - vm * st
            cy = yc + um * st + vm * ct
            out.append((cx, cy, float(aa), float(bb), np.degrees(th)))
        return out

    ntile = max(1, int(np.ceil((span + 1.4) / (2. * max_radius))))
    for _ in range(8):
        ells = _fit(ntile)
        if not ells or max(max(e[2], e[3]) for e in ells) <= max_radius * 1.15:
            break
        ntile += 1
    return ells



def _footprint_polygons(grown, simplify=0.0, fill_holes=True):
    """Trace a footprint's outline as DS9 polygons (0-based vertices).

    An ellipse has to be inflated to cover the corners of what it encloses -- on
    a real thumbnail that masked 41% more pixels than were detected -- which the
    median filter then has to interpolate over. A polygon follows the detected
    pixels instead.

    Two things about the contour tracer matter here. Matplotlib (>= 3.8) returns
    ONE compound Path per level, holding every ring as a separate subpath, so
    reading `.vertices` straight off it concatenates the rings: two disjoint
    blobs, or a blob and the hole inside it, come back as a single polygon whose
    outline jumps between them. DS9 has no hole syntax and no multi-part
    polygon, so each ring has to be split out at its MOVETO and written on its
    own. Holes are filled first, which both removes the ambiguity (every
    remaining ring is an outer boundary) and is what you want anyway -- an
    unmasked pixel in the middle of a source is not useful to ZAP.

    Args:
        grown (bool ndarray): the footprint, [ny, nx].
        simplify (float): drop vertices that lie within this distance of the
            segment between their neighbours [pixels]. 0 keeps every vertex.
        fill_holes (bool): fill interior holes before tracing.

    Returns:
        list of ndarray [(n, 2)]: polygon vertices (x, y), 0-based. One entry
        per connected blob.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as _MplPath
    if fill_holes:
        grown = ndimage.binary_fill_holes(grown)
    fig = plt.figure()
    try:
        cs = plt.contour(grown.astype(float), levels=[0.5])
        paths = (cs.get_paths() if hasattr(cs, 'get_paths') else
                 [pp for c in cs.collections for pp in c.get_paths()])
        segs = []
        for p in paths:
            v = np.asarray(p.vertices)
            codes = getattr(p, 'codes', None)
            if codes is None:
                segs.append(v)
                continue
            starts = np.flatnonzero(np.asarray(codes) == _MplPath.MOVETO)
            if len(starts) == 0:
                segs.append(v)
                continue
            for a, b in zip(starts, list(starts[1:]) + [len(v)]):
                segs.append(v[a:b])
    finally:
        plt.close(fig)
    out = []
    for v in segs:
        if len(v) < 3:
            continue
        if simplify > 0:
            keep = [0]
            for i in range(1, len(v) - 1):
                p0, p1, p2 = v[keep[-1]], v[i], v[i + 1]
                d = np.abs(np.cross(p2 - p0, p1 - p0)) / max(np.hypot(*(p2 - p0)), 1e-9)
                if d > simplify:
                    keep.append(i)
            keep.append(len(v) - 1)
            v = v[keep]
        out.append(v)
    return out


def automask_bcrmsk(imfname, mode='spots', thresh=250., nsig=5., zmed=5, grow=1,
                    edge_trim=2, minpix=4, regfn=None, mask_out=None, write_mask=True,
                    plot=False, verbose=True, **strip_kw):
    """Make the astroscrappy protection mask for the DRP from an _intf frame.

    KCWI_DRP's RemoveCosmicRays passes <base>_intf_bcrmsk.fits to astroscrappy as
    `inmask`, the pixels detect_cosmics is to ignore. On a bright source
    astroscrappy misfires: the core is narrow and steep in exactly the way its
    Laplacian filter looks for, so parts of the object are "cleaned" away.
    Measured against the DRP's own flags on eight 600 s frames of a bright
    quasar, its flagging rate rises from 0.08% of pixels at sky level to 1.4%
    above 2500 counts -- a seventeen-fold excess over a rate that should not
    depend on brightness at all, i.e. about 150 object pixels destroyed per
    frame.

    Every pixel inside this mask is also a pixel where a real cosmic ray now
    survives, to be caught later by the per-frame rejection in kcwi_stack
    (crr_split) or by a crmsk, so the mask wants to be as small as it can be
    while still covering the pixels bright enough to be at risk. Since the risk
    tracks brightness, the default is to mask the bright spots themselves
    rather than whole traces:

        1. each slice's sky spectrum (the median across its columns at every
           wavelength) is subtracted, leaving the object;
        2. a short running median along wavelength (zmed) removes anything
           narrower than the instrumental resolution, so that a cosmic ray the
           DRP missed on the first pass cannot be protected on the second;
        3. what is left above `thresh` counts (and above nsig x the local
           scatter) is masked, grown by `grow` pixels.

    On those eight frames this covers 80% of the false flags with 0.83% of the
    frame masked, against 76% and 2.51% for full-height strips (mode='strips'),
    and gives up 3x fewer real cosmic rays -- better on both counts.

    Needs an _intf frame, so run the DRP once, build the masks, set
    `bcrmsk = True` in kcwi_blue.cfg and run it again.

    Args:
        imfname (str): the _intf.fits frame.
        mode (str): 'spots' (the default, above) or 'strips' for the
            full-height traces of the sky mask, tightened; 'strips' takes the
            automask_intf keywords through **strip_kw.
        thresh (float): counts above the slice's sky spectrum at which a pixel
            is worth protecting. 250 is where the measured flagging rate first
            runs away; it is an absolute level, so scale it if your exposures
            or binning differ a lot from 600-990 s at 2x2. Higher is safer for
            cosmic rays and covers less of the object: 1000 gives 53% coverage
            for a third of the losses, 150 gives 83% for a quarter more.
        nsig (float): the pixel must also be this far above the local scatter,
            so a short exposure is not over-masked.
        zmed (int): running median along wavelength [pixels]. Must be wider
            than a cosmic ray and narrower than a real spectral feature.
        grow (int): pixels of growth around each spot.
        minpix (int): smallest spot kept.
        regfn (str): region file. Default <base>_bcrmsk.reg -- deliberately not
            <base>.reg, which is the sky mask's.
        mask_out (str): FITS mask. Default <base>_intf_bcrmsk.fits, the name
            RemoveCosmicRays looks for.
        write_mask (bool): write that FITS. On by default: unlike the sky mask,
            the FITS is what the DRP actually reads.

    Returns:
        (regfn, nboxes, frac)
    """
    base = imfname[:-len('_intf.fits')] if imfname.endswith('_intf.fits') \
        else imfname.replace('.fits', '')
    regfn = regfn or (base + '_bcrmsk.reg')
    mask_out = mask_out or (base + '_intf_bcrmsk.fits')
    if mode == 'strips':
        kw = dict(nsig=20., grow_nsig=10., grow=0, maxfrac=0.25, ybin=100,
                  wing_nsig=10., seed_nsig=15., bright_nsig=100., bright_reach=0)
        kw.update(strip_kw)
        return automask_intf(imfname, edge_trim=edge_trim, adaptive=True,
                             regfn=regfn, write_mask=write_mask, mask_out=mask_out,
                             plot=plot, verbose=verbose, **kw)
    if mode != 'spots':
        raise ValueError("mode must be 'spots' or 'strips'")

    data = fits.getdata(imfname)
    if data.ndim != 2:
        raise ValueError('{} is not a 2D _intf image'.format(imfname))
    ny, nx = data.shape
    slices = find_slices(data, edge_trim=edge_trim)
    if verbose:
        print('{}: {} x {}, {} slices'.format(os.path.basename(imfname), nx, ny, len(slices)))
    mask = np.zeros((ny, nx), dtype=bool)
    for (x0, x1) in slices:
        sub = data[:, x0:x1].astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            resid = sub - np.nanmedian(sub, axis=1)[:, None]      # drop the sky spectrum
            if zmed and zmed > 1:
                # anything narrower than this is a cosmic ray, not the object
                resid = ndimage.median_filter(resid, size=(int(zmed), 1), mode='nearest')
            sig = _robust_sigma(resid[np.isfinite(resid)])
        lim = thresh
        if np.isfinite(sig) and sig > 0:
            lim = max(thresh, nsig * sig)
        f = np.isfinite(resid) & (resid > lim)
        if grow > 0:
            f = ndimage.binary_dilation(f, np.ones((3, 3), bool), iterations=int(grow))
        mask[:, x0:x1] = f

    lab, n = ndimage.label(mask, structure=np.ones((3, 3), bool))
    polys = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None:
            continue
        pix = lab[sl] == i
        if int(pix.sum()) < minpix:
            mask[sl] &= ~pix
            continue
        # the exact pixel footprint, not its bounding box: a box around a thin
        # diagonal spot covers far more than was detected, and the region file
        # has to rebuild the same mask when --convert reads it back
        y0, x0 = sl[0].start, sl[1].start
        # pad, or a footprint that fills its bounding box has no 0.5 crossing
        # for the contour tracer to follow and comes back empty
        got = _footprint_polygons(np.pad(pix, 1))
        if got:
            for v in got:
                polys.append([(x + x0 - 1, y + y0 - 1) for (x, y) in v])
        else:
            ys, xs = np.where(pix)
            polys.append([(xs.min() + x0 - 0.5, ys.min() + y0 - 0.5),
                          (xs.max() + x0 + 0.5, ys.min() + y0 - 0.5),
                          (xs.max() + x0 + 0.5, ys.max() + y0 + 0.5),
                          (xs.min() + x0 - 0.5, ys.max() + y0 + 0.5)])
    frac = mask.mean()

    if plot:
        _plot_spots(imfname, data, slices, mask, thresh, verbose=verbose)

    if len(polys) == 0:
        if verbose:
            print('  nothing above {:g} counts: no region file written'.format(thresh))
        if write_mask:
            fits.PrimaryHDU(np.uint8(mask)).writeto(mask_out, overwrite=True)
        return None, 0, 0.

    lines = ['polygon(' + ','.join('{:.2f}'.format(c + 1) for xy in v for c in xy) + ')'
             for v in polys]
    _write_regions(regfn, lines, os.path.basename(imfname))
    if verbose:
        print('  {} spot(s) above {:g} counts, {:.2%} of the frame masked -> {}'.format(
            len(polys), thresh, frac, os.path.basename(regfn)))
    if write_mask:
        fits.PrimaryHDU(np.uint8(mask)).writeto(mask_out, overwrite=True)
        if verbose:
            print('  wrote {}'.format(os.path.basename(mask_out)))
    return regfn, len(polys), frac

def automask_thum(imfname, nsig=4., grow=1, grow_nsig=2., minpix=2, maxfrac=0.4,
                  edge=2, max_radius=6., shape='polygon', regfn=None, write_mask=False,
                  plot=False, verbose=True, data=None, cover=None, coord='image',
                  mask_out=None, edgemask=None, label=None):
    """Detect sources in a white-light thumbnail and write the median-filter mask.

    The thumbnail has one column per slice, so each column's own level is removed
    before detection (that is the striping seen in the image). Growth is kept
    small on purpose: every masked pixel has to be interpolated over by the
    median filter, so a generous mask costs accuracy.

    Args:
        imfname (str): the _icube.thum.fits image.
        nsig (float): detection threshold [MAD of the de-striped image].
        grow (int): pixels added around each detection after it has been grown
            down to grow_nsig (keep this at 0-2).
        grow_nsig (float): the mask follows connected pixels above this level.
        minpix (int): smallest detection kept, in pixels above nsig.
        maxfrac (float): the most of one slice (column) that may be masked. The
            faintest sources are dropped first when a column would exceed it.
        edge (int): drop detections within this many pixels of the edge of the
            field, where the thumbnail is unreliable.
        max_radius (float): largest semi-axis of any one ellipse [pixels]. A
            longer source is tiled into several ellipses rather than truncated.
        shape (str): 'polygon' (default: follows the detected pixels, so the
            median filter interpolates over as little as possible) or 'ellipse'
            (easier to adjust in DS9, but has to be inflated to cover what it
            encloses -- about 40% more pixels on a real thumbnail).
        regfn (str): output region file. Default: the image name with .fits
            replaced by .reg, which is what kcwi_makemask_medfilter expects.
        write_mask (bool): also write <base>_icube.mask.fits (uint8, 1 = masked).
        plot (bool): write a quick-look PNG.
        data, cover (ndarray): supply the image and its usable area directly,
            instead of reading them from imfname (used by automask_wl).
        coord (str): 'image' or 'physical' header line for the region file.
        mask_out (str): name for the mask FITS, overriding the default.
        edgemask (bool ndarray): pixels KSkyWizard already treats as edge; they
            are written as 1 and the detections as 2, as its own code does.
        label (str): name shown in the messages and on the plot.

    Returns:
        (regfn, nsources, frac)
    """
    if data is None:
        data = fits.getdata(imfname).astype(float)
    if data.ndim != 2:
        raise ValueError('{} is not a 2D image'.format(imfname))
    ny, nx = data.shape
    if cover is None:
        cover = np.isfinite(data) & (data != 0)
    label = label or os.path.basename(imfname)
    if verbose:
        print('{}: {} x {} ({} slices, {:.0%} usable)'.format(
            label, nx, ny, nx, cover.mean()))

    # de-stripe: remove each slice's own level
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        col = np.nanmedian(np.where(cover, data, np.nan), axis=0)
    col = np.where(np.isfinite(col), col, 0.)
    resid = np.where(cover, data - col[None, :], np.nan)
    sig = _robust_sigma(resid[cover])
    if not np.isfinite(sig) or sig <= 0:
        if verbose:
            print('  cannot estimate the noise; skipping')
        return None, 0, 0.

    # the usable field: coverage eroded by `edge` pixels
    inner = cover
    if edge > 0:
        inner = ndimage.binary_erosion(cover, np.ones((3, 3), bool), iterations=int(edge))

    snr = np.where(cover, resid / sig, 0.)
    # two thresholds: a core above nsig marks a real source, and the mask then
    # follows the connected pixels above grow_nsig. That keeps the mask on the
    # source's own footprint instead of a blind dilation, which matters because
    # the median filter has to interpolate over everything masked.
    core = (snr > nsig) & inner
    cand = (snr > grow_nsig) & inner
    lab, n = ndimage.label(cand, structure=np.ones((3, 3), bool))
    if n == 0 or not core.any():
        if verbose:
            print('  nothing detected: no region file written '
                  '(lower --nsig if a source is expected)')
        return None, 0, 0.

    # keep components holding at least minpix core pixels, strongest first
    objs = []
    for i in range(1, n + 1):
        pix = (lab == i)
        ncore = int((pix & core).sum())
        if ncore < minpix:
            continue
        objs.append({'pix': pix, 'npix': int(pix.sum()), 'peak': float(np.nanmax(snr[pix]))})
    if len(objs) == 0:
        if verbose:
            print('  nothing detected above the core threshold: no region file written')
        return None, 0, 0.
    objs.sort(key=lambda o: -o['peak'])

    # Describe each detection as the ellipse that encloses its pixels, and accept
    # it only while no slice goes over its budget. The budget is applied to the
    # ellipse, because the mask the DRP ends up using is built from the region
    # file, not from the detected pixels.
    yy_g, xx_g = np.mgrid[0:ny, 0:nx]
    mask = np.zeros((ny, nx), dtype=bool)
    colmax = max(1, int(np.floor(maxfrac * ny)))
    kept, skipped = [], []
    strut = np.ones((3, 3), bool)
    for o in objs:
        grown = o['pix']
        if grow > 0:
            grown = ndimage.binary_dilation(grown, strut, iterations=int(grow))
        grown &= cover
        if shape == 'polygon':
            polys = _footprint_polygons(grown)
            em = np.zeros((ny, nx), dtype=bool)
            from matplotlib.path import Path
            pts = np.column_stack([xx_g.ravel(), yy_g.ravel()])
            for v in polys:
                em |= Path(v).contains_points(pts).reshape(ny, nx)
            em = (em | grown) & cover
            # An accepted detection is masked as a solid shape. `cover` drops
            # pixels KSkyWizard calls edge (mean FLAGS > edge_var), and on a
            # bright source the core is exactly where the flags are raised --
            # on a real Q1700-QSO frame four pixels beside the peak sat at
            # 1.1-1.9 and were dropped, leaving a hole through the middle of
            # the QSO in the mask, the per-slice budget and the written
            # _zap_mask.fits while the region polygon covered it. Filling
            # interior holes closes that without pushing the outer boundary
            # into the edge region.
            em = ndimage.binary_fill_holes(em)
            o['polys'] = polys
            o['ell'] = (float(np.mean([v[:, 0].mean() for v in polys])) if polys else 0.,
                        float(np.mean([v[:, 1].mean() for v in polys])) if polys else 0.,
                        0., 0., 0.)
            o['ells'] = []
        else:
            ells = _cover_with_ellipses(grown, snr, max_radius=max_radius)
            em = np.zeros((ny, nx), dtype=bool)
            for (cx, cy, aa, bb, ang) in ells:
                th = np.radians(ang)
                ug = (xx_g - cx) * np.cos(th) + (yy_g - cy) * np.sin(th)
                vg = -(xx_g - cx) * np.sin(th) + (yy_g - cy) * np.cos(th)
                em |= ((ug / aa)**2 + (vg / bb)**2 <= 1.0)
            em &= cover
            em = ndimage.binary_fill_holes(em)   # see the note in the polygon branch
            o['ells'] = ells
            o['polys'] = []
            o['ell'] = ells[0] if ells else (0., 0., 0., 0., 0.)
        missed = int((grown & ~em).sum())
        if missed and verbose:
            print('  [warning] {} detected pixel(s) fall outside the region shapes '
                  'near (slice {:.0f}, row {:.0f})'.format(missed, o['ell'][0], o['ell'][1]))
        trial = mask | em
        if (trial.sum(axis=0) > colmax).any():
            skipped.append(o)
            continue
        mask = trial
        kept.append(o)
    if skipped and verbose:
        where = ', '.join('(slice {:.0f}, row {:.0f})'.format(o['ell'][0], o['ell'][1])
                          for o in skipped[:6])
        print('  {} detection(s) skipped to keep every slice under {:.0%}: {}{}'.format(
            len(skipped), maxfrac, where, ' ...' if len(skipped) > 6 else ''))
        print('      raise --maxfrac if these should be masked too')

    lines = ['ellipse({:.2f},{:.2f},{:.2f},{:.2f},{:.1f})'.format(
        e[0] + 1, e[1] + 1, e[2], e[3], e[4]) for o in kept for e in o['ells']]
    for o in kept:
        for v in o.get('polys', []):
            lines.append('polygon(' + ','.join('{:.2f}'.format(c + 1)
                                               for xy in v for c in xy) + ')')

    if regfn is None:
        regfn = imfname.replace('.fits', '.reg')
    _write_regions(regfn, lines, os.path.basename(imfname), coord=coord)
    frac = mask.mean()
    if verbose:
        print('  {} source(s) in {} shape(s), {:.2%} masked (worst slice {:.0%}) -> {}'.format(
            len(kept), len(lines), frac, mask.sum(axis=0).max() / float(ny),
            os.path.basename(regfn)))

    if write_mask:
        outfile = mask_out or imfname.replace('thum', 'mask')
        if edgemask is None:
            out = np.uint8(mask)
        else:                                    # KSkyWizard: edge = 1, source = 2
            out = np.zeros(mask.shape, dtype=np.uint8)
            out[mask] = 2
            out[edgemask] = 1
        fits.PrimaryHDU(out).writeto(outfile, overwrite=True)
        if verbose:
            print('  wrote {}'.format(os.path.basename(outfile)))

    if plot:
        _plot_thum(imfname, data, resid / sig, mask, kept, skipped, nsig, maxfrac,
                   verbose=verbose, label=label, pngfn=(regfn.replace('.reg', '') +
                                                        '_automask.png'))
    return regfn, len(kept), frac



# --------------------------------------------------------------------------
# turning an existing region file into a mask (no detection)
# --------------------------------------------------------------------------

def _parse_regions(regstr, shape):
    """Rasterize a DS9 region file without the `regions` package.

    Handles box, circle, ellipse, polygon and point, with rotation, in image or
    physical coordinates (identical for KCWI frames, and what
    kcwi_masksky_ds9.py assumes). Regions are combined with OR, as that script
    does; '-' exclusions are not honoured, so draw only inclusive shapes.
    """
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx] + 1.0                   # 1-based, like DS9
    mask = np.zeros(shape, dtype=bool)
    nshape = 0
    for line in regstr.splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('global'):
            continue
        m = re.match(r'^([-+]?)\s*([a-zA-Z]+)\s*\(([^)]*)\)', line)
        if m is None:
            continue
        kind = m.group(2).lower()
        try:
            p = [float(v) for v in m.group(3).replace(' ', '').split(',')]
        except ValueError:
            continue
        if kind == 'box' and len(p) >= 4:
            xc, yc, w, h = p[:4]
            ang = np.radians(p[4]) if len(p) > 4 else 0.
            u = (xx - xc) * np.cos(ang) + (yy - yc) * np.sin(ang)
            v = -(xx - xc) * np.sin(ang) + (yy - yc) * np.cos(ang)
            reg = (np.abs(u) <= w / 2.) & (np.abs(v) <= h / 2.)
        elif kind == 'circle' and len(p) >= 3:
            reg = ((xx - p[0])**2 + (yy - p[1])**2) <= p[2]**2
        elif kind == 'ellipse' and len(p) >= 4:
            xc, yc, a, b = p[:4]
            ang = np.radians(p[4]) if len(p) > 4 else 0.
            u = (xx - xc) * np.cos(ang) + (yy - yc) * np.sin(ang)
            v = -(xx - xc) * np.sin(ang) + (yy - yc) * np.cos(ang)
            reg = (u / max(a, 1e-6))**2 + (v / max(b, 1e-6))**2 <= 1.
        elif kind == 'polygon' and len(p) >= 6:
            from matplotlib.path import Path
            pts = np.array(p).reshape(-1, 2)
            inside = Path(pts).contains_points(
                np.column_stack([xx.ravel(), yy.ravel()]))
            reg = inside.reshape(shape)
        elif kind == 'point' and len(p) >= 2:
            reg = (np.abs(xx - p[0]) < 0.5) & (np.abs(yy - p[1]) < 0.5)
        else:
            continue
        mask |= reg
        nshape += 1
    return mask, nshape


def regions_to_mask(regfn, shape, use_regions_pkg=True, verbose=True, header=None):
    """Build a boolean mask from a DS9 region file.

    Uses kcwikit.kcwi.kcwi_regions (astropy `regions` or `pyregion`) when one
    of those packages is installed, so the result matches kcwi_masksky_ds9.py
    exactly, including sky coordinates (needs `header`), '-' exclusions and
    shapes outside the image; otherwise falls back to the parser above.
    """
    with open(regfn, 'r') as f:
        regstr = f.read()
    if use_regions_pkg:
        try:
            from kcwikit.kcwi.kcwi_regions import regions_to_mask as _shared_regions_to_mask
            mask = _shared_regions_to_mask(regfn, header=header, shape=shape)
            nshape = sum(1 for l in regstr.splitlines()
                         if l.strip() and not l.lstrip().startswith('#')
                         and '(' in l)
            return mask, nshape
        except ImportError:
            if verbose:
                print('  (neither regions nor pyregion is installed; using the built-in parser)')
    return _parse_regions(regstr, shape)


def convert_region(imfname, regfn=None, thum=None, kind=None, edge_var=1.,
                   verbose=True):
    """Write the mask FITS from an existing region file, leaving the .reg alone.

    Use this after hand-editing a region file in DS9: nothing is re-detected and
    the region file is only read.

    Args:
        imfname (str): the _intf.fits frame or _icube.thum.fits image whose
            dimensions the mask takes.
        regfn (str): the region file. Default: <base>.reg for an _intf frame,
            or the image name with .fits -> .reg for a thumbnail.
        thum (bool): thumbnail mode; inferred from the file name when None.
        kind (str): 'intf', 'thum', 'wl' or 'bcrmsk'; inferred from the file name
            when None ('bcrmsk' never is -- pass it explicitly, or use
            --bcrmsk --convert).
        edge_var (float): for 'wl', the second-extension value above which a
            pixel is edge, as in KSkyWizard.

    Returns:
        (outfile, nshapes, frac) or (None, 0, 0.) if the region file is missing.
    """
    base = os.path.basename(imfname)
    if kind is None:
        kind = ('wl' if imfname.endswith('_wlimg.fits') else
                'thum' if 'thum' in base else 'intf')
    if thum is not None and kind == 'intf' and thum:
        kind = 'thum'
    if regfn is None:
        regfn = (imfname.replace('_wlimg.fits', '.reg') if kind == 'wl' else
                 imfname.replace('.fits', '.reg') if kind == 'thum' else
                 imfname.replace('_intf.fits', '_bcrmsk.reg') if kind == 'bcrmsk' else
                 imfname.replace('_intf.fits', '.reg'))
    if not os.path.exists(regfn):
        if verbose:
            print('{}: {} not found, skipping'.format(
                os.path.basename(imfname), os.path.basename(regfn)))
        return None, 0, 0.
    hdul = fits.open(imfname)
    shape = hdul[0].data.shape
    mask, nshape = regions_to_mask(regfn, shape, verbose=verbose, header=hdul[0].header)
    if kind == 'wl':
        outfile = imfname.replace('_wlimg.fits', '_zap_mask.fits')
        edgemask = ~np.isfinite(hdul[0].data)
        if len(hdul) > 1 and hdul[1].data is not None and hdul[1].data.shape == shape:
            edgemask |= ~np.isfinite(hdul[1].data) | (hdul[1].data > edge_var)
        out = np.zeros(shape, dtype=np.uint8)       # KSkyWizard: edge 1, source 2
        out[mask] = 2
        out[edgemask] = 1
    elif kind == 'thum':
        outfile = imfname.replace('thum', 'mask')
        out = np.uint8(mask)
    elif kind == 'bcrmsk':
        # RemoveCosmicRays tests `bcrinmask == 1`, so anything else is unmasked
        outfile = imfname.replace('_intf.fits', '_intf_bcrmsk.fits')
        out = np.uint8(mask)
    else:
        outfile = imfname.replace('_intf', '_smsk')
        out = np.uint8(mask)
    hdul.close()
    if os.path.abspath(outfile) == os.path.abspath(imfname):
        print('  refusing to overwrite the input image: cannot name the output for '
              '{}'.format(base))
        return None, 0, 0.
    fits.PrimaryHDU(out).writeto(outfile, overwrite=True)
    if verbose:
        print('{}: {} region(s) from {} -> {} ({:.2%} masked)'.format(
            os.path.basename(imfname), nshape, os.path.basename(regfn),
            os.path.basename(outfile), mask.mean()))
    return outfile, nshape, mask.mean()



def automask_wl(imfname, edge_var=1., min_usable=0.2, trim_nsig=5., nsig=4.,
                shape='polygon', max_radius=10., grow_nsig=3., maxfrac=1.0,
                edge=0, write_mask=False, regfn=None, verbose=True, **kw):
    """Make the KSkyWizard region file from a white-light image.

    KSkyWizard writes <prefix>_<index>_wlimg.fits and expects you to draw the
    source regions by hand and save them as <prefix>_<index>.reg beside it; it
    then builds the ZAP mask with pyregion, marking its own edge pixels (where
    the second extension exceeds `edge_var`) as 1 and the region pixels as 2.
    This writes that region file automatically, in physical coordinates, as
    polygons that follow the detected pixels (pass shape='ellipse' for the
    rounded shapes the GUI itself writes).

    Detection is the same as in thumbnail mode, except that the usable area
    comes from KSkyWizard's own edge criterion, which keeps the detector edges
    out (on a real frame an artifact strip along the boundary row would
    otherwise register as a source).

    Args:
        imfname (str): the _wlimg.fits image.
        edge_var (float): second-extension value above which a pixel is edge,
            as in KSkyWizard's `edgemask = mhdu[1].data > 1`. That extension is
            the mean of the DRP FLAGS along wavelength, so a frame whose flags
            are raised at most wavelengths (a CR mask, say) can exceed it
            everywhere.
        min_usable (float): if the flag cut leaves less than this fraction of the
            finite pixels, it is ignored and detection runs on all of them.
        trim_nsig (float): how far off the interior level the median of a
            boundary row has to be for it to be dropped as illumination ramp
            rather than kept as sky. See `_trim_edge_ramp`.
        edge (int): detections this close to the edge of the usable area are
            dropped. 0 here, unlike in thumbnail mode: the level trim already
            removes the ramp, and a blind erosion on top of it was throwing
            away a real source two rows inside the field.
        maxfrac (float): most of one slice that may be masked. 1.0 (the default
            here) means no limit: ZAP estimates the background per pixel with a
            PCA, so unlike the median filter it does not need part of each slice
            left unmasked.
        nsig (float): detection threshold. On a real frame the three obvious
            sources came out at 25-59 sigma and two faint ones that had been
            masked by hand at 3.3-3.8, so use --nsig 3 to catch those.
        shape (str): 'polygon' (the default) or 'ellipse'. Polygons follow the
            detected pixels; on three real Q1700-QSO white-light frames they
            masked 8.4/7.6/7.8% of the field against 13.4/12.3/12.6% for
            ellipses, for the same detections -- an ellipse has to be inflated
            to cover the corners of what it encloses, and an elongated source
            gets tiled into several overlapping ellipses on top of that. Use
            'ellipse' to match what the KSkyWizard GUI writes by hand.
        max_radius (float): largest semi-axis of one ellipse [pixels]; 10 suits
            the white-light grid, where 6 split each source into many tiles.
            Ignored when shape='polygon'.
        grow_nsig (float): growth threshold; 3 keeps neighbouring sources from
            bridging into one long chain, which 2 did on a real frame.
        write_mask (bool): also write <base>_zap_mask.fits in KSkyWizard's
            format (edge = 1, source = 2), which the GUI would otherwise build.

    Returns:
        (regfn, nsources, frac)
    """
    hdul = fits.open(imfname)
    data = hdul[0].data.astype(float)
    finite = np.isfinite(data)
    cover = finite.copy()
    edgemask = ~finite
    if len(hdul) > 1 and hdul[1].data is not None and hdul[1].data.shape == data.shape:
        # KSkyWizard's second extension is the mean of the DRP FLAGS along
        # wavelength, and it calls a pixel edge when that exceeds edge_var
        bad = ~np.isfinite(hdul[1].data) | (hdul[1].data > edge_var)
        edgemask = edgemask | bad
        kept = finite & ~bad
        if kept.sum() >= min_usable * finite.sum():
            cover = kept
        elif verbose:
            # whole-frame flags (a CR mask applied to every wavelength, say) can
            # put the mean FLAGS above edge_var everywhere; detecting on the
            # finite pixels is better than skipping the frame
            print('{}: mean FLAGS > {:g} over {:.0%} of the frame; detecting on all '
                  'finite pixels instead (raise --edge-var to use the flag cut)'.format(
                      os.path.basename(imfname), edge_var, 1 - kept.sum() / max(finite.sum(), 1)))
    hdul.close()
    # the outermost rows along the slice are the illumination ramp, not sky:
    # drop them by level, because the flag cut above removes them only when it
    # is usable and `edge` only when it happens to be deep enough
    cover, trimmed = _trim_edge_ramp(data, cover, nsig=trim_nsig)
    if verbose and any(sum(t) for t in trimmed):
        print('{}: trimmed {} row(s) and {} column(s) of edge ramp'.format(
            os.path.basename(imfname), sum(trimmed[0]), sum(trimmed[1])))
    base = imfname[:-len('_wlimg.fits')] if imfname.endswith('_wlimg.fits') \
        else imfname.replace('.fits', '')
    return automask_thum(imfname, data=data, cover=cover, coord='physical',
                         nsig=nsig, shape=shape, max_radius=max_radius, edge=edge,
                         grow_nsig=grow_nsig, maxfrac=maxfrac, write_mask=write_mask,
                         regfn=regfn or (base + '.reg'),
                         mask_out=base + '_zap_mask.fits', edgemask=edgemask,
                         label=os.path.basename(imfname), verbose=verbose, **kw)


# --------------------------------------------------------------------------
# quick-look plots
# --------------------------------------------------------------------------

def _zscale(img, contrast=0.25):
    v = img[np.isfinite(img) & (img != 0)]
    if v.size == 0:
        return 0., 1.
    lo, hi = np.percentile(v, [2, 98])
    med = np.median(v)
    return med - (med - lo) / contrast * 0.25, med + (hi - med) / contrast * 0.25


def _plot_intf(imfname, data, slices, boxes, snr, nsig, verbose=True, pngfn=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ny, nx = data.shape
    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    vmin, vmax = _zscale(data)
    axs[0].imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax,
                  aspect='auto', interpolation='nearest')
    from matplotlib.patches import Rectangle
    for (bx0, bx1, by0, by1) in boxes:
        axs[0].add_patch(Rectangle((bx0 - 0.5, by0 - 0.5), bx1 - bx0, by1 - by0,
                                   color='red', alpha=0.25, lw=0))
        axs[1].axvspan(bx0 - 0.5, bx1 - 0.5, color='red', alpha=0.12, lw=0)
    for (x0, x1) in slices:
        axs[1].axvline(x0 - 0.5, color='0.8', lw=0.5)
    axs[0].set_ylabel('wavelength [pix]')
    axs[0].set_title('{} -- {} masked column group(s)'.format(
        os.path.basename(imfname), len(boxes)), fontsize=10)
    axs[1].plot(np.arange(nx), snr, color='C0', lw=0.8)
    axs[1].axhline(nsig, color='red', ls='--', lw=0.8, label='{:.0f} sigma'.format(nsig))
    axs[1].set_ylim(-3, max(10, min(np.nanmax(snr) * 1.1, 60)))
    axs[1].set_xlim(0, nx)
    axs[1].set_xlabel('x [pix]')
    axs[1].set_ylabel('column S/N')
    axs[1].legend(fontsize=8, loc='upper right')
    pngfn = pngfn or imfname.replace('_intf.fits', '_automask.png')
    fig.savefig(pngfn, dpi=110, bbox_inches='tight')
    plt.close(fig)
    if verbose:
        print('  plot -> {}'.format(os.path.basename(pngfn)))


def _plot_spots(imfname, data, slices, mask, thresh, verbose=True):
    """Quick look at the astroscrappy protection mask: frame plus what is covered."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ny, nx = data.shape
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                           gridspec_kw={'height_ratios': [4, 1]})
    lo, hi = _zscale(data)
    ax[0].imshow(data, origin='lower', cmap='Greys_r', vmin=lo, vmax=hi,
                 aspect='auto', interpolation='nearest')
    if mask.any():
        ax[0].contour(mask.astype(float), levels=[0.5], colors=['#c1666b'], linewidths=0.6)
    ax[0].set_ylabel('wavelength [pix]')
    ax[0].set_title('{} -- {:.2%} masked above {:g} counts'.format(
        os.path.basename(imfname), mask.mean(), thresh))
    ax[1].plot(np.arange(nx), 100. * mask.mean(axis=0), color='#1b4965', lw=0.8)
    for (x0, x1) in slices:
        ax[1].axvspan(x0, x1, color='0.93', zorder=0)
    ax[1].set_xlabel('x [pix]')
    ax[1].set_ylabel('% of column\nmasked')
    ax[1].set_xlim(0, nx)
    fig.tight_layout()
    pngfn = imfname.replace('_intf.fits', '_bcrmsk.png')
    fig.savefig(pngfn, dpi=120)
    plt.close(fig)
    if verbose:
        print('  plot -> {}'.format(os.path.basename(pngfn)))


def _plot_thum(imfname, data, snr, mask, objs, skipped, nsig, maxfrac, verbose=True,
               label=None, pngfn=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ny, nx = data.shape
    # equal aspect: these frames are tall and narrow (34 x 96 on the white-light
    # grid) and stretching them to fill a square axes made a source that is
    # elongated along the slice look round, which is not what ds9 shows
    fig, axs = plt.subplots(1, 3, figsize=(11, 5),
                            gridspec_kw={'width_ratios': [1, 1, 1.4]})
    vmin, vmax = _zscale(data)
    axs[0].imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax,
                  aspect='equal', interpolation='nearest')
    axs[0].set_title('thumbnail', fontsize=9)
    axs[1].imshow(snr, origin='lower', cmap='gray', vmin=-2, vmax=max(5, nsig * 3),
                  aspect='equal', interpolation='nearest')
    from matplotlib.patches import Ellipse, Polygon
    # Draw the shapes that go into the region file, not a contour of the
    # rasterised mask: the contour follows pixel edges and reports any hole in
    # the mask as a second closed loop, which is what made a small ring appear
    # in the middle of the QSO.
    for o in objs:
        for v in o.get('polys', []):
            axs[1].add_patch(Polygon(v, closed=True, fill=False, ec='red', lw=0.9))
        for (xc, yc, aa, bb, ang) in o.get('ells', []):
            axs[1].add_patch(Ellipse((xc, yc), 2 * aa, 2 * bb, angle=ang,
                                     fill=False, ec='red', lw=0.9))
    for o in skipped:
        for v in o.get('polys', []):
            axs[1].add_patch(Polygon(v, closed=True, fill=False, ec='cyan',
                                     ls='--', lw=0.9))
        for (xc, yc, aa, bb, ang) in o.get('ells', []):
            axs[1].add_patch(Ellipse((xc, yc), 2 * aa, 2 * bb, angle=ang,
                                     fill=False, ec='cyan', ls='--', lw=0.9))
    axs[1].set_title('de-striped S/N: {} source(s) masked{}'.format(
        len(objs), ', {} skipped (cyan)'.format(len(skipped)) if skipped else ''), fontsize=9)
    frac = mask.sum(axis=0) / float(ny)
    axs[2].bar(np.arange(nx), frac, color='C0')
    if maxfrac < 1.:
        axs[2].axhline(maxfrac, color='red', ls='--', lw=0.8,
                       label='budget {:.0%}'.format(maxfrac))
    axs[2].set_xlabel('slice')
    axs[2].set_ylabel('masked fraction')
    axs[2].set_ylim(0, max(maxfrac * 1.5 if maxfrac < 1. else 0.1,
                           frac.max() * 1.2 if frac.max() > 0 else 0.1))
    if maxfrac < 1.:
        axs[2].legend(fontsize=8)
    axs[2].set_title('per-slice budget' if maxfrac < 1. else
                     'masked fraction per slice (no budget)', fontsize=9)
    axs[2].set_box_aspect(1.6)
    fig.suptitle(label or os.path.basename(imfname), fontsize=10)
    fig.tight_layout()
    pngfn = pngfn or imfname.replace('.fits', '_automask.png')
    fig.savefig(pngfn, dpi=110, bbox_inches='tight')
    plt.close(fig)
    if verbose:
        print('  plot -> {}'.format(os.path.basename(pngfn)))


# --------------------------------------------------------------------------
# frame lists
# --------------------------------------------------------------------------

def frames_from_json(jsonfn, thum=False, wl=False, allfiles=False, verbose=True):
    """Build the list of frames from a grouping JSON, as the open_ds9 scripts do.

    The JSON only ever selects which exposures to work on here. Each frame is
    treated on its own, so every one gets its own region file; applying one
    group's region file to the whole group is what kcwi_makemask_medfilter is
    for, and duplicating it here just gave two ways to do the same thing.

    Args:
        jsonfn (str): {group: [raw file names]} as written by kcrm_group_frames.
        thum (bool): return <base>_icube.thum.fits instead of <base>_intf.fits.
        wl (bool): return <base>_wlimg.fits (KSkyWizard white-light images).
        allfiles (bool): every exposure listed in the JSON. The default is the
            first file of each group, which is the frame the region file is
            conventionally drawn on.

    Returns:
        list of str.
    """
    with open(jsonfn, 'r') as f:
        groups = json.load(f)
    suffix = ('_wlimg.fits' if wl else '_icube.thum.fits') if (thum or wl) else '_intf.fits'
    out = []
    for group, files in groups.items():
        if not files:
            continue
        for fn in (files if allfiles else files[:1]):
            stem = fn[:-5] if fn.endswith('.fits') else fn
            cand = stem + suffix
            if os.path.exists(cand):
                out.append(cand)
            elif verbose:
                print('{} not found (group {})'.format(cand, group))
    return out


# --------------------------------------------------------------------------

def parser_init():
    p = argparse.ArgumentParser(
        description='Detect sources in KCWI frames and write the DS9 region files '
                    'and masks that kcwi_makemask / kcwi_makemask_medfilter expect.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('images', nargs='*', help='_intf.fits frames or _icube.thum.fits images')
    p.add_argument('-j', '--json', default=None,
                   help='grouping JSON; uses the first file of each group')
    p.add_argument('-a', '--all', action='store_true',
                   help='with --json, process every exposure listed rather than only '
                        'the first file of each group; each frame is treated on its '
                        'own, so it gets its own region file and its own mask')
    # --per-frame / --ungroup were a separate flag for one release; --all now
    # does exactly that, so they are kept working but out of the help
    p.add_argument('--per-frame', '--ungroup', dest='all', action='store_true',
                   help=argparse.SUPPRESS)
    p.add_argument('-t', '--thum', action='store_true',
                   help='thumbnail mode (median-filter mask); implied by a *thum*.fits name')
    p.add_argument('-b', '--bcrmsk', action='store_true',
                   help="astroscrappy protection mask: write <base>_intf_bcrmsk.fits "
                        "(and <base>_bcrmsk.reg) from <base>_intf.fits, the inmask "
                        "KCWI_DRP's RemoveCosmicRays passes to detect_cosmics")
    p.add_argument('-w', '--wl', action='store_true',
                   help='white-light mode: write the KSkyWizard region file <base>.reg '
                        'from <base>_wlimg.fits; implied by a *_wlimg.fits name')
    p.add_argument('--min-usable', type=float, default=0.2,
                   help='wl: ignore the flag cut if it leaves less than this fraction '
                        'of the finite pixels')
    p.add_argument('--edge-var', type=float, default=1.,
                   help='wl: second-extension value above which a pixel counts as edge, '
                        "as in KSkyWizard's own edgemask")
    p.add_argument('--bcr-mode', choices=['spots', 'strips'], default='spots',
                   help='bcrmsk: mask the bright spots (default) or full-height traces')
    p.add_argument('--bcr-thresh', type=float, default=250.,
                   help='bcrmsk spots: counts above the slice sky at which a pixel is '
                        'protected (higher = fewer real cosmic rays given up)')
    p.add_argument('--bcr-zmed', type=int, default=5,
                   help='bcrmsk spots: running median along wavelength [pix], wider than '
                        'a cosmic ray and narrower than a real feature')
    p.add_argument('--trim-nsig', type=float, default=5.,
                   help='wl: how far off the interior level a boundary row has to be '
                        'to be dropped as illumination ramp')
    p.add_argument('-n', '--nsig', type=float, default=None,
                   help='detection threshold [sigma] (default 8 for intf, 4 for thum)')
    p.add_argument('-g', '--grow', type=int, default=None,
                   help='growth [pixels] (default 3 for intf, 1 for thum, 2 for wl, '
                        '1 for --bcrmsk spots, 0 for --bcrmsk strips)')
    p.add_argument('--edge-trim', type=int, default=2,
                   help='intf: columns dropped at each end of a slice')
    p.add_argument('--rect', action='store_true',
                   help='intf: full-height rectangles instead of a mask whose width '
                        'follows the brightness with wavelength')
    p.add_argument('--ybin', type=int, default=None,
                   help='intf: rows per wavelength bin for the adaptive width '
                        '(default 200; 100 for --bcrmsk)')
    p.add_argument('--wing-nsig', type=float, default=None,
                   help='intf: threshold for extending the mask sideways [sigma] '
                        '(default 3; 10 for --bcrmsk)')
    p.add_argument('--seed-nsig', type=float, default=None,
                   help='intf: threshold a seeded position must reach to be masked '
                        '(default 5; 15 for --bcrmsk)')
    p.add_argument('--bright-nsig', type=float, default=None,
                   help='intf: a trace above this seeds neighbouring slices, where its '
                        'scattered light leaks')
    p.add_argument('--bright-reach', type=int, default=None,
                   help='intf: how many slices away a seed is carried (0 disables; '
                        'default 3, and 0 for --bcrmsk)')
    p.add_argument('--maxfrac', type=float, default=None,
                   help='most of one slice that may be masked (default 1.0 = no limit for '
                        'intf and for wl, where ZAP works per pixel; 0.4 for thum; '
                        '0.25 for --bcrmsk strips)')
    p.add_argument('--bkg-window', type=float, default=11.,
                   help='intf: running-median window for the slice background [columns]')
    p.add_argument('--minpix', type=int, default=2,
                   help='thum: smallest detection kept [pixels above --nsig]')
    p.add_argument('--grow-nsig', type=float, default=None,
                   help='the detection follows connected pixels/columns above this '
                        'level [sigma] (default 3 for intf, 2 for thum)')
    p.add_argument('--edge', type=int, default=None,
                   help='ignore detections this close to the edge of the field '
                        '(default 2 for thum; 0 for wl, where the ramp is trimmed by level)')
    p.add_argument('--max-radius', type=float, default=None,
                   help='largest semi-axis of one ellipse; longer sources are tiled '
                        'into several [pixels] (default 6 for thum, 10 for wl)')
    p.add_argument('--shape', choices=['polygon', 'ellipse'], default=None,
                   help='region shape for thum/wl modes (default: polygon for both). '
                        'polygon follows the detected pixels, so less is interpolated '
                        'over; ellipse is easier to adjust in DS9')
    p.add_argument('-p', '--plot', action='store_true',
                   help='write a quick-look PNG per frame')
    p.add_argument('-s', '--smsk', '--mask', dest='mask', action='store_true',
                   help='also write the mask FITS directly')
    p.add_argument('-c', '--convert', action='store_true',
                   help='do not detect anything: build the mask FITS from the existing '
                        'region file. Use this after hand-editing a .reg in DS9; the '
                        'region file is only read, never rewritten')
    p.add_argument('--regfrom', default=None,
                   help='with --convert: take the regions from this file instead of the '
                        "default name; the same file is applied to every image given")
    p.add_argument('--no-clobber', action='store_true',
                   help='skip a frame whose region file already exists, so hand edits '
                        'are never overwritten')
    p.add_argument('-q', '--quiet', action='store_true')
    return p


def main():
    args = parser_init().parse_args()
    images = list(args.images)
    if args.json:
        images += frames_from_json(args.json, thum=args.thum, wl=args.wl,
                                   allfiles=args.all, verbose=not args.quiet)
    if not images:
        print('no frames given: pass file names or --json FILE')
        return

    if args.convert:
        # region file -> mask, no detection, nothing written over the .reg
        for im in images:
            if not os.path.exists(im):
                print('skipping {}: not found'.format(im))
                continue
            is_wl = args.wl or im.endswith('_wlimg.fits')
            is_bcr = args.bcrmsk and not is_wl
            is_thum = (args.thum or ('thum' in os.path.basename(im))) and not is_wl \
                and not is_bcr
            regfn = args.regfrom
            if regfn is None and is_wl:
                regfn = im.replace('_wlimg.fits', '.reg')
            if regfn is None and is_bcr:
                regfn = im.replace('_intf.fits', '_bcrmsk.reg')
            # thum and _intf fall through: convert_region derives the frame's
            # own .reg name
            convert_region(im, regfn=regfn, edge_var=args.edge_var,
                           kind='wl' if is_wl else 'bcrmsk' if is_bcr
                                else 'thum' if is_thum else 'intf',
                           verbose=not args.quiet)
        return

    for im in images:
        if not os.path.exists(im):
            print('skipping {}: not found'.format(im))
            continue
        is_wl = args.wl or im.endswith('_wlimg.fits')
        is_bcr = args.bcrmsk and not is_wl
        is_thum = (args.thum or ('thum' in os.path.basename(im))) and not is_wl \
            and not is_bcr
        regfn = (im.replace('_wlimg.fits', '.reg') if is_wl else
                 im.replace('_intf.fits', '_bcrmsk.reg') if is_bcr else
                 im.replace('.fits', '.reg') if is_thum else
                 im.replace('_intf.fits', '.reg'))
        if args.no_clobber and os.path.exists(regfn):
            if not args.quiet:
                print('{}: {} exists, left alone (--no-clobber)'.format(
                    os.path.basename(im), os.path.basename(regfn)))
            continue
        if is_wl:
            automask_wl(im, edge_var=args.edge_var, min_usable=args.min_usable,
                        trim_nsig=args.trim_nsig,
                        nsig=4. if args.nsig is None else args.nsig,
                        grow=2 if args.grow is None else args.grow,
                        minpix=args.minpix,
                        grow_nsig=2. if args.grow_nsig is None else args.grow_nsig,
                        maxfrac=1.0 if args.maxfrac is None else args.maxfrac,
                        edge=0 if args.edge is None else args.edge,
                        max_radius=10. if args.max_radius is None else args.max_radius,
                        shape='polygon' if args.shape is None else args.shape,
                        write_mask=args.mask, plot=args.plot, verbose=not args.quiet)
        elif is_thum:
            automask_thum(im,
                          nsig=4. if args.nsig is None else args.nsig,
                          grow=1 if args.grow is None else args.grow,
                          minpix=args.minpix,
                          grow_nsig=2. if args.grow_nsig is None else args.grow_nsig,
                          maxfrac=0.4 if args.maxfrac is None else args.maxfrac,
                          edge=2 if args.edge is None else args.edge,
                          max_radius=6. if args.max_radius is None else args.max_radius,
                          shape='polygon' if args.shape is None else args.shape,
                          write_mask=args.mask, plot=args.plot, verbose=not args.quiet)
        elif is_bcr:
            if args.bcr_mode == 'strips':
                sk = dict(bkg_window=args.bkg_window)
                for k, v in [('nsig', args.nsig), ('grow_nsig', args.grow_nsig),
                             ('grow', args.grow), ('maxfrac', args.maxfrac),
                             ('ybin', args.ybin), ('wing_nsig', args.wing_nsig),
                             ('seed_nsig', args.seed_nsig), ('bright_nsig', args.bright_nsig),
                             ('bright_reach', args.bright_reach)]:
                    if v is not None:
                        sk[k] = v
                automask_bcrmsk(im, mode='strips', edge_trim=args.edge_trim,
                                plot=args.plot, verbose=not args.quiet, **sk)
            else:
                automask_bcrmsk(im, mode='spots', thresh=args.bcr_thresh,
                                nsig=5. if args.nsig is None else args.nsig,
                                zmed=args.bcr_zmed,
                                grow=1 if args.grow is None else args.grow,
                                minpix=args.minpix, edge_trim=args.edge_trim,
                                plot=args.plot, verbose=not args.quiet)
        elif '_intf.fits' in im:
            automask_intf(im,
                          nsig=8. if args.nsig is None else args.nsig,
                          grow_nsig=3. if args.grow_nsig is None else args.grow_nsig,
                          grow=3 if args.grow is None else args.grow,
                          edge_trim=args.edge_trim, adaptive=not args.rect,
                          ybin=200 if args.ybin is None else args.ybin,
                          wing_nsig=3. if args.wing_nsig is None else args.wing_nsig,
                          seed_nsig=5. if args.seed_nsig is None else args.seed_nsig,
                          bright_nsig=20. if args.bright_nsig is None else args.bright_nsig,
                          bright_reach=3 if args.bright_reach is None else args.bright_reach,
                          maxfrac=1.0 if args.maxfrac is None else args.maxfrac,
                          bkg_window=args.bkg_window,
                          write_mask=args.mask, plot=args.plot, verbose=not args.quiet)
        else:
            print('skipping {}: not an _intf.fits frame or a thumbnail'.format(im))


if __name__ == '__main__':
    main()
