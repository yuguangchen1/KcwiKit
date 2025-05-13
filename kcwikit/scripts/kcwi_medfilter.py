import argparse
import os
from itertools import repeat
import multiprocessing

_nproc = max(int(multiprocessing.cpu_count() - 2), 1)


#Third-party Imports
from astropy.io import fits
from astropy import table
import numpy as np
from scipy.interpolate import CubicSpline, interp1d


from tqdm import tqdm
import pdb



def parser_init():
    description = 'Conduct median filtering and subtract the median-filtered background.'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('file', type=str, nargs='+',
            help='File name')
    parser.add_argument('-p', '--par', dest='parfn',
            type=str, default='medfilter.ppar',
            help='Parameter file')
    parser.add_argument('--zbin', dest='zbin',
            type=int, help='Window size in the z-direction.')
    parser.add_argument('--ybin', dest='ybin',
            type=int, help='Window size in the y-direction.')
    parser.add_argument('--ytrim', dest='ytrim',
            type=int, help='Number of pixels to be trimmed on the edge.')
    parser.add_argument('--override', dest='override',
            action='store_true',
            help='Act on processed files.')
    parser.add_argument('--nproc', dest='nproc',
            type=int, default=_nproc,
            help='Number of processors to use.')


    args = parser.parse_args()
    return args


_default_par = {'fn':'', 'zbin':100, 'ybin':16, 'ytrim':4, 'bkgtype':0}

_pre = 'kcwi_medfilter.py'

def kcwi_medfilter_readpar(parfn):
    # Read ppar file

    with open(parfn, 'r') as f:
        lines = f.readlines()

    ppar = []
    for line in lines:
        if line.strip() == '':
            continue

        elements = line.split('=')
        key = elements[0].strip()
        value = elements[1].strip()

        if key=='FN':
            ppar.append(_default_par.copy())
            ppar[-1]['fn'] = value
        elif key=='ZBIN':
            ppar[-1]['zbin'] = int(value)
        elif key=='YBIN':
            ppar[-1]['ybin'] = int(value)
        elif key=='YTRIM':
            ppar[-1]['ytrim'] = int(value)
        elif key=='BKGTYPE':
            ppar[-1]['bkgtype'] = int(value)
        elif key=='FNNUM':
            print(_pre + ': [Warning] FNNUM deprecated in ppar file. Please use FN.')
        else:
            print(_pre + ': [Warning] ppar key not recognized - ' + key)

    return ppar


def kcwi_medfilter_master_par(args, ppar=None):
    # Generate the master parameter table

    args_par = vars(args)

    # Set same amount of par as the number of FITS files
    if isinstance(args_par['file'], str):
        args_par['file'] = [args_par['file']]
    master_par = [_default_par.copy() for i in range(len(args_par['file']))]

    # Change fnnum
    for i, par in enumerate(master_par):
        par['fn'] = args_par['file'][i]

    # convert to table
    master_par = table.Table(master_par)

    # Change master_par based on individual ppar item
    if ppar is not None:
        for i, par in enumerate(ppar):

            index = np.where((master_par['fn'] == par['fn']))[0]
            if len(index) == 0:
                print(_pre + ': [Warning] No matching file for ppar entry - ' + par['fn'])
            else:
                master_par['zbin'][index] = par['zbin']
                master_par['ybin'][index] = par['ybin']
                master_par['ytrim'][index] = par['ytrim']

    # override ppar with args
    if args_par['zbin'] is not None:
        master_par['zbin'] = args_par['zbin']
    if args_par['ybin'] is not None:
        master_par['ybin'] = args_par['ybin']
    if args_par['ytrim'] is not None:
        master_par['ytrim'] = args_par['ytrim']

    return master_par


def kcwi_medfilter_actonone(args, par):
    if os.path.isfile(par['fn']):

        # file names
        cubefn = par['fn']

        # allow second pass for cubed or cubes
        if 'cubed' in cubefn or 'cubes' in cubefn:
            flagpass = 2
        else:
            flagpass = 1

        # all versions of masks
        # 2D masks for icube.fits
        mimgfn1 = cubefn.replace('.fits', '.mask.fits')
        mimgfn2 = cubefn.replace('.fits', '_2d.mask.fits')
        # 3D masks for icubed or icubes
        mstackfn = cubefn.replace('.fits', '.stackmask.fits')

        # read cube
        hdl = fits.open(cubefn)
        cube = hdl[0].data.copy()
        hdr = hdl[0].header

        # Skip non-object
        if hdr['IMTYPE'] != 'OBJECT':
            print(_pre + ': [Warning] Not a science cube. Skipping - ' + \
                cubefn)

            return

        print(_pre + ': Processing - ' + cubefn)

        # Check header
        if 'MEDFILT' in hdr:
            if hdr['MEDFILT'] >= flagpass:
                if not vars(args)['override']:
                    print(_pre + ': [Warning] Processed before. Set -override to force processing - ' + \
                            cubefn)
                return

        # Gather cubes
        cube0 = cube.copy()
        shape = cube.shape
        if len(hdl) > 1:
            # python DRP
            mcube = hdl['MASK'].data.copy()
            fcube = hdl['FLAGS'].data.copy()
        else:
            # IDL DRP
            mcube = fits.open(cubefn.replace('icube', 'mcube'))[0].data.copy()
            fcube = mcube

        # interpolation flag
        flagcube = np.zeros_like(cube0, dtype=int)


        # read masks
        if os.path.isfile(mimgfn1):
            mimg1 = fits.open(mimgfn1)[0].data
            print("Continuum mask found: ", mimgfn1)
        else:
            mimg1 = np.zeros(shape)

        if os.path.isfile(mimgfn2):
            mimg2 = fits.open(mimgfn2)[0].data
            print("Line mask found: ", mimgfn1)
        else:
            mimg2 = None

        if os.path.isfile(mstackfn):
            mstack = fits.open(mstackfn)[0].data
        else:
            mstack = None

        # wave
        wave = (np.arange(shape[0]) - hdr['CRPIX3'] + 1) * hdr['CD3_3'] + hdr['CRVAL3']
        cwave = hdr['BCWAVE']

        # skipping flag
        if par['bkgtype'] == 1:
            print(_pre + ': BKGTYPE is 1. Skipping - ' + cubefn)
            return

        # padded
        flagcube[cube == 0] = 1
        cube[cube==0] = np.nan

        # masking
        cube[mcube != 0] = np.nan; flagcube[mcube != 0] = -100
        cube[fcube != 0] = np.nan; flagcube[fcube != 0] = -100

        # temp fix for out-of-FoV
        flagcube[fcube > 100] = 1

        # emission line mask
        if mimg2 is not None:
            # resemble cube - copied from kcwi_flatten_cube
            mcube2 = np.zeros(shape)
            trm = 0
            for ii in range(shape[2]):
                ix0 = 0 + trm
                ix1 = shape[1] - trm
                ox0 = ii * shape[1] + trm
                ox1 = (ii + 1) * shape[1] - trm
                mcube2[:, ix0:ix1, ii] = mimg2[:, ox0:ox1]
            cube[(mcube2 != 0) & (flagcube != 1)] = np.nan
            flagcube[(mcube2 != 0) & (flagcube != 1)] = -100

        # continuum mask
        qq = (mimg1 != 0)
        if np.sum(qq) != 0:
            for kk in range(shape[0]):
                cube[kk][qq & (flagcube[kk, :, :]!=1)] = np.nan
                flagcube[kk][qq & (flagcube[kk, :, :]!=1)] = -100

        # 2nd pass mask
        if mstack is not None:

            # masked pixels
            index = (mstack == 1)
            cube[index & (flagcube != 1)] = np.nan
            flagcube[index & (flagcube != 1)] = -100

            # out of FoV pixels
            index = ~np.isfinite(mstack)
            cube[index & (flagcube != 1)] = np.nan
            flagcube[index & (flagcube != 1)] = -200

            # temp fix the edges in flag cubes after DAR
            flagcube[:, 0, :] = 1
            flagcube[:, -1, :] = 1
            flagcube[:, :, 0] = 1
            flagcube[:, :, -1] = 1


        # trim
        for kk in range(shape[0]):
            for ii in range(shape[2]):
                line = cube[kk, :, ii]
                flagline = flagcube[kk, :, ii]

                # edges of unpadded value
                qgood = np.where((flagline != -200) & (flagline != 1) )[0]
                if len(qgood) > 0:
                    qmin, qmax = np.min(qgood), np.max(qgood)

                    # lower
                    smallline = line[:qmin + par['ytrim']]
                    smallflag = flagline[:qmin + par['ytrim']]
                    smallline[smallflag != 1] = np.nan
                    smallflag[smallflag != 1] = -200
                    # upper
                    smallline = line[qmax - par['ytrim'] + 1:]
                    smallflag = flagline[qmax - par['ytrim'] + 1:]
                    smallline[smallflag != 1] = np.nan
                    smallflag[smallflag != 1] = -200


        """
        # lower
        smallcube = cube[:, :par['ytrim'], :]
        smallflagcube = flagcube[:, :par['ytrim'], :]
        smallcube[smallflagcube != 1] = np.nan
        smallflagcube[smallflagcube != 1] = -200
        # upper
        smallcube = cube[:, -par['ytrim']:, :]
        smallflagcube = flagcube[:, -par['ytrim']:, :]
        smallcube[smallflagcube != 1] = np.nan
        smallflagcube[smallflagcube != 1] = -200
        """


        # construct the raw median cube
        medcube = np.zeros(shape)
        for ii in tqdm(range(shape[2]), desc=par['fn']):
            for jj in range(shape[1]):
                for kk in range(shape[0]):

                    # skip bad pixels
                    if flagcube[kk, jj, ii]!=0:
                        continue

                    # median window
                    yrange = ( int(max(par['ytrim'], (jj - par['ybin'] / 2))), \
                            int(min( (shape[1] - par['ytrim']), jj + par['ybin'] / 2)) )
                    zrange = ( int(max(0, (kk - par['zbin'] / 2))), \
                            int(min( (shape[0] - 1), (kk + par['zbin'] / 2))) )
                    medcube[kk, jj, ii] = np.nanmedian(cube[zrange[0]:zrange[1], yrange[0]:yrange[1], ii])

        # interpolate bad data
        """
        hdu_tmp = fits.PrimaryHDU(medcube)
        hdu_tmp.writeto('tmp.fits', overwrite=True)
        hdu_tmp = fits.PrimaryHDU(flagcube)
        hdu_tmp.writeto('tmp_flag.fits', overwrite=True)

        medcube = fits.open('tmp.fits')[0].data
        flagcube = fits.open('tmp_flag.fits')[0].data
        """

        # interpolate along the slices
        for ii in range(shape[2]):
            for kk in range(shape[0]):

                line0 = medcube[kk, :, ii]
                line = line0  # soft copy just to have same variable names as IDL
                flagline0 = flagcube[kk, :, ii]
                flagline = flagline0

                qv = (flagline == 0)
                if np.sum(qv) == 0:
                    # No good data
                    flagline0[(~qv) & (flagline != 1)] = -300 # interp in wavelength

                else:
                    # find good data limits
                    qm = (np.min(np.where(qv)[0]), np.max(np.where(qv)[0]))

                    # cubic interpolation
                    q100 = np.where(flagline == -100)[0]
                    if len(q100) > 0:

                        # how many need extrapolation?
                        q = np.where( (q100 < qm[0]) | (q100 > qm[1]) )
                        # convert to nearest neighbor
                        flagline[q100[q]] = -200

                        q100 = np.where(flagline == -100)[0]
                        # need >= 4 good values for cubic spline
                        if np.sum(qv) > 3:
                            cs = CubicSpline(np.where(qv)[0], line[qv])
                            line0[q100] = cs(q100)
                        else:
                            flagline0[q100] = -200


                    # nearest neighbor
                    q200 = np.where(flagline == -200)[0]
                    if len(q200) > 0:
                        if np.sum(qv) > 1:
                            nn = interp1d(np.where(qv)[0], line[qv], kind='nearest',
                                    bounds_error=False, fill_value='extrapolate')
                            line0[q200] = nn(q200)
                        elif np.sum(qv)==1:
                            line0[q200] = line[qv][0]


        # wavelength interpolation
        index3 = (flagcube == -300)
        if np.sum(index3) > 0:
            for ii in range(shape[2]):
                for jj in range(shape[1]):
                    tmpmed0 = medcube[:, jj, ii]
                    tmpmed = tmpmed0
                    tmpflag0 = flagcube[:, jj, ii]
                    tmpflag = tmpflag0
                    qlin3 = np.where(tmpflag == -300)[0]

                    if len(qlin3)==0:
                        continue

                    # extrapolate with nearest neighbor
                    qlin0 = np.where(tmpflag == 0)[0]
                    q12 = np.where( (tmpflag==-100) | (tmpflag==-200) )[0]
                    if len(qlin0) == 0:
                        # numbers from other wavelengths are interpolated as well.
                        # use nearest neighbor
                        if len(q12) > 1:
                            nn = interp1d(q12, tmpmed[q12], kind='nearest',
                                    bounds_error=False, fill_value='extrapolate')
                            tmpmed[qlin3] = nn(qlin3)
                        elif len(q12)==1:
                            tmpmed[qlin3] = tmpmed[q12][0]
                        else:
                            tmpmed[qlin3] = 0
                        qlin0 = q12


                    if len(q12) > 0:
                        qq = np.where(qlin3 < np.min(qlin0))[0]
                        tmpmed[qlin3[qq]] = tmpmed[np.min(qlin0)]
                        tmpflag[qlin3[qq]] = 0
                        qq = np.where(qlin3 > np.max(qlin0))[0]
                        tmpmed[qlin3[qq]] = tmpmed[np.max(qlin0)]
                        tmpflag[qlin3[qq]] = 0

                        qlin3 = np.where(tmpflag == -300)[0]
                        li = interp1d(qlin0, tmpmed[qlin0], kind='linear',
                                bounds_error=False, fill_value='extrapolate')
                        tmpmed[qlin3] = li(qlin3)


        # recover padding
        medcube[flagcube==1] = 0

        # subtract
        cube1 = cube0 - medcube

        # output
        hdr['MEDFILT'] = flagpass
        hdr['HISTORY'] = 'Median filtered, level = {}'.format(flagpass)
        hdr['HISTORY'] = 'YTRIM = {0:d}, YBIN = {1:d}, ZBIN = {2:d}'.format(par['ytrim'],
                par['ybin'], par['zbin'])

        # Overwrite the original cube
        hdl[0].data = cube1
        hdl.writeto(par['fn'], overwrite=True)

        # Median cube
        hdu_med = fits.PrimaryHDU(medcube, header=hdr)
        hdu_med.writeto(par['fn'].replace('.fits', '.med.fits'), overwrite=True)


    else:
        print(_pre + ': File not found - ' + par['fn'])

    return




def median_filter(args, master_par):

    # set up progress bar
    tqdm.set_lock(multiprocessing.RLock())

    # multiprocessing
    print(_pre + ': Using {0:d} processors'.format(vars(args)['nproc']))
    pool = multiprocessing.Pool(processes=vars(args)['nproc'], initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),))
    pool.starmap(kcwi_medfilter_actonone, zip(repeat(args), master_par))

    # Loop through all files
    #for i, par in enumerate(master_par):
    #    kcwi_medfilter_actonone(args, par)


    return


def main():
    #arg_parser = parser_init()
    #args = arg_parser.parse_args()

    args = parser_init()

    # Read ppar file
    if os.path.isfile(vars(args)['parfn']):
        ppar = kcwi_medfilter_readpar(vars(args)['parfn'])
    else:
        ppar = None

    # Generate a master parameter table
    master_par = kcwi_medfilter_master_par(args, ppar)
    median_filter(args, master_par)

if __name__ == '__main__':

    main()
