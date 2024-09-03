import argparse
import os

#Third-party Imports
from astropy.io import fits
import astropy.table
from astropy import io
import numpy as np
from scipy.interpolate import interp1d


import pdb


def parser_init():
    description = 'Combining inverse sensitivity curves from multiple exposures of standard stars.'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('listfile', type=str, nargs=1,
            help='List file containing all exposures.')
    parser.add_argument('--noplot', dest='noplot',
            action='store_true',
            help='Creat plot?')

    args = parser.parse_args()
    return args


_pre = 'kcwi_combinestd.py'

def combinestd(listfile, noplot=False):

    # import plotting libraries
    if not noplot:
        import colorcet as cc
        from bokeh.plotting import figure, show, output_file
        from bokeh.layouts import column
        from bokeh.models import Range1d
        from bokeh.models import Span
        # from bokeh.io import export_svg
        from bokeh.io import export_png


    # convert to string
    if isinstance(listfile, list):
        listfile = listfile[0]


    # read the list
    listtab = io.ascii.read(listfile, format='no_header', comment='\s*#')


    # setup bokeh plot
    if not noplot:
        p1 = figure(title="Effective aperture", x_axis_label='Wavelength', y_axis_label='EA')
        p2 = figure(title='Inverse sensitivity', x_axis_label='Wavelength', y_axis_label='IS')

        palette = [cc.rainbow[int(i / len(listtab) * 255)] for i in range(len(listtab))]



    # loop through all FITS files
    for i, invsensfn in enumerate(listtab['col1']):

        # read file
        hdu = fits.open(invsensfn)[0]
        data = hdu.data
        hdr= hdu.header

        # save the first file
        if i==0:
            data0 = data.copy()
            hdr0 = hdr.copy()

        # check instrument setup
        ifunam = hdr['IFUNAM']
        gratnam = hdr['BGRATNAM']
        cwave = np.round(hdr['BCWAVE'])
        pwave = np.round(hdr['BPWAVE'])
        if i==0:
            ifunam0 = ifunam
            gratnam0 = gratnam
            cwave0 = cwave
            pwave0 = pwave

        if ifunam!=ifunam0 or gratnam!=gratnam0 or cwave!=cwave0 or pwave!=pwave0:
            print(_pre + ': [Warning] Inconsistent instrument config:')
            print('   0 - IFU = {}, GRAT = {}, CWAVE = {}, PWAVE = {}'.format(ifunam0, gratnam0, cwave0, pwave0))
            print('  ' + invsensfn + ' - IFU = {}, GRAT = {}, CWAVE = {}, PWAVE = {}'.format(\
                    ifunam, gratnam, cwave, pwave))

        # Check wavelength
        wave = (np.arange(hdr['NAXIS1']) - hdr['CRPIX1'] + 1) * hdr['CDELT1'] + hdr['CRVAL1']
        dw = wave[1] - wave[0]
        if i==0:
            wave0 = wave
            dw0 = dw

        # resample
        if wave0[0]!=wave[0] or dw0!=dw or len(wave0)!=len(wave):
            li = interp1d(wave, data[0, :], kind='linear', bounds_error=False, fill_value='extrapolate')
            invsens_0 = li(wave0)
            li = interp1d(wave, data[1, :], kind='linear', bounds_error=False, fill_value='extrapolate')
            invsens_1 = li(wave0)
        else:
            invsens_0 = data[0, :]
            invsens_1 = data[1, :]

        # Plot
        if not noplot:
            p2.line(wave0, invsens_0, color=palette[i])
            p2.line(wave0, invsens_1, color=palette[i], line_width=2,legend_label=os.path.basename(invsensfn))


        # Store all data
        if i==0:
            invsens_0_all = np.zeros((len(listtab), len(wave0)))
            invsens_1_all = np.zeros((len(listtab), len(wave0)))

        invsens_0_all[i, :] = invsens_0
        invsens_1_all[i, :] = invsens_1


        # Plot EA
        if not noplot:
            hdu_ea = fits.open(invsensfn.replace('_invsens', '_ea'))[0]
            data_ea = hdu_ea.data

            qgood = (wave > hdr0['WAVGOOD0']) & (wave < hdr0['WAVGOOD1'])

            p1.line(wave, data_ea[0, :], color=palette[i])
            p1.line(wave, data_ea[1, :], color=palette[i], line_width=2,\
                    legend_label =f"{os.path.basename(invsensfn)} {np.max(data[1, qgood]):.3e}")


    # Averaging
    if len(listtab)!=1:
        data0[0, :] = np.nanmean(invsens_0_all, axis=0)
        data0[1, :] = np.nanmean(invsens_1_all, axis=0)
    else:
        data0[0, :] = invsens_0_all.flatten()
        data0[1, :] = invsens_1_all.flatten()


    # write
    hdu_new = fits.PrimaryHDU(data0, header=hdr)
    hdu_new.writeto(listfile.replace('.list', '_invsens.fits'), overwrite=True)



    if not noplot:
        # plot combined
        p2.line(wave0, data0[0, :], color='black')
        p2.line(wave0, data0[1, :], color='black', line_width=2, legend_label='Combined')


        # wavegood boundary
        wavgood0_vline = Span(location=hdr0['WAVGOOD0'],
                              dimension='height', line_color='red',
                              line_dash='dashed', line_width=2)
        p1.add_layout(wavgood0_vline)
        p2.add_layout(wavgood0_vline)
        wavgood1_vline = Span(location=hdr0['WAVGOOD1'],
                              dimension='height', line_color='red',
                              line_dash='dashed', line_width=2)
        p1.add_layout(wavgood1_vline)
        p2.add_layout(wavgood1_vline)

        # Setting range
        p1.y_range = Range1d(0, 2.5e5)
        p2.y_range = Range1d(0, 1e-16)
        col = column(p1, p2)
        output_file(filename = listfile.replace('list','html'))
        show(col)
        export_png(col, filename=listfile.replace('list','png'))


    return

def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    combinestd(**vars(args))

if __name__=='__main__':

    main()
