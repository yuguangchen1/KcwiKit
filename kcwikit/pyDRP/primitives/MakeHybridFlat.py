from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader, strip_fname, get_master_name
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

from bokeh.plotting import figure
import numpy as np
from scipy.interpolate import LSQBivariateSpline
import os
import warnings



class MakeHybridFlat(BaseImg):
    """
    Generate hybrid illumination correction from mflat, mdome, and mtwif.

    Uses 2d-spline fits to correct internal or dome flat illumination based
    on the twilight flat.

    This primitive writes out the following files and entries are made in the 
    proc file:

        * MTWIN - a \*_mtwin.fits file and an MTWIN entry
        * MTWDO - a \*_mtwdo.fits file and an MTWDO entry

    """

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Check if the master flats have already been created.
        """

        # check for flat stack to use in generating master flat
        self.logger.info("Checking precondition for MakeHybridFlat")

        tab_twin = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata,
            target_type='MTWIN', nearest=True)
        
        tab_twdo = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata,
            target_type='MTWDO', nearest=True)
        
        if len(tab_twin) > 0 and len(tab_twdo) > 0:
            self.logger.info("Already have enough hybrid flat")
            return False
        else:
            # Are there enough original master flats? 

            tab_twif = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MTWIF', nearest=True)
            
            if len(tab_twif) < 1:
                self.logger.info("No twilight flats, expecting 1")
                return False
            
            tab_flat = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MFLAT', nearest=True)

            tab_dome = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MDOME', nearest=True)
            
            if len(tab_flat) < 1 and len(tab_dome) < 1:
                self.logger.info("No internal or dome flats, expecting 1")
                return False

            return True

    def lsq_bspline_smooth_from_xy_maps(
        self,
        image,
        xmap,
        ymap,
        x_knot_spacing,
        y_knot_spacing,
        mask=None,
        weights=None,
        kx=3,
        ky=3,
        bbox=None,
        return_spline=False,
    ):
        """
        Fit a 2D least-squares B-spline surface z = f(x, y) to an image using
        physical coordinate maps, with independent knot spacing in x and y.

        Parameters
        ----------
        image : 2D ndarray
            Noisy image values.
        xmap, ymap : 2D ndarray
            Physical x/y coordinates for each pixel, same shape as image.
        x_knot_spacing, y_knot_spacing : float
            Desired interior knot spacing in physical units.
            Larger spacing => smoother result along that axis.
        mask : 2D bool ndarray, optional
            True for valid pixels to include.
        weights : 2D ndarray, optional
            Per-pixel weights. Typically 1/sigma.
        kx, ky : int, optional
            Spline degree in x and y.
        bbox : [xmin, xmax, ymin, ymax], optional
            Fit domain. If None, inferred from data.
        return_spline : bool, optional
            If True, also return the fitted spline object.

        Returns
        -------
        smoothed : 2D ndarray
            Smoothed image evaluated at the original coordinates.
        spline : LSQBivariateSpline, optional
            Returned if return_spline=True.
        """

        image = np.asarray(image, dtype=float)
        xmap = np.asarray(xmap, dtype=float)
        ymap = np.asarray(ymap, dtype=float)

        if image.shape != xmap.shape or image.shape != ymap.shape:
            raise ValueError("image, xmap, and ymap must have the same shape")

        if x_knot_spacing <= 0 or y_knot_spacing <= 0:
            raise ValueError("x_knot_spacing and y_knot_spacing must be positive")

        finite = np.isfinite(image) & np.isfinite(xmap) & np.isfinite(ymap)
        if mask is None:
            valid = finite
        else:
            valid = finite & np.asarray(mask, dtype=bool)

        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            if weights.shape != image.shape:
                raise ValueError("weights must have the same shape as image")
            valid &= np.isfinite(weights) & (weights > 0)

        x = xmap[valid].ravel()
        y = ymap[valid].ravel()
        z = image[valid].ravel()
        w = None if weights is None else weights[valid].ravel()

        if x.size < (kx + 1) * (ky + 1):
            raise ValueError("Not enough valid points for the requested spline degree")

        if bbox is None:
            xmin, xmax = np.min(x), np.max(x)
            ymin, ymax = np.min(y), np.max(y)
        else:
            xmin, xmax, ymin, ymax = bbox

        # Interior knots only
        tx = np.arange(xmin + x_knot_spacing, xmax - x_knot_spacing / 2, x_knot_spacing)
        ty = np.arange(ymin + y_knot_spacing, ymax - y_knot_spacing / 2, y_knot_spacing)

        # Need enough room for spline degree
        if len(tx) < max(0, kx - 1):
            tx = np.array([])
        if len(ty) < max(0, ky - 1):
            ty = np.array([])

        spline = LSQBivariateSpline(x, y, z, tx=tx, ty=ty, w=w, kx=kx, ky=ky)

        smoothed = np.full_like(image, np.nan, dtype=float)
        smoothed[valid] = spline.ev(xmap[valid], ymap[valid])

        if return_spline:
            return smoothed, spline
        return smoothed
    

    def _perform(self):

        if not self.config.instrument.makehybridflat:
            self.logger.info("Skipping hybrid flat by request")
        else:

            self.logger.info("Creating hybrid illumination correction")

            camera = self.action.args.ccddata.header['CAMERA'].upper()

            # get twilight frame
            tab_twif = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MTWIF', nearest=True)
            
            twiffn = get_master_name(tab_twif, 'mtwif')
            self.logger.info("Reading image: %s" % twiffn)
            ccd_twif = kcwi_fits_reader(
                os.path.join(self.config.instrument.cwd, 'redux', twiffn))[0]
            img_twif = ccd_twif.data
            
            # get root for maps
            tab = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata, target_type='MARC',
                target_group=self.action.args.groupid)
            if len(tab) <= 0:
                self.logger.error("Geometry not solved!")
                return self.action.args

            mroot = strip_fname(tab['filename'][-1])

            # Wavelength map image
            wmf = mroot + '_wavemap.fits'
            self.logger.info("Reading image: %s" % wmf)
            wavemap = kcwi_fits_reader(
                os.path.join(self.config.instrument.cwd, 'redux', wmf))[0].data

            # Slice map image
            slf = mroot + '_slicemap.fits'
            self.logger.info("Reading image: %s" % slf)
            slicemap = kcwi_fits_reader(os.path.join(
                self.config.instrument.cwd, 'redux', slf))[0].data

            # Position map image
            pof = mroot + '_posmap.fits'
            self.logger.info("Reading image: %s" % pof)
            posmap = kcwi_fits_reader(os.path.join(
                self.config.instrument.cwd, 'redux', pof))[0].data

            # make twin
            tab_flat = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MFLAT', nearest=True)

            if len(tab_flat) > 0:

                # read mflat image
                flatfn = get_master_name(tab_flat, 'mflat')
                self.logger.info("Reading image: %s" % flatfn)
                ccd_flat = kcwi_fits_reader(
                    os.path.join(self.config.instrument.cwd, 'redux', flatfn))[0]
                img_flat = ccd_flat.data

                img_div = img_twif / img_flat

                wave_range = np.nanmax(wavemap[wavemap > 0]) - np.nanmin(wavemap[wavemap > 0])
                pos_range = np.nanmax(posmap[posmap > 0]) - np.nanmin(posmap[posmap > 0])

                img_div_fit = np.full_like(img_div, np.nan, dtype=float)

                mask = (~np.isnan(img_div)) & (img_twif > 0) & (img_flat > 0) & (img_div != 1) & \
                        (posmap > 4) & (posmap < np.nanmax(posmap) - 4)
                img_div[~mask] = np.nan

                for islice in np.unique(slicemap):
                    if islice < 0:
                        continue

                    if 'RED' in camera:
                        xknots = 400
                    else:
                        xknots = 1000

                    # add iterative sigma rejection
                    niter = 3
                    nsig = 3                    
                    clip_mask = np.ones_like(img_div, dtype=bool)
                    for iter in range(niter):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=UserWarning)
                            img_tmp, spline = self.lsq_bspline_smooth_from_xy_maps(
                                img_div,
                                xmap=wavemap,
                                ymap=posmap,
                                mask=(slicemap == islice) & mask & clip_mask,
                                x_knot_spacing=wave_range / xknots,
                                y_knot_spacing=pos_range / 8,
                                kx=3,
                                ky=1,
                                return_spline=True
                            )
                    
                        clip_mask = np.abs(img_tmp - img_div) < nsig * np.nanstd(img_tmp)

                        
                    # regenerate map
                    recover_mask = (slicemap == islice)
                    smoothed = np.full_like(img_div, np.nan, dtype=float)
                    smoothed[recover_mask] = spline.ev(wavemap[recover_mask], posmap[recover_mask])

                    img_div_fit[(slicemap == islice)] = smoothed[(slicemap == islice)]

                    if self.config.instrument.plot_level >= 1:

                        figfn = 'mtwin_{0}_{1}'.format(strip_fname(flatfn), islice)

                        xplt = wavemap[mask & (slicemap == islice)]
                        yplt = img_div[mask & (slicemap == islice)]
                        yfit = img_div_fit[mask & (slicemap == islice)]

                        good_data = np.isfinite(xplt) & np.isfinite(yplt) & np.isfinite(yfit)
                        yrange = np.nanpercentile(yplt[good_data], [1, 99])
                        p = figure(
                            title=figfn + ' Hybrid Illumination Correction',
                            x_axis_label='Wave (A)',
                            y_axis_label='Flux (e-)',
                            y_range=yrange,
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        
                        p.circle(xplt[good_data], yplt[good_data], size=1, line_alpha=0., fill_color='black',
                            legend_label='Data')
                        p.circle(xplt[good_data], yfit[good_data], size=1, line_alpha=0., fill_color='red',
                            legend_label='Fit')
                        p.legend.location = "top_left"
                        bokeh_plot(p, self.context.bokeh_session)

                        save_plot(p, filename=figfn+".png")
                        

                # create hybrid image
                img_div_fit[~np.isfinite(img_div_fit)] = 1
                img_twin = img_flat * img_div_fit
                img_twin[img_flat == 1] = 1

                ccd_twin = ccd_flat.copy()
                ccd_twin.data = img_twin
                ccd_twin.uncertainty.array = ccd_flat.uncertainty.array * img_div_fit

                newfn = strip_fname(tab_flat['filename'][-1]) + '_mtwin.fits'

                log_string = MakeHybridFlat.__module__
                ccd_twin.header['HISTORY'] = log_string
                ccd_twin.header['IMTYPE'] = 'TWIN'
                ccd_twin.header['MASTFLAT'] = (True, 'master flat image?')
                ccd_twin.header['WAVMAPF'] = wmf
                ccd_twin.header['SLIMAPF'] = slf
                ccd_twin.header['POSMAPF'] = pof


                self.logger.info("Writing image: %s" % newfn)
                kcwi_fits_writer(ccd_twin, output_file=newfn,
                                output_dir=self.config.instrument.output_directory)
                self.context.proctab.update_proctab(frame=ccd_twin, suffix='mtwin',
                                                    newtype='MTWIN',
                                                    filename=tab_flat['filename'][-1])
                self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

            # twdo flats
            tab_dome = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type='MDOME', nearest=True)
            
            if len(tab_dome) > 0:

                # read mflat image
                domefn = get_master_name(tab_dome, 'mdome')
                self.logger.info("Reading image: %s" % domefn)
                ccd_dome = kcwi_fits_reader(
                    os.path.join(self.config.instrument.cwd, 'redux', domefn))[0]
                img_dome = ccd_dome.data

                img_div = img_twif / img_dome

                wave_range = np.nanmax(wavemap[wavemap > 0]) - np.nanmin(wavemap[wavemap > 0])
                pos_range = np.nanmax(posmap[posmap > 0]) - np.nanmin(posmap[posmap > 0])

                img_div_fit = np.full_like(img_div, np.nan, dtype=float)

                mask = (~np.isnan(img_div)) & (img_twif > 0) & (img_flat > 0) & (img_div != 1) & \
                        (posmap > 4) & (posmap < np.nanmax(posmap) - 4)
                img_div[~mask] = np.nan

                for islice in np.unique(slicemap):
                    if islice < 0:
                        continue

                    if 'RED' in camera:
                        xknots = 400
                    else:
                        xknots = 1000

                    # add iterative sigma rejection
                    niter = 3
                    nsig = 3
                    clip_mask = np.ones_like(img_div, dtype=bool)
                    for iter in range(niter):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=UserWarning)
                            img_tmp, spline = self.lsq_bspline_smooth_from_xy_maps(
                                img_div,
                                xmap=wavemap,
                                ymap=posmap,
                                mask=(slicemap == islice) & mask & clip_mask,
                                x_knot_spacing=wave_range / xknots,
                                y_knot_spacing=pos_range / 8,
                                kx=3,
                                ky=1,
                                return_spline=True
                            )

                            clip_mask = np.abs(img_tmp - img_div) < nsig * np.nanstd(img_tmp)

                        
                    # regenerate map
                    recover_mask = (slicemap == islice)
                    smoothed = np.full_like(img_div, np.nan, dtype=float)
                    smoothed[recover_mask] = spline.ev(wavemap[recover_mask], posmap[recover_mask])

                    img_div_fit[(slicemap == islice)] = smoothed[(slicemap == islice)]

                    if self.config.instrument.plot_level >= 1:

                        figfn = 'mtwdo_{0}_{1}'.format(strip_fname(domefn), islice)

                        xplt = wavemap[mask & (slicemap == islice)]
                        yplt = img_div[mask & (slicemap == islice)]
                        yfit = img_div_fit[mask & (slicemap == islice)]

                        good_data = np.isfinite(xplt) & np.isfinite(yplt) & np.isfinite(yfit)
                        yrange = np.nanpercentile(yplt[good_data], [1, 99])
                        p = figure(
                            title=figfn + ' Hybrid Illumination Correction',
                            x_axis_label='Wave (A)',
                            y_axis_label='Flux (e-)',
                            y_range=yrange,
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        
                        p.circle(xplt[good_data], yplt[good_data], size=1, line_alpha=0., fill_color='black',
                            legend_label='Data')
                        p.circle(xplt[good_data], yfit[good_data], size=1, line_alpha=0., fill_color='red',
                            legend_label='Fit')
                        p.legend.location = "top_left"
                        bokeh_plot(p, self.context.bokeh_session)

                        save_plot(p, filename=figfn+".png")
                        

                # create hybrid image
                img_div_fit[~np.isfinite(img_div_fit)] = 1
                img_twdo = img_dome * img_div_fit
                img_twdo[img_dome == 1] = 1

                ccd_twdo = ccd_dome.copy()
                ccd_twdo.data = img_twdo
                ccd_twdo.uncertainty.array = ccd_dome.uncertainty.array * img_div_fit

                newfn = strip_fname(tab_dome['filename'][-1]) + '_mtwdo.fits'

                log_string = MakeHybridFlat.__module__
                ccd_twdo.header['HISTORY'] = log_string
                ccd_twdo.header['IMTYPE'] = 'TWDO'
                ccd_twdo.header['MASTFLAT'] = (True, 'master flat image?')
                ccd_twdo.header['WAVMAPF'] = wmf
                ccd_twdo.header['SLIMAPF'] = slf
                ccd_twdo.header['POSMAPF'] = pof


                self.logger.info("Writing image: %s" % newfn)
                kcwi_fits_writer(ccd_twdo, output_file=newfn,
                                output_dir=self.config.instrument.output_directory)
                self.context.proctab.update_proctab(frame=ccd_twdo, suffix='mtwdo',
                                                    newtype='MTWDO',
                                                    filename=tab_dome['filename'][-1],)
                self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

            return self.action.args

