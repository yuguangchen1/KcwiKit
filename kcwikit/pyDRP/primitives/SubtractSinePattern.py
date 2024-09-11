from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, kcwi_fits_reader
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

from bokeh.plotting import figure
import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import time
import os


class SubtractSinePattern(BasePrimitive):
    """Subtract periodic pattern from raw image (adapted from Pypeit by YC)"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        return True

    def sec2slice(self, subarray, one_indexed=False, include_end=False, require_dim=None, binning=None):
        """
        Convert a string representation of an array subsection (slice) into
        a list of slice objects.

        Args:
            subarray (str):
                The string to convert.  Should have the form of normal slice
                operation, 'start:stop:step'.  The parser ignores whether or
                not the string has the brackets '[]', but the string must
                contain the appropriate ':' and ',' characters.
            one_indexed (:obj:`bool`, optional):
                The string should be interpreted as 1-indexed.  Default
                is to assume python indexing.
            include_end (:obj:`bool`, optional):
                **If** the end is defined, adjust the slice such that
                the last element is included.  Default is to exclude the
                last element as with normal python slicing.
            require_dim (:obj:`int`, optional):
                Test if the string indicates the slice along the proper
                number of dimensions.
            binning (:obj:`str`, optional):
                Assume the slice is for an unbinned array and adjust the
                returned slice for this binning in each dimension.  If two
                dimensional, the format of this string must be, e.g., `1,2`
                for unbinned rows and a factor of 2 binning along columns.

        Returns:
            tuple: A tuple of slice objects, one per dimension of the
            prospective array.

        Raises:
            TypeError:
                Raised if the input `subarray` is not a string.
            ValueError:
                Raised if the string does not match the required
                dimensionality or if the string does not look like a
                slice.
        """

        # Check it's a string
        if not isinstance(subarray, str):
            raise TypeError('Can only parse string-based subarray sections.')
        # Remove brackets if they're included
        sections = subarray.strip('[]').split(',')
        # Check the dimensionality
        ndim = len(sections)
        _binning = [1]*ndim if binning is None else np.array(binning.split(',')).astype(int)
        if len(_binning) != ndim:
            raise ValueError('Incorrect binning dimensions (found {0}, expected {1}).'.format(
                                len(_binning), ndim))
        if require_dim is not None and ndim != require_dim:
            raise ValueError('Number of slices ({0}) in {1} does not match '.format(ndim, subarray) + 
                            'required dimensions ({0}).'.format(require_dim))
        # Convert the slice of each dimension from a string to a slice
        # object
        slices = []
        for s,b in zip(sections,_binning):
            flipped = False
            # Must be able to find the colon
            if ':' not in s:
                raise ValueError('Unrecognized slice string: {0}'.format(s))
            # Initial conversion
            _s = [ None if x == '' else int(x) for x in s.split(':') ]
            if len(_s) > 3:
                raise ValueError('String as too many sections.  Must have format \'start:stop:step\'.')
            if len(_s) < 3:
                # Include step
                _s += [ None ]
            # Must check order first so "include_last" and "one_indexed" are correctly applied
            # Check that the first two elements of the slice are ordered correctly
            if _s[0] is not None and _s[1] is not None:
                if _s[0] > _s[1]:
                    flipped = True
                    _s = [_s[1], _s[0], _s[2]]
            if one_indexed:
                # Decrement to convert from 1- to 0-indexing
                _s = [ None if x is None else x-1 for x in _s ]
            if include_end and _s[1] is not None:
                # Increment to include last 
                _s[1] += 1
            _s = [ None if ss is None else ss//b for ss in _s ]
            if flipped:
                if _s[0] == 0:
                    _s = [_s[1]-1, None, -1]
                else:
                    _s = [_s[1]-1, _s[0]-1, -1]
            # Append the new slice
            slices += [slice(*_s)]

        return tuple(slices)

    def get_sec(self):
        # Always assume normal FITS header formatting
        head0 = self.action.args.ccddata.header
        numamps = head0['NVIDINP']

        one_indexed = True
        include_last = True
        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(self.action.args.ccddata.data.shape, dtype=int)
            for i in range(numamps):
                # Get the data section
                sec = head0[section+"{0:1d}".format(i+1)]

                datasec = self.sec2slice(sec, one_indexed=one_indexed,
                                            include_end=include_last, require_dim=2)#, binning=binning)
                # Flip the datasec
                datasec = datasec[::-1]

                # Assign the amplifier
                pix_img[datasec] = i+1

            # Finish
            if section == 'DSEC':
                rawdatasec_img = pix_img.copy()
            elif section == 'BSEC':
                oscansec_img = pix_img.copy()
        return rawdatasec_img, oscansec_img
    
    """
    def get_sec(self):
        # testing red
        head0 = self.action.args.ccddata.header
        
        bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd

        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(self.action.args.ccddata.data.shape, dtype=int)

            if head0['CAMERA'] == 'BLUE':
                numamps = head0['NVIDINP']
                one_indexed = True
                include_last = True
                for i in range(numamps):
                    # Get the data section
                    sec = head0[section+"{0:1d}".format(i+1)]

                    datasec = self.sec2slice(sec, one_indexed=one_indexed,
                                                include_end=include_last, require_dim=2)#, binning=binning)
                    # Flip the datasec
                    datasec = datasec[::-1]

                    # Assign the amplifier
                    pix_img[datasec] = i+1

                # Finish
                if section == 'DSEC':
                    rawdatasec_img = pix_img.copy()
                elif section == 'BSEC':
                    oscansec_img = pix_img.copy()
            
            elif head0['CAMERA']=='RED':
                one_indexed = False
                include_last = True
                for i, ia in enumerate(amps):
                    iac = ia - aoff
                    # Get the data section
                    if section == 'DSEC':
                        data = dsec
                    elif section == 'BSEC':
                        data = bsec
                    
                    if ia <= 1:
                        # this may be wrong for other amplifiers
                        sec = '[{2:d}:{3:d},{0:d}:{1:d}]'.format(*data[iac])
                    else:
                        sec = '[{3:d}:{2:d},{0:d}:{1:d}]'.format(*data[iac])

                    datasec = self.sec2slice(sec, one_indexed=one_indexed,
                                                include_end=include_last, require_dim=2)#, binning=binning)
                    # Flip the datasec
                    datasec = datasec[::-1]

                    # Assign the amplifier
                    pix_img[datasec] = i+1

                # Finish
                if section == 'DSEC':
                    rawdatasec_img = pix_img.copy()
                elif section == 'BSEC':
                    oscansec_img = pix_img.copy()
                

        return rawdatasec_img, oscansec_img
        """


    def pattern_frequency(self, frame, axis=1):
        """
        Using the supplied 2D array, calculate the pattern frequency
        along the specified axis.

        Args:
            frame (`numpy.ndarray`_):
                2D array to measure the pattern frequency
            axis (:obj:`int`, optional):
                Which axis should the pattern frequency be measured?

        Returns:
            :obj:`float`: The frequency of the sinusoidal pattern.
        """
        # For axis=0, transpose
        arr = frame.copy()
        if axis == 0:
            arr = frame.T

        # Calculate the output image dimensions of the model signal
        # Subtract the DC offset
        arr -= np.median(arr, axis=1)[:, np.newaxis]
        # Find significant deviations and ignore those rows
        mad = 1.4826*np.median(np.abs(arr))
        ww = np.where(arr > 10*mad)
        # Create a mask of these rows
        msk = np.sort(np.unique(ww[0]))

        # Compute the Fourier transform to obtain an estimate of the dominant frequency component
        amp = np.fft.rfft(arr, axis=1)
        idx = (np.arange(arr.shape[0]), np.argmax(np.abs(amp), axis=1))

        # Construct the variables of the sinusoidal waveform
        amps = (np.abs(amp))[idx] * (2.0 / arr.shape[1])
        phss = np.arctan2(amp.imag, amp.real)[idx]
        frqs = idx[1]

        # Use the above to as initial guess parameters in chi-squared minimisation
        cosfunc = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * p[1] * xarr + p[2])
        xdata = np.linspace(0.0, 1.0, arr.shape[1])
        # Calculate the amplitude distribution
        amp_dist = np.zeros(arr.shape[0])
        frq_dist = np.zeros(arr.shape[0])
        # Loop over all rows to new independent values that can be averaged
        for ii in range(arr.shape[0]):
            if ii in msk:
                continue
            try:
                popt, pcov = curve_fit(cosfunc, xdata, arr[ii, :], p0=[amps[ii], frqs[ii], phss[ii]],
                                    bounds=([-np.inf, frqs[ii]-1, -np.inf],
                                            [+np.inf, frqs[ii]+1, +np.inf]))
            except ValueError:
                self.logger.warning(f'Input data invalid for pattern frequency fit of row {ii+1}/{arr.shape[0]}')
                continue
            except RuntimeError:
                self.logger.warning(f'Pattern frequency fit failed for row {ii+1}/{arr.shape[0]}')
                continue
            amp_dist[ii] = popt[0]
            frq_dist[ii] = popt[1]
        ww = np.where(amp_dist > 0.0)
        use_amp = np.median(amp_dist[ww])
        use_frq = np.median(frq_dist[ww])
        # Calculate the frequency distribution with a prior on the amplitude
        frq_dist = np.zeros(arr.shape[0])
        for ii in range(arr.shape[0]):
            if ii in msk:
                continue
            try:
                popt, pcov = curve_fit(cosfunc, xdata, arr[ii, :], p0=[use_amp, use_frq, phss[ii]],
                                    bounds=([use_amp * 0.99999999, use_frq-1, -np.inf],
                                            [use_amp * 1.00000001, use_frq+1, +np.inf]))
            except ValueError:
                self.logger.warning(f'Input data invalid for pattern frequency fit of row {ii+1}/{arr.shape[0]}')
                continue
            except RuntimeError:
                self.logger.warning(f'Pattern frequency fit failed for row {ii+1}/{arr.shape[0]}')
                continue
            frq_dist[ii] = popt[1]
        # Ignore masked values, and return the best estimate of the frequency
        ww = np.where(frq_dist > 0.0)
        medfrq = np.median(frq_dist[ww])

        return medfrq/(arr.shape[1]-1)


    def calc_pattern_freq(self, frame, rawdatasec_img, oscansec_img):
        """
        Calculate the pattern frequency using the overscan region that covers
        the overscan and data sections. Using a larger range allows the
        frequency to be pinned down with high accuracy.

        NOTE: The amplifiers are arranged as follows:

        |   (0,ny)  --------- (nx,ny)
        |           | 3 | 4 |
        |           ---------
        |           | 1 | 2 |
        |     (0,0) --------- (nx, 0)

        .. todo::

            PATTERN FREQUENCY ALGORITHM HAS NOT BEEN TESTED WHEN BINNING != 1x1

        Parameters
        ----------
        frame : `numpy.ndarray`_
            Raw data frame to be used to estimate the pattern frequency.
        rawdatasec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            data pixels (0 is no data, non-zero values indicate the amplifier
            number).
        oscansec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            overscan pixels (0 is no data, non-zero values indicate the
            amplifier number).
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file.

        Returns
        -------
        hdu : `astropy.io.fits.HDUList`_
            The input HDUList, with header updated to include the frequency
            of each amplifier.
        """

        # Make a copy of te original frame
        raw_img = frame.copy()

        # Get a unique list of the amplifiers
        unq_amps = np.sort(np.unique(oscansec_img[np.where(oscansec_img >= 1)]))
        num_amps = unq_amps.size

        # Loop through amplifiers and calculate the frequency
        for amp in unq_amps:
            # Grab the pixels where the amplifier has data
            pixs = np.where((rawdatasec_img == amp) | (oscansec_img == amp))
            rmin, rmax = np.min(pixs[1]), np.max(pixs[1])
            # Deal with the different locations of the overscan regions in 2- and 4- amp mode
            if num_amps == 2:
                cmin = 1+np.max(pixs[0])
                frame = raw_img[cmin:, rmin:rmax].astype(np.float64)
            elif num_amps == 4:
                if amp in [1, 2]:
                    pixalt = np.where((rawdatasec_img == amp+2) | (oscansec_img == amp+2))
                    cmin = 1+np.max(pixs[0])
                    cmax = (np.min(pixalt[0]) + cmin)//2  # Average of the bottom of the top amp, and top of the bottom amp
                else:
                    pixalt = np.where((rawdatasec_img == amp-2) | (oscansec_img == amp-2))
                    cmax = 1+np.min(pixs[0])
                    cmin = (np.max(pixalt[0]) + cmax)//2
                frame = raw_img[cmin:cmax, rmin:rmax].astype(np.float64)
            # Calculate the pattern frequency
            freq = self.pattern_frequency(frame)

            self.logger.info("Pattern frequency of amplifier {0:d}/{1:d} = {2:f}".format(amp, num_amps, freq))
            # Add the frequency to the zeroth header
            self.action.args.ccddata.header['PYPFRQ{0:02d}'.format(amp)] = freq

        # Return the updated HDU
        return

    def rect_slice_with_mask(self, image, mask, mask_val=1):
        """
        Generate rectangular slices from a mask image.

        Args:
            image (`numpy.ndarray`_):
                Image to mask
            mask (`numpy.ndarray`_):
                Mask image
            mask_val (:obj:`int`, optional):
                Value to mask on

        Returns:
            :obj:`tuple`: The image at mask values and a 2-tuple with the
            :obj:`slice` objects that select the masked data.
        """
        pix = np.where(mask == mask_val)
        slices = (slice(np.min(pix[0]), np.max(pix[0])+1), slice(np.min(pix[1]), np.max(pix[1])+1))
        return image[slices], slices

    def subtract_pattern(self, rawframe, datasec_img, oscansec_img, frequency=None, axis=1, debug=False):
        """
        Subtract a sinusoidal pattern from the input rawframe. The algorithm
        calculates the frequency of the signal, generates a model, and subtracts
        this signal from the data. This sinusoidal pattern noise was first
        identified in KCWI, but the source of this pattern noise is not currently
        known.

        Args:
            rawframe (`numpy.ndarray`_):
                Frame from which to subtract overscan
            numamplifiers (:obj:`int`):
                Number of amplifiers for this detector.
            datasec_img (`numpy.ndarray`_):
                An array the same shape as rawframe that identifies
                the pixels associated with the data on each amplifier.
                0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
            oscansec_img (`numpy.ndarray`_):
                An array the same shape as rawframe that identifies
                the pixels associated with the overscan region on each
                amplifier.
                0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
            frequency (:obj:`float`, :obj:`list`, optional):
                The frequency (or list of frequencies - one for each amplifier)
                of the sinusoidal pattern. If None, the frequency of each amplifier
                will be determined from the overscan region.
            axis (:obj:`int`, optional):
                Which axis should the pattern subtraction be applied?
            debug (:obj:`bool`, optional):
                Debug the code (True means yes)

        Returns:
            `numpy.ndarray`_: The input frame with the pattern subtracted
        """

        # Copy the data so that the subtraction is not done in place
        frame_orig = rawframe.copy()
        outframe = rawframe.copy()
        tmp_oscan = oscansec_img.copy()
        tmp_data = datasec_img.copy()
        if axis == 0:
            frame_orig = rawframe.copy().T
            outframe = rawframe.copy().T
            tmp_oscan = oscansec_img.copy().T
            tmp_data = datasec_img.copy().T

        # Amplifiers
        amps = np.sort(np.unique(tmp_data[tmp_data > 0])).tolist()

        # Estimate the frequency in each amplifier (then average over all amps)
        if frequency is None:
            frq = np.zeros(len(amps))
            for aa, amp in enumerate(amps):
                pixs = np.where(tmp_oscan == amp)
                #pixs = np.where((tmp_oscan == amp) | (tmp_data ==  amp))
                cmin, cmax = np.min(pixs[0]), np.max(pixs[0])
                rmin, rmax = np.min(pixs[1]), np.max(pixs[1])
                frame = frame_orig[cmin:cmax, rmin:rmax].astype(np.float64)
                frq[aa] = self.pattern_frequency(frame)
            frequency = np.mean(frq)

        # Perform the overscan subtraction for each amplifier
        for aa, amp in enumerate(amps):
            # Get the frequency to use for this amplifier
            if isinstance(frequency, list):
                # if it's a list, then use a different frequency for each amplifier
                use_fr = frequency[aa]
            else:
                # float
                use_fr = frequency

            # Extract overscan
            overscan, os_slice = self.rect_slice_with_mask(frame_orig, tmp_oscan, amp)
            # Extract overscan+data
            oscandata, osd_slice = self.rect_slice_with_mask(frame_orig, tmp_oscan+tmp_data, amp)
            # Subtract the DC offset
            overscan = overscan - np.median(overscan, axis=1)[:, np.newaxis]

            # Convert frequency to the size of the overscan region
            self.logger.info("Subtracting detector pattern with frequency = {0:f}".format(use_fr))
            use_fr *= (overscan.shape[1]-1)

            # Get a first guess of the amplitude and phase information
            amp = np.fft.rfft(overscan, axis=1)
            idx = (np.arange(overscan.shape[0]), np.argmax(np.abs(amp), axis=1))
            # Convert result to amplitude and phase
            amps = (np.abs(amp))[idx] * (2.0 / overscan.shape[1])
            phss = np.arctan2(amp.imag, amp.real)[idx]

            # Use the above to as initial guess parameters in chi-squared minimisation
            cosfunc = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * p[1] * xarr + p[2])
            xdata, step = np.linspace(0.0, 1.0, overscan.shape[1], retstep=True)
            xdata_all = (np.arange(osd_slice[1].start, osd_slice[1].stop) - os_slice[1].start) * step
            model_pattern = np.zeros_like(oscandata)
            val = np.zeros(overscan.shape[0])
            # Get the best estimate of the amplitude
            for ii in range(overscan.shape[0]):
                try:
                    popt, pcov = curve_fit(cosfunc, xdata, overscan[ii, :], p0=[amps[ii], use_fr, phss[ii]],
                                        bounds=([-np.inf, use_fr * 0.99999999, -np.inf], [+np.inf, use_fr * 1.00000001, +np.inf]))
                except ValueError:
                    self.logger.warning("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                    continue
                except RuntimeError:
                    self.logger.warning("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                    continue
                val[ii] = popt[0]
                model_pattern[ii, :] = cosfunc(xdata_all, *popt)
            use_amp = np.median(val)
            # Get the best estimate of the phase, and generate a model
            for ii in range(overscan.shape[0]):
                try:
                    popt, pcov = curve_fit(cosfunc, xdata, overscan[ii, :], p0=[use_amp, use_fr, phss[ii]],
                                        bounds=([use_amp * 0.99999999, use_fr * 0.99999999, -np.inf],
                                                [use_amp * 1.00000001, use_fr * 1.00000001, +np.inf]))
                except ValueError:
                    self.logger.warning("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                    continue
                except RuntimeError:
                    self.logger.warning("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                    continue
                model_pattern[ii, :] = cosfunc(xdata_all, *popt)
            outframe[osd_slice] -= model_pattern

        debug = False
        if debug:
            from astropy.io import fits
            hdu = fits.PrimaryHDU(rawframe)
            hdu.writeto("tmp/tst_raw.fits", overwrite=True)
            hdu = fits.PrimaryHDU(outframe)
            hdu.writeto("tmp/tst_sub.fits", overwrite=True)
            hdu = fits.PrimaryHDU(np.array(rawframe, dtype=float) - np.array(outframe, dtype=float))
            hdu.writeto("tmp/tst_mod.fits", overwrite=True)

        # Transpose if the input frame if applied along a different axis
        if axis == 0:
            outframe = outframe.T
        # Return the result
        return outframe
    

    def _perform(self):

        # Header keyword to update
        key = 'SINESUB'
        keycom = "was sine pattern subtracted?"
        suffix = 'ins'

        out_file = os.path.join(self.config.instrument.output_directory, \
            os.path.basename(self.action.args.name))
        if suffix is not None:
            (main_name, extension) = os.path.splitext(out_file)
        out_file = main_name + "_" + suffix + extension

        # Skip if requested
        if self.action.args.ccddata.header['CAMERA'] == 'RED':
            return self.action.args
        if not self.config.instrument.subtract_sine:
            self.logger.info("Skipping sine-pattern subtraction by request")
            self.action.args.ccddata.header[key] = (False, keycom)

            return self.action.args
        elif os.path.isfile(out_file) and not self.config.instrument.clobber:
            self.logger.info(f"Reading existing file: {out_file}")
            self.action.ccddata = kcwi_fits_reader(out_file)[0]

            return self.action.args
        else:
            # perform

            # locate oscan and raw frame
            rawdatasec_img, oscansec_img = self.get_sec()

            # calc pattern freq
            self.calc_pattern_freq(self.action.args.ccddata.data, rawdatasec_img, oscansec_img)

            # fit and subtract
            frequency = []
            try:
                # Grab a list of all the amplifiers
                amps = np.sort(np.unique(oscansec_img[np.where(oscansec_img > 0)]))
                for amp in amps:
                    frequency.append(self.action.args.ccddata.header['PYPFRQ{0:02d}'.format(amp)])
                # Final check to make sure the list isn't empty (which it shouldn't be, anyway)
                if len(frequency) == 0:
                    frequency = None
            except KeyError:
                frequency = None
            # Subtract the pattern and overwrite the current image
            self.action.args.ccddata.data = self.subtract_pattern(self.action.args.ccddata.data, \
                                                    rawdatasec_img, oscansec_img,
                                                    frequency=frequency)

        log_string = SubtractSinePattern.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out ins image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="ins")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix=suffix,
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        self.logger.info(log_string)

        return self.action.args

    # END: SubtractSinePatter()
