#package for GUI setup
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import scrolledtext
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import sys
import subprocess
import threading
from glob import glob


#package for astro analyis
from glob import glob
from astropy.io import fits
import os 
import kcwi_tools 
import pyregion
import re
from scipy.interpolate import interp1d, splrep, splev
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import List
from pypeit.par import pypeitpar
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph



"""
Author: Zhuyun Zhuang, with the help from chatGPT! (free version!)
Date: 05/2024
Usage: python kcwi_zap_viewer.py

Description: Load data from a directory and allow user to interactively
perform sky subtraction of KCWI/KCRM data with ZAP, and flux calibration. 
"""
####Some user-defined setup. People should update it based on their own needs######
#only set it for test purporse. When it browse the directory, it will start from your favorite directory storing the data :)
initial_dir = '/scr/zzhuang/keck_obs/kcwi/2023sep23/red/redux/noskysub_drp'

#Set it to the place where you put the TelFit file from pypeit. Can download it via "pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits"
telgridfile = '/scr/zzhuang/telluric/TelFit_MaunaKea_3100_26100_R20000.fits'

#pick this region to generate the white-lighted image because sky lines are much stronger elsewhere. 
# TODO: Can also make it as an input or variable parameter in the GUI
wlimg_wave_range = [6380, 7200] 



def telluric_correct(infile_path: str,  star_ra: float, star_dec: float, 
                     spectrograph = 'keck_kcrm'):
    """
    Telluric correct the spectrum of a given standard star.
    Original written by Milan Sharma Mandigo-Stoba for DBSP_DRP at https://github.com/finagle29/DBSP_DRP/blob/main/dbsp_drp/telluric.py;
    Adatpted by Zhuyun Zhuang for the KCWI data

    Args:
        coadd (str): Coadd filename.
        star_ra (float), star_dec (float): the RA and Dec of the STD
        spectrograph (str): PypeIt name of spectrograph.
    """
    spectrograph = load_spectrograph(spectrograph)
    par = spectrograph.default_pypeit_par()

    #only used for standard star
    par['telluric']['objmodel'] = 'star'
    par['telluric']['star_ra'] = star_ra
    par['telluric']['star_dec'] = star_dec

    # if par['telluric']['telgridfile'] is None:
    #     if par['sensfunc']['IR']['telgridfile'] is not None:
    #         par['telluric']['telgridfile'] = par['sensfunc']['IR']['telgridfile']
    par['telluric']['telgridfile'] = telgridfile
    par['telluric']['teltype'] = 'grid'

    # Parse the output filename
    outfile = re.sub('.fits', '_tellcorr.fits', infile_path)
    modelfile = re.sub('.fits', '_tellmodel.fits', infile_path)

    try:
        TelStar = telluric.star_telluric(infile_path, par['telluric']['telgridfile'],
                                        modelfile, outfile,
                                        star_type=par['telluric']['star_type'],
                                        star_mag=par['telluric']['star_mag'],
                                        star_ra=par['telluric']['star_ra'],
                                        star_dec=par['telluric']['star_dec'],
                                        func=par['telluric']['func'],
                                        model=par['telluric']['model'],
                                        polyorder=par['telluric']['polyorder'],
                                        only_orders=par['telluric']['only_orders'],
                                        teltype=par['telluric']['teltype'], tell_npca=par['telluric']['tell_npca'],
                                        mask_hydrogen_lines=par['sensfunc']['mask_hydrogen_lines'],
                                        mask_helium_lines=par['sensfunc']['mask_helium_lines'],
                                        hydrogen_mask_wid=par['sensfunc']['hydrogen_mask_wid'],
                                        delta_coeff_bounds=par['telluric']['delta_coeff_bounds'],
                                        minmax_coeff_bounds=par['telluric']['minmax_coeff_bounds'],
                                        pix_shift_bounds=par['telluric']['pix_shift_bounds'],
                                        maxiter=par['telluric']['maxiter'],
                                        popsize=par['telluric']['popsize'],
                                        tol=par['telluric']['tol'])
    except ValueError:
        print(f"[ERROR] Telluric correction of {os.path.base(infile_path)} FAILED!")


def kcwi_correct_extin(img0, hdr0):
    """Atmospheric extinction correction from official KCWI DRP"""
    img = img0.copy()
    hdr = hdr0.copy()

    # get airmass
    air = hdr['AIRMASS']

    # read extinction data
    full_path = re.sub('py/kcwi_tools.py', 'data/extin/snfext.fits', kcwi_tools.__file__)
    if os.path.exists(full_path):
        hdul = fits.open(full_path)
        exwl = hdul[1].data['LAMBDA']
        exma = hdul[1].data['EXT']
        # get object wavelengths
        sz = img.shape
        dw = hdr['CD3_3']
        w0 = hdr['CRVAL3']
        owls = np.arange(sz[0]) * dw + w0
        # linear interpolation
        exint = interp1d(exwl, exma, kind='cubic', bounds_error=False,
                         fill_value='extrapolate')
        # resample extinction curve
        oexma = exint(owls)
        # convert to flux ratio
        flxr = 10.**(oexma * air * 0.4)
        if len(sz) == 3:
            # apply to cube
            for ix in range(sz[2]):
                for iy in range(sz[1]):
                    img[:, iy, ix] *= flxr
        else:
            # apply to vector
            img *= flxr

        flrmn = np.nanmean(flxr)
        hdr['HISTORY'] = 'kcwi_correct_extin'
        hdr['EXTCOR'] = (True, 'extinction corrected?')
        hdr['AVEXCOR'] = (flrmn, 'average extin. correction (flux ratio)')
        # if logger:
        #     logger.info("Extinction corrected")
        # else:
        print("Extinction corrected")
            
        return img
    else:
        print("Extinction data file (%s) not found!" % full_path)


def reassign_skyseg(x_segments, widths, new_segment, new_width, min_span = 100):
    """
    Reassign the cfwidth to skyseg required by MUSE.  Primarily written by ChatGPT.
    Using skyseg defined by MUSE directly would cause weird fluctuations if a strong emission line of the science object
    falls in the region where the sky lines are sparse. This function is to replace the cfwidth of the
    regions w/ strong emission lines to a narrower value

    Parameters:
    x_segments (arr, with shape of n): An array of segment grid points.
    widths: (arr, with shape of n-1): The cfwidth of each segment.
    new_segment (arr): The new segment as [a1, a2].
    new_width (float): The weight for the new segment.
    min_span (int): the minimial length of each segment

    Returns:
    list of tuple: The updated list of segments with their weights.
    """
    a1, a2 = new_segment

    # Add new segment to the list of segments and sort
    x_segments.extend([a1, a2])
    x_segments = sorted(set(x_segments))

    # Initialize new widths list
    updated_widths = []

    # Assign widths to the new segments
    for i, (start, end) in enumerate(zip(x_segments[:-1], x_segments[1:])):
        if a1 <= start < a2:
            updated_widths.append(new_width)
        else:
            width = widths[min(i, len(widths) - 1)]  # Ensure valid index access
            updated_widths.append(width)

    # Merge adjacent segments and ensure each spans at least min_span
    final_segments = []
    final_widths = []

    i = 0
    while i < len(x_segments) - 1:
        start, end = x_segments[i], x_segments[i + 1]
        width = updated_widths[i]

        # Merge adjacent segments
        while i < len(x_segments) - 1 and end - start < min_span:
            next_start, next_end = x_segments[i + 1], x_segments[i + 2]
            if next_start - end < min_span:
                end = next_end
                width = new_width
                i += 1
            else:
                break

        final_segments.append((start, end))
        final_widths.append(width)
        i += 1

    # Convert final_segments to 1D array
    final_x_segments = []
    for start, end in final_segments:
        if not final_x_segments or final_x_segments[-1] != start:
            final_x_segments.append(start)
        final_x_segments.append(end)

    return final_x_segments, final_widths
        
class KCWIViewerApp:
    def __init__(self, root):
        """
        Initialize the GUI window
        """
        self.root = root
        self.root.title("KCWI Data Viewer")
        self.last_focused_entry = None

        # self.menu = tk.Menu(root)
        # root.config(menu=self.menu)

        # ############## Add a menu ###########
        # self.file_menu = tk.Menu(self.menu, tearoff=0)
        # self.menu.add_cascade(label="File", menu=self.file_menu)

        # self.file_menu.add_command(label="input directory", command=self.browse_input_directory)
        # self.file_menu.add_command(label="output directory", command=self.browse_output_directory)


        # self.menubar = tk.Menu(root)
        # self.filemenu = tk.Menu(self.menubar, tearoff = 0)
        # self.filemenu.add_command(label = 'input directory', command = self.browse_input_directory)
        # self.filemenu.add_command(label = 'output directory', command = self.browse_output_directory)

        ############ Input directory #################
        self.input_dir_label = tk.Label(root, text="Input Directory:")
        self.input_dir_label.grid(row=0, column=0)
        self.input_dir_entry = tk.Entry(root)
        self.input_dir_entry.grid(row=0, column=1)
        self.browse_button = tk.Button(root, text="Browse", command=self.browse_input_directory)
        self.browse_button.grid(row=0, column=2)

        ############# Output directory############
        self.output_dir_label = tk.Label(root, text="Output Directory:")
        self.output_dir_label.grid(row=0, column=3)
        self.output_dir_entry = tk.Entry(root)
        self.output_dir_entry.grid(row=0, column=4)
        self.browse_button = tk.Button(root, text="Browse", command=self.browse_output_directory)
        self.browse_button.grid(row=0, column=5)

        ############# Science index input############
        self.index_label = tk.Label(root, text="Science Frame No.:")
        self.index_label.grid(row=1, column=0)
        self.index_entry = tk.Entry(root)
        self.index_entry.grid(row=1, column=1)
        self.index_entry.bind("<Return>", self.update_index)

        ############# Buttons for increasing and decreasing index for science frame############
        self.increase_button = tk.Button(root, text="Previous", command=self.decrease_index)
        self.increase_button.grid(row=1, column=2)
        self.decrease_button = tk.Button(root, text="Next", command=self.increase_index)
        self.decrease_button.grid(row=1, column=3)

        #############Check box to select whether an off-field sky is used############
        self.use_index2_var = tk.BooleanVar()
        self.use_index2_checkbox = tk.Checkbutton(root, text="Use Off-field Sky Frame No.", variable=self.use_index2_var, command=self.toggle_index2_entry)
        self.use_index2_checkbox.grid(row=1, column=4)

        ############# sky index input############
        self.index2_entry = tk.Entry(root)
        self.index2_entry.grid(row=1, column=5)
        self.index2_entry.bind("<Return>", self.update_index2)

        ############# Naming structure selection - used for loading the data cube for science data
        # self.structure_label = tk.Label(root, text="Science Data Product type:")
        # self.structure_label.grid(row=2, column=0)
        self.structure_var = tk.StringVar()
        self.structure_options = ["icubed", "icubes"]
        self.structure_var.set('icubed')
        self.structure_menu = tk.OptionMenu(root, self.structure_var, *self.structure_options)
        self.structure_menu.grid(row=2, column=0)

        #############  Load raw data (DRP-reduced cube) button############
        self.load_button = tk.Button(root, text="Load Raw Cube", command=lambda: self.load_data('raw'))
        self.load_button.grid(row=2, column=1)

        #############  Save cropped data (good wavelength region) button############
        self.load_button = tk.Button(root, text="Save Cropped Cube", command=self.save_cropped_data)
        self.load_button.grid(row=2, column=2)

        ############# Load cropped data button############
        self.load_crop_button = tk.Button(root, text="Load Cropped Cube", command=lambda: self.load_data('cropped'))
        self.load_crop_button.grid(row=2, column=3)

        ############# input the redshift of a given source #########
        self.redshift_label = tk.Label(root, text = 'Redshift:')
        self.redshift_label.grid(row = 2, column =4)
        self.redshift_entry = tk.Entry(root)
        self.redshift_entry.grid(row = 2, column = 5)
        self.redshift_entry.bind("<Return>", self.update_redshift)


        ############# Input the mask frame No. (can be different from the science one since we usally take the at least three frames) ############
        # This function is only used for convenience, so that for each pair we only need to create one 
        self.mask_index_label = tk.Label(root, text="ZAP Mask Frame No.:")
        self.mask_index_label.grid(row=3, column=0)
        self.mask_entry = tk.Entry(root)
        self.mask_entry.grid(row=3, column=1)
        self.mask_entry.bind("<Return>", self.update_mindex)
        self.update_mask_button = tk.Button(root, text = 'Update ZAP mask', command=self.update_zap_mask)
        self.update_mask_button.grid(row =3, column = 2)

        #############Check box to select whether multiple skysegment used, and whether have an additional sky seg near Halpha############
        self.use_multi_skyseg = tk.BooleanVar()
        self.use_multi_skyseg_checkbox = tk.Checkbutton(root, text="Use multiple skyseg in ZAP", variable=self.use_multi_skyseg, 
                                                        onvalue = True, offvalue = False)
        self.use_multi_skyseg_checkbox.grid(row=3, column=3)
        self.use_multi_skyseg.set(True)
        self.use_Ha_seg = tk.BooleanVar()
        self.use_Ha_seg_checkbox = tk.Checkbutton(root, text='Additional Sky Seg near Halpha', variable=self.use_Ha_seg,
                                                  onvalue = True, offvalue = False)
        self.use_Ha_seg_checkbox.grid(row = 3, column = 4)
        self.use_Ha_seg.set(True)


        #############Run button for the ZAP ############
        self.run_zap_button = tk.Button(root, text = 'Run ZAP', command = self.run_zap_precondition)
        self.run_zap_button.grid(row =3, column =5)
        

        ############# Setup the flux calibration part for stds ###############
        self.stddir = re.sub('py/kcwi_tools.py', 'data/stds', kcwi_tools.__file__)  #the base directory to read in the flux-calibrated spectrum of a given std
        self.std_index_label = tk.Label(root, text="Standard star invsens (DRP): ")
        self.std_index_label.grid(row = 4, column = 0)
        self.std_entry = tk.Entry(root)
        self.std_entry.grid(row = 4, column = 1)
        self.load_std_button = tk.Button(root, text = 'Browse DRP invsens', command = lambda: self.load_invsens('DRP'))
        self.load_std_button.grid(row = 4, column = 2)
        self.load_std_update_button = tk.Button(root, text = 'Browse updated invsens', command = lambda: self.load_invsens('updated'))
        self.load_std_update_button.grid(row = 4, column = 3)
        self.save_std_button = tk.Button(root, text = 'Save updated invsens', command = self.save_updated_invsens)
        self.save_std_button.grid(row = 4, column = 4)
        # self.std_bspline = tk.Label(root, text = 'B-Spline pars [breakpoints, polyorder]:')
        # self.std_bspline.grid(row = 4, column = 3)
        # self.std_bspline_entry = tk.Entry(root)
        # self.std_bspline_entry.grid(row = 4, column = 5)
        # self.std_bspline_entry.bind("<Return>", self.update_bspline_pars)

        # Plotting area
        self.figure = Figure(figsize = (12, 5))
        self.ax = self.figure.add_subplot(111)
        # self.figure.tight_layout()
        # Global font settings
        plt.rc('font', family='serif', size=12)

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)
        self.canvas.get_tk_widget().focus_force()


        # Add navigation toolbar
        self.toolbarFrame = tk.Frame(master=root)
        self.toolbarFrame.grid(row=6,column=0, columnspan = 6)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        # self.toolbar.home_toggle(True)

        #add the interactive 
        # self.canvas.get_tk_widget().bind('<Button-1>', lambda event: self.canvas.focus_set()) 
        self.canvas.get_tk_widget().focus_set()#only activate the interactive mode when left click the left canvas first
        self.canvas.mpl_connect('button_press_event', self.on_click) #add all the interactive functions with mouse button 
        self.canvas.mpl_connect('key_press_event', self.on_type) #add all the interactive functions with keyboard press

        # Initialize indices
        self.index = -1 #science frame number
        self.index2 = -1 #sky frame number; if specified
        self.redshift = 0.0
        self.index2_entry.config(state=tk.DISABLED)
        self.update_index_entries()

        # Output text widget
        # self.output_text = tk.Text(root, height=5, width=100, font=("Arial", 12))
        self.output_text = scrolledtext.ScrolledText(root, height=7, width=140, font=("Arial", 12), wrap="word")
        self.output_text.grid(row=7, column=0, columnspan=6)
        sys.stderr = RedirectText(self.output_text)
        sys.stdout = RedirectText(self.output_text)
        self.redirect_text = RedirectText(self.output_text)
        self.output_text.bind("<Return>", self.run_zap)

        #text widget for adjusting sky segment
        # self.output_text = scrolledtext.ScrolledText(root, height = 5, width = 40, font = ('Arial', 12))
        # self.output_text.grid(row =7, column = 4, columnspan = 2)
        # self.output_text.bind('<Return>', self.ontype_skyseg)
        # self.skyseg_input_active = False

        ######### final configuration
        # self.menubar.add_cascade(label="File", menu=self.filemenu)
        # root.config(menu = self.menubar)
        root.columnconfigure(1, minsize=100, weight=0)

    def insert_text(self, text, newline = True):

        #finish this line
        if newline:
            self.output_text.insert(tk.END, text + '\n')

        #allow other texts to be input into the current line
        else:
            self.output_text.insert(tk.END, text)
        self.output_text.see(tk.END)

    def focus_in(self, event):
        """
        A variable to keep track of which entry widget was last focused, so that the code would only update one variable at one time.
        """
        self.last_focused_entry = event.widget


    def toggle_index2_entry(self):
        """Decide whether index2 (off-field sky frame should be loaded)"""

        self.index2_entry.config(state=tk.NORMAL if self.use_index2_var.get() else tk.DISABLED)
        #when the off-field sky option is turned off
        if not self.use_index2_var.get():
            self.index2_entry.delete(0, tk.END)
            self.index2 = -1
        else:
            self.index2_entry.delete(0, tk.END)

    def increase_index(self):
        """"Increase the science frame number"""
        self.index += 1
        self.update_index_entries()

    def decrease_index(self):
        """Decrease the science frame number"""
        self.index = max(0, self.index - 1)
        self.update_index_entries()

    def update_index_entries(self):
        #if update the science frame number
        if self.last_focused_entry == self.index_entry and self.index > 0:
            self.index_entry.delete(0, tk.END)
            self.index_entry.insert(tk.END, str(self.index))
            self.insert_text(f"[INFO] Set the science frame: {self.prefix}_{self.index:05d}")
            self.index_entry.selection_clear()

        #if update the sky frame number
        if self.last_focused_entry == self.index2_entry and self.index2 > 0:
            self.index2_entry.delete(0, tk.END)
            self.index2_entry.insert(tk.END, str(self.index2))
            self.insert_text(f"[INFO] Set the sky frame for science frame {self.prefix}_{self.index:05d}: {self.prefix}_{self.index2:05d} ")
            self.index2_entry.selection_clear()

        #if update the mask frame number
        if self.last_focused_entry == self.mask_entry:
            self.mask_entry.delete(0, tk.END)
            self.mask_entry.insert(tk.END, str(self.mindex))
            self.insert_text(f"[INFO] Set the ZAP mask frame for science frame {self.prefix}_{self.index:05d}: {self.prefix}_{self.mindex:05d} ")

        #update the redshift
        if self.last_focused_entry == self.redshift_entry:
            self.redshift_entry.delete(0, tk.END)
            self.redshift_entry.insert(tk.END, str(self.redshift))
            self.insert_text(f"[INFO] Set the redshift to z = {self.redshift:0.3f} ")

         # Set focus back to the master window to allow canvas to capture keys
        self.root.focus_set()


    def update_index(self, event):
        try:
            self.index = int(self.index_entry.get())
        except ValueError:
            pass
        self.update_index_entries()

    def update_index2(self, event):
        """
        Update the sky frame number if the "Use off-field sky" box is selected
        """
        try:
            self.index2 = int(self.index2_entry.get()) #get the sky frame number
        except ValueError:
            pass #no sky frame
        self.update_index_entries()

    def update_mindex(self, event):
        """
        Update the ZAP mask frame number
        """
        try:
            #get the mask frame number        
            self.mindex = int(self.mask_entry.get())
        except ValueError:
            self.insert_text(f"[ERROR] Need to specify the frame used for generating the object mask!")
        self.update_index_entries()

    def update_redshift(self, event):
        """
        Update the input redshift
        """
        try:
            #get the mask frame number        
            self.redshift = float(self.redshift_entry.get())
        except ValueError:
            self.redshift = 0.0
            # self.insert_text(f"[INFO] Redshift not set. Adopt z = 0")
        self.update_index_entries()

    def browse_input_directory(self):
        """Select the input directory"""

        #specify the input directory
        directory = filedialog.askdirectory(initialdir=initial_dir)
        self.input_dir_entry.delete(0, tk.END)
        self.input_dir_entry.insert(tk.END, directory)
        self.base = self.input_dir_entry.get()
        self.insert_text(f"[INFO] Loading the data from {self.base}")

        #find the common string name for a given date
        self.prefix = os.path.basename(glob(f'{directory}/kr*.fits')[0])[:8]

    def browse_output_directory(self):
        """Select the output directory"""
        directory = filedialog.askdirectory()
        self.output_dir_entry.delete(0, tk.END)
        self.output_dir_entry.insert(tk.END, directory)
        self.output = self.output_dir_entry.get()
        # print(f'[INFO] Set the output directory to {self.output}')
        self.insert_text(f"[INFO] Set the output directory to {self.output}")

    def load_data(self, datatype):
        """
        Load the DRP-reduced data cube
        """
        #get the cube type
        ctype = self.structure_var.get()  

        #determine where to load the data
        if datatype == 'raw':
            base = self.base
        elif datatype == 'cropped':
            base = self.output
        else:
            self.insert_text(f'[ERROR] datatype not recognized! Need to be "raw" or "cropped" ')

        # load the science frame
        if self.index > 0:
            self.scihdu = fits.open(f'{base}/{self.prefix}_{self.index:05d}_{ctype}.fits')
            self.objname = self.scihdu[0].header['OBJECT']
            if datatype == 'raw':
                self.scihdu = self.crop_cube(self.scihdu) #crop the data cube to good wavelength region
                self.insert_text(f"[INFO] Loading the DRP-reduced science frame {self.prefix}_{self.index:05d}")
            elif datatype == 'cropped':
                self.insert_text(f"[INFO] Loading the cropped DRP-reduced science frame {self.prefix}_{self.index:05d}")

            self.scihdr = self.scihdu[0].header
            self.obswave = (np.arange(self.scihdr['NAXIS3']) + 1 - self.scihdr['CRPIX3']) * self.scihdr['CD3_3'] + self.scihdr['CRVAL3']

        else:
            self.insert_text(f"[ERROR] Wrong science frame! Need to set it to a positive integer. Check the KCWI log for the frame number!")

        # self.z = 0
        if self.index2 > 0:
            self.skyhdu = fits.open(f'{base}/{self.prefix}_{self.index2:05d}_{ctype}.fits')
            if datatype == 'raw':
                self.skyhdu = self.crop_cube(self.skyhdu) #crop the data cube to good wavelength region
                self.insert_text(f"[INFO] Loading the DRP-reduced sky frame {self.prefix}_{self.index2:05d}")
            else:
                self.insert_text(f"[INFO] Loading the cropped DRP-reduced sky frame {self.prefix}_{self.index:05d} for {self.prefix}_{self.index:05d}")
            self.plot_spectrum(self.scihdu) #plot the spectrum of the central region (x = [12, 22], y = [43, 53]; 10x10 box) for a quick look. 
        else:
            self.skyhdu = None
            self.plot_spectrum(self.scihdu, hdu_sky=self.skyhdu) #plot the spectrum of the central region (x = [12, 22], y = [43, 53]; 10x10 box) for a quick look. 

    def save_cropped_data(self):
        """Save the cropped datacube within the good wavelength region, white-lighted image, and preliminary mask for ZAP"""

        #for science frame
        hdu_list = []
        indices = []
        if self.index > 0:
            hdu_list.append(self.scihdu)
            indices.append(self.index)
        if self.index2 > 0:
            hdu_list.append(self.skyhdu)
            indices.append(self.index2)
        
        if len(hdu_list) > 0:
            for i, hdu in  enumerate(hdu_list):
                index = indices[i]

                self.insert_text(f"[INFO] Saving the cropped datacube for {self.prefix}_{index:05d}")

                #cropped datacube
                hdu.writeto(f'{self.output}/{self.prefix}_{index:05d}_icubed.fits', overwrite = True)

                #white-lighted image
                wlimg_index = np.where((self.obswave >= wlimg_wave_range[0]) & (self.obswave <= wlimg_wave_range[1]))[0]
                wlimg = np.sum(hdu[0].data[wlimg_index], axis = 0)
                hdr2d = kcwi_tools.collapse_header(self.scihdu[0].header)
                wlhdu = fits.PrimaryHDU(wlimg, header = hdr2d)
                mask = np.mean(self.scihdu['FLAGS'].data, axis = 0)
                mhdu = fits.ImageHDU(mask, header = hdr2d)
                hdulist = fits.HDUList([wlhdu, mhdu])
                hdulist.writeto(f'{self.output}/{self.prefix}_{index:05d}_wlimg.fits', overwrite = True)

                #preliminary ZAP mask - only for the sky cube
                if i > 0:
                    edgemask = mask > 1
                    allmask = np.zeros_like(mask)
                    allmask[edgemask] = 1
                    maskhdu = fits.PrimaryHDU(allmask, header = hdr2d)
                    maskhdu.writeto(f'{self.output}/{self.prefix}_{index:05d}_zap_mask.fits', overwrite = True)

    def update_zap_mask(self):
        """Update the ZAP mask for a given science frame"""

        mindex = self.mindex

        mhdu = fits.open(f'{self.output}/{self.prefix}_{self.index:05d}_wlimg.fits')
        edgemask = mhdu[1].data > 1

        region_path = f'{self.output}/{self.prefix}_{mindex:05d}.reg'
        if not os.path.exists(region_path):
            self.insert_text(f"[ERROR] Region file {self.prefix}_{mindex:05d}.reg not exists!")
        else:
            self.insert_text(f"[INFO] Reading region file {self.prefix}_{mindex:05d}.reg for {self.prefix}_{self.index:05d}")
            r = pyregion.open(region_path)
            region_mask  = r.get_mask(hdu=mhdu[0])
            allmask = np.zeros_like(mhdu[0].data)
            allmask[edgemask | region_mask] = 1
            maskhdu = fits.PrimaryHDU(allmask, header = mhdu[0].header)

            #save ZAP mask
            filename = f'{self.output}/{self.prefix}_{self.index:05d}_zap_mask.fits'
            if os.path.exists(filename):
                confirm = messagebox.askyesno("File Exists",
                                              f'ZAP mask {self.prefix}_{self.index:05d}_zap_mask.fits already exists. Do you want to overwrite it?')
                if not confirm:
                    return
            self.insert_text(f"[INFO] Saving ZAP mask to {self.prefix}_{self.index:05d}_zap_mask.fits")
            maskhdu.writeto(filename, overwrite = True)

    
    def crop_cube(self, hdu):
        """Crop the HDUList to only consist of the region"""
        
        hdr = hdu[0].header
        wave = (np.arange(hdr['NAXIS3']) + 1 - hdr['CRPIX3']) * hdr['CD3_3'] + hdr['CRVAL3']
        wavegood_idx = np.where((wave >= hdr['WAVGOOD0']) & (wave <= hdr['WAVGOOD1']))[0]
        
        for h in hdu:
            h.data = h.data[wavegood_idx]
        hdu[0].header['CRVAL3'] = wave[wavegood_idx][0]
        return hdu


    def plot_spectrum(self, hdu, z = 0, xrange = [12, 22], yrange = [43, 53], hdu_sky = None, restore_limit = False):
        """
        plot the spectrum of a given region
        
        Args:
            hdu (fits.HDUList) of a given data cube
            z (float): the redshift of the object
            xrange (list): the start and ending of the x-axis from which the spec is extracted
            yrange (list): the start and ending of the x-axis from which the spec is extracted
        """
        #restore the limit of the original plot, so that the plot would stay in the zoomed-in version
        if restore_limit:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

        self.ax.clear()

        spec = np.sum(hdu[0].data[:, yrange[0]:yrange[1], xrange[0]: xrange[1]], axis = (1,2))
        err = np.sqrt(np.sum(hdu['UNCERT'].data[:, yrange[0]:yrange[1], xrange[0]: xrange[1]]**2, axis = (1,2)))
        self.ax.step(self.obswave / (1+z), spec, color ='k', lw = 1, label = 'sci spec')
        self.ax.step(self.obswave / (1+z), err, color ='r', lw = 1, label = 'sci err')
        if hdu_sky is not None:
            skyspec = np.sum(hdu_sky.data[:, yrange[0]:yrange[1], xrange[0]: xrange[1]], axis = (1,2))
            self.ax.step(self.obswave / (1+z), skyspec, color ='lightskyblue', lw = 1, label = 'sky spec')
        self.ax.legend()
        self.ax.set_title(f'{self.prefix}_{self.index:05d}  -  {self.objname}')
        # self.figure.tight_layout()
        self.ax.set_xlabel('Wavelength (A)')
        self.ax.set_ylabel('Flux (Electrons)')

        if restore_limit:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
            
        self.canvas.draw()


    def run_zap_precondition(self):
        """
        Initialize the ZAP setup, used for the connection with ZAP button
        """

        self.zap = {}

        ### Only use one sky segment
        if not self.use_multi_skyseg.get():
            self._zap_skyseg = [] #use one segment, default setting of ZAP
            self._zap_cfwidth = 300 #cfwidth = 300, default setting of ZAP

        else:
            #Use the sky segment from the old version of MUSE. See Table 1 in Soto+16 for details 
            self._zap_skyseg = [0, 5400, 5850, 6440, 6750, 7200, 7700, 8265, 8602, 8731, 9275, 10000] 
            self._zap_cfwidth = [300] * (len(self._zap_skyseg) - 1) #cfwidth = 300, default setting of ZAP

            if self.use_Ha_seg.get():
                print('yes')
                Ha_comp = [6540, 6595] #spectral range of the Halpha+NII complex. 
                self._zap_skyseg, self._zap_cfwidth = reassign_skyseg(self._zap_skyseg, self._zap_cfwidth, 
                                                        (Ha_comp[0]*(1+self.redshift), Ha_comp[1]*(1+self.redshift)), 5)
        
        #store the default skyseg and cfwidth in case we need to restore them later
        # self._zap_skyseg = skyseg.copy()
        # if np.isscalar(cfwidth):
        #     self._zap_cfwidth = cfwidth
        # else:
        #     self._zap_cfwidth = cfwidth.copy()
        self.zap['skyseg'] = self._zap_skyseg.copy()
        self.zap['cfwidth'] = self._zap_cfwidth

        #print the default sky segment                                                
        # self.output_text.delete(1.0, tk.END) #clear the widget
        self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])

        self.skyseg_input_active = True
        # self.insert_text(f"Type 'q' if satisfied with the segment, type 'r' to reset the segmengt, type 'd' to set your own sky segment, or type the new segment as 'start end cfwidth' [e.g., 7000 7200 50]\n: ", newline = False)
        # self.insert_text(f"Type 'q' if satisfied with the segment")
        # self.insert_text(f"Type 'r' to reset the default segmengt")
        # self.insert("Type the new segment as 'start end cfwidth' [e.g., 7000 7200 50] if you want to update the sky segment")
        # self.insert("If you want to define your own sky segment, please type the new setting as [segment]; [cfwidth] (e.g., [7000, 8000, 9000]; [50, 60]) ")
        # self.insert(': ', newline = False)

        # self.ontype_skyseg_update()

    def run_zap(self, event):
        """
        Main program to run ZAP
        """
        if self.skyseg_input_active:
            # #get the content of the last line
            line = self.output_text.get("end-1c linestart", "end-1c")

            #the reset the sky segment
            if 'r' in line:
                self.zap['skyseg'] = self._zap_skyseg
                self.zap['cfwidth'] = self._zap_cfwidth
                print(len(self.zap['skyseg']), len(self.zap['cfwidth']))
                print(self._zap_skyseg)
                print(self._zap_cfwidth)
                self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
               # textwidget.delete("end-1c linestart", "end")

                # self.output_text

            #the user-defined sky segment
            elif '[' in line and ']' in line and ';' in line:
                #load the content of the last line
                line = self.output_text.get("end-1c linestart", "end-1c")
                line = re.sub('[\[:,\]]','', line)
                print(line)
                seg, w = line.split(';')

                try:
                    seg = [float(n) for n in seg.split()]
                    cfwidth = [float(n) for n in w.split()]
                    if len(seg) - 1 == len(cfwidth):
                        if np.abs(seg[0]) > 1e-6 or np.abs(seg[-1] - 10000) > 1e-6:
                            self.insert('Please remeber to include 0 and 10000 as the first and last element of your sky segment list!')
                        self.zap['skyseg'] = seg
                        self.zap['cfwidth'] = cfwidth
                        self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])

                    else:
                        self.insert_text(f"\nInvalid input!\n: ", newline = False)

                        # self.insert_text(f"Invalid input! Please input the new sky segments and cfwidths! The format need to be similar as [7000, 8000, 9000]; [50, 60]. The first one is new segment array.\n: ", newline = False)

                except:
                    self.insert_text(f"Invalid input!\n: ", newline = False)

                    # self.insert_text(f"Invalid input! Please input the new sky segments and cfwidths! The format need to be similar as [7000, 8000, 9000]; [50, 60]. The first one is new segment array.\n: ", newline = False)


            #take the new input of sky segment
            elif 'q' not in line:
                #load the content of the last line
                line = self.output_text.get("end-1c linestart", "end-1c")
                
                try:
                    #read three numbers
                    inputs = [float(n) for n in line[1:].split()]
                    self.zap['skyseg'], self.zap['cfwidth'] = reassign_skyseg(self.zap['skyseg'], self.zap['cfwidth'], 
                                                            (inputs[0], inputs[1]), inputs[2])
                    self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
                    # self.insert_text(f"Type 'q' if satisfied with the segment, type 'r' to reset the segmengt, or type the new segment as 'start end cfwidth' [e.g., 7000 7200 50].\n: ", newline = False)
                except:
                    self.insert_text(f"Invalid input!\n: ", newline = False)
            
            #running ZAP
            elif 'q' in line:
                self.insert_text('[INFO] Running ZAP...')
                self.skyseg_input_active = False

            # Prevent the default behavior of the Return key
            return "break"
            
        
    
    def print_sky_seg(self, skyseg, cfwidth = 300):
        """
        print the sky segment into the text widget
        Args:
            skyseg: If not using multiple sky segment for the PCA analysis, it is set as [] by default (same as the official MUSE ZAP);
                    If using multiple sky segment, it takes a list of the endpoints in each interval (including the start and the end wavelength)  - N elements          
            cfwidth: A scalar (for one sky segment) or a list with (N-1 elements representing the cfwidth in each segment)
        """
        nseg = len(skyseg)

        self.insert_text('\n ##########Sky segment for ZAP##########')

        #one cfwidth all segments
        if np.isscalar(cfwidth):
            if nseg == 0:
                self.insert_text('(0, 10000) cfwidth = %i'%cfwidth)
                self.insert_text(f'Using 1 sky segment!')

            else:
                for i in range(nseg - 1):
                    self.insert_text(f'(%i, %i) cfwidth = %i'%(skyseg[i], skyseg[i+1], cfwidth))
                self.insert_text(f'Using {nseg:d} sky segment!')


        #multiple sky segments
        else:
            if nseg - 1 != len(cfwidth):
                self.insert_text(f'Length mismatch! The length of cfwidth needs to be N-1 for N sky segments!')
                return

            for i in range(nseg - 1):
                self.insert_text(f'(%i, %i) cfwidth = %i'%(skyseg[i], skyseg[i+1], cfwidth[i]))
            self.insert_text(f'Using {nseg:d} sky segment!')
        
        if nseg == 0:
            ####TODO get the single sky element case done

        #multiple sky segment
        else:
            lines = (
                    "Type 'q' if satisfied with the segment, " 
                    "type 'r' to reset the default segmengt, "
                    "type the new segment as 'start end cfwidth' [e.g., 7000 7200 50] if you want to update the sky segment, "
                    "or if you want to define your own sky segment, please type the new setting as '[list of segment]; [list of cfwidth]' (e.g., [0, 7000, 10000]; [50, 60]).")
            # self.insert_text(f"Type 'q' if satisfied with the segment, type 'r' to reset the default segmengt")
            # self.insert_text(f"Type 'r' to reset the default segmengt")
            # self.insert_text("Type the new segment as 'start end cfwidth' [e.g., 7000 7200 50] if you want to update the sky segment")
            # self.insert_text("If you want to define your own sky segment, please type the new setting as '[list of segment]; [list of cfwidth]' (e.g., [7000, 8000, 9000]; [50, 60]) ")
            self.insert_text(lines)
            self.insert_text(': ', newline = False) 

            
    def load_invsens(self, type):
        """
        Load the  _invsens.fits of a standard star.
        """

        filename = filedialog.askopenfilename(initialdir=self.output) #REAL code
        # filename = '/scr/zzhuang/keck_obs/kcwi/2023sep23/red/redux/kr230923_00174_invsens.fits' #for test purpose only
        self.std_entry.delete(0, tk.END)
        self.std_entry.insert(tk.END, filename)

        #load the invsens file 
        hdu = fits.open(self.std_entry.get())
        self.std = {}

        hdr = hdu[0].header
        wvl = (np.arange(hdr['NAXIS1']) + 1 - hdr['CRPIX1']) * hdr['CDELT1'] + hdr['CRVAL1']
        #Add 3A on each side to give the spline-fit some buffer 
        good = np.where((wvl >= hdr['WAVGOOD0'] - 3) & (wvl<= hdr['WAVGOOD1'] + 3))[0]
        self.std['invsens_hdr'] = hdr

        #raw invsens from the DRP
        if type == 'DRP':
            self.std['wave'] = wvl[good]
            #the first row is the raw ratio of real flux and count to be fitted for; second row the best-fit invsens; third row is the raw electron/s 
            self.std['invsens_data'] = hdu[0].data[0, good]
            self.std['invsens_model_drp'] = hdu[0].data[1, good]
            self.std['counts'] = hdu[0].data[2, good]
            self.std['name'] = hdr['OBJECT'].lower()
            self.std['frame'] = re.sub('_invsens.fits', '', os.path.basename(self.std_entry.get()))
            flag = np.full(len(self.std['wave']), 1, dtype = int)
            self.std['flag'] = self.mask_skyline_region(self.std['wave'], flag)

            #sens func and telluric model to be updated
            self.std['invsens_model'] = None
            self.std['tellmodel'] = None

        #the updated version, so no need to crop the data
        elif type == 'updated':
            self.std['wave'] = wvl
            #the first row is the raw ratio of real flux and count to be fitted for; second row the best-fit invsens; third row is the raw electron/s 
            self.std['invsens_data'] = hdu[0].data[0]
            self.std['invsens_model_drp'] = hdu[0].data[1]
            self.std['counts'] = hdu[0].data[2]
            self.std['name'] = hdr['OBJECT'].lower()
            self.std['invsens_model'] = hdu[0].data[3]
            self.std['flag'] = hdu[0].data[4]
            if len(hdu[0].data) > 5:
                self.std['tellmodel'] = hdu[0].data[5]
            else:
                self.std['tellmodel'] = None
            self.std['frame'] = re.sub('_invsens_updated.fits', '', os.path.basename(self.std_entry.get()))



        #setup the B-Spline fit parameters
        self.std['bspline_bkpt'] = 150 #breakpoints
        self.std['bspline_polyorder'] = 3 #polynomial order between interval
        
        self.insert_text(
                                f"[INFO] Loading the {self.std_entry.get()}")


        #load the flux-calibrated spec for comparison
        try:
            std = fits.getdata('{0}/{1}.fits'.format(self.stddir, self.std['name']))
        except ValueError:
            self.insert_text(f"Cannot find the std spec for {self.std_name} in {self.stddir}")

        use = np.where((std['WAVELENGTH'] >= hdr['WAVGOOD0'] - 5) & (std['WAVELENGTH'] <= hdr['WAVGOOD1'] + 5))[0] #only load the spectrum in the good wavelength region
        self.std['spec_calib'] = np.column_stack((std['WAVELENGTH'][use], std['FLUX'][use]))

        # self.std['flag'] = self.mask_skyline_region(self.std['wave']) #flag indicated if a region is masked out
        # self.std['use_ind'] = np.array([], dtype = int) #a list of indices for the single points used for fitting
        #Initialize the region for later interactive selection

        self.std['region_start'] = None
        self.std['region_end'] = None

        #set the focuse to the canvas page
        self.canvas.get_tk_widget().focus_set()

        self.plot_std()

    def save_updated_invsens(self):

        """
        Save the updated invsens file to the output directory
        """
        
        ##########save the updated invsens to the output directory
        newhdr = self.std['invsens_hdr'].copy()
        newhdr['NAXIS1'] = len(self.std['wave'])
        newhdr['CRVAL1'] = self.std['wave'][0]
        #still similar to the format of the original invsens.fits; add the updated invsens model and flag to the last two rows
        newdata = np.vstack((self.std['invsens_data'], self.std['invsens_model_drp'], self.std['counts'],
                             self.std['invsens_model'], self.std['flag']))

        #add the tellmodel @ AM =1 to the file
        if self.std['tellmodel'] is not None:
            newdata = np.vstack((newdata, self.std['tellmodel']))

        newhdu = fits.PrimaryHDU(newdata, header = newhdr)

        frame = self.std['frame']
        filename = f'{self.output}/{frame}_invsens_updated.fits'
        if os.path.exists(filename):
            confirm = messagebox.askyesno("File Exists",
                                            f'Updated invsens file {frame}_invsens_updated.fits already exists. Do you want to overwrite it?')
            if not confirm:
                return
        self.insert_text(f"[INFO] Saving the updated invsens file to {frame}_invsens_updated.fits")
        newhdu.writeto(filename, overwrite = True)
        

    def plot_std(self, restore_limit = False):
        """
        Plot the best-fit STD 
        """
        #restore the limit of the original plot, so that the plot would stay in the zoomed-in version
        if restore_limit:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
        
         #plot the flux-calibrated spec
        self.ax.clear()
        self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model_drp'], 
                        color = 'r', lw =1, label = 'best-fit model from DRP', where = 'mid', alpha = 0.5) #raw count x invsens = flux-calibrated spec
        self.ax.plot(self.std['spec_calib'][:,0], self.std['spec_calib'][:,1], color = 'k', label = 'flux-calibrated spec of std')


        use_region = np.where(self.std['flag'] == 1)[0]
        use_point = np.where(self.std['flag'] == 2)[0]

        #plot the updated flux-calibrated model
        if self.std['invsens_model'] is not None:
            self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model'], color = 'cyan', 
                         lw =1, label = 'best-fit new model', where = 'mid') #raw count x invsens = flux-calibrated spec
            self.ax.plot(self.std['wave'][use_region], (self.std['counts'] * self.std['invsens_model'])[use_region], 'x', 
                        color = 'lightgreen', ms = 5, label = 'selected for fitting') #selected regions 
            if len(use_point) > 0:
                self.ax.plot(self.std['wave'][use_point], (self.std['counts'] * self.std['invsens_model'])[use_point], 'o', 
                            color = 'darkgreen', ms = 10, label = 'selected for fitting') #selected pixels 

             #plot the telluric-corrected model
            if self.std['tellmodel'] is not None:
                telluric = self.std['tellmodel']**(self.std['invsens_hdr']['AIRMASS']) #convert the model at AM=1.0 to the real AM
                self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model'] / telluric, color = 'royalblue',
                            where = 'mid', lw = 1, label = 'telluric-corrected std spec')
        
        #plot the DRP-reduced flux-calibrated model
        else:
            self.ax.plot(self.std['wave'][use_region], (self.std['counts'] * self.std['invsens_model_drp'])[use_region], 'x', 
                    color = 'lightgreen', ms = 5, label = 'selected for fitting') #selected regions 
            if len(use_point) > 0:
                self.ax.plot(self.std['wave'][use_point], (self.std['counts'] * self.std['invsens_model_drp'])[use_point], 'o', 
                            color = 'darkgreen', ms = 10, label = 'selected for fitting') #selected pixels 


        self.ax.set_title('%s - %s'%(self.std['frame'], self.std['name']))
        self.ax.set_xlabel('Wavelength (A)')
        self.ax.set_ylabel('flam')
        self.ax.legend()

        if restore_limit:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        self.canvas.draw()
        

    def mask_skyline_region(self, wave, flag):
        """
        Mask out the region with dense sky line and telluric absorption from flux calibration fitting
        """
        regions = [[6274, 6302],[6864.00, 6950.00], [7160.00, 7385.00], [7590, 7691],[8102, 8375],
                   [8943, 9225], [9300, 9400],
                   [6514.00, 6600.00] #Halpha region
                   ]

        mask = False
        for r in regions:
            mask |= (wave >= r[0]) & (wave <= r[1])

        flag[mask] = 0

        return flag
        # [6864, 6935], [7164, 7345], [7591, 7694], [8131]]

    def run_telluric(self): 

        """
        Run the telluric correction on the standard star
        """

        if not hasattr(self, 'std'):
            self.insert_text(f"[ERROR] The updated invsens file not exists! Run the sensfunc fit first with the key'f' first!")
            return
        
        ################# create the "fake" Pypeit output for telluric correction ############
        frame = self.std['frame']
        self.insert_text(f"[INFO] Running telluric correction on {frame}...")

        newhdr = self.std['invsens_hdr'].copy()
        newhdr['NAXIS1'] = len(self.std['wave'])
        newhdr['CRVAL1'] = self.std['wave'][0]
        std_flux = self.std['counts'] * self.std['invsens_model']  / 1e-17 #pypeit takes the unit as 1e-17 erg/s/cm2/A
        SNR = 1000 #fake SNR of STD (required by Pypeit input)
        std_ivar = 1. / (std_flux / SNR)**2
        col1 = fits.Column(name = 'wave', format = '1D', array = self.std['wave'])
        col2 = fits.Column(name = 'flux', format = '1D', array = std_flux )
        col3 = fits.Column(name = 'ivar', format = '1D', array = std_ivar)

        #mask used for telluric correction, bad pixels with mask = 0
        mask = np.full(len(std_flux), 1, dtype = int)
        #the std spec in both DRP and pypeit seems to have some problems here; mask it out for g19b2b, should check for other stds
        # mask[(self.std['wave'] >=6310) & self.std['wave'] <= 6380] = 0 
        col4 = fits.Column(name = 'mask', format = '1K', array = mask)
        hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
        # hdu.header['PYP_SPEC'] = 'keck_lris_red'
        hdulist = fits.HDUList([fits.PrimaryHDU(header = newhdr), hdu])
        hdulist[1].header['PYP_SPEC'] = 'keck_kcrm'

        #need to update the header for pypeit input
        keys_1 = { 'DMODCLS': 'OneSpec ', 'DMODVER': '1.0.2   ', 'FLUXED': True, 
        'CHECKSUM': 'CG5BCD59CD5ACD59', 'DATASUM': '60558086'
       }
        keys_2 = {'PYP_SPEC': 'keck_kcrm', 'PYPELINE': 'SlicerIFU',  'TARGET':newhdr['OBJECT'],
        'DISPNAME': self.std['invsens_hdr']['RGRATNAM'], 'decker': self.std['invsens_hdr']['IFUNAM'],  
        'binning': self.std['invsens_hdr']['BINNING'], 'FILENAME': '%s.fits'%frame,
        'NCOMP': 1,
        # 'AIRMASS':1.0 #modify the airmass to one because the std has been extinction corrected by DRP
        }

        for k in keys_1.keys():
            hdulist[1].header[k] = keys_1[k]
        for k in keys_2.keys():
            hdulist[0].header[k] = keys_2[k]
        std_spec1d_path = f'{self.output}/{frame}_standard_spec1d.fits'
        hdulist.writeto(f'{std_spec1d_path}', overwrite = True)

        #Generate the pypeit input file
        coord = SkyCoord(ra = self.std['invsens_hdr']['RA'], dec = self.std['invsens_hdr']['DEC'], unit = (u.hourangle, u.deg))
        # lines = ['[telluric]',
        #          f'\ttelgridfile = {telgridfile}',
        #          '\tobjmodel = star',
        #          f'\tstar_ra = {coord.ra.value:0.5f}',
        #          f'\tstar_dec = {coord.dec.value:0.5f}',
        #          '\tteltype = grid']

        # tellfile_out = (f'{self.output}/{frame}_standard.tell')
        # with open(tellfile_out, 'w') as file:
        #     file.writelines(lines)

        # Run the pypeit_tellfit command for telluric correction
        # command = f'pypeit_tellfit -t {tellfile_out} {std_spec1d_path}'
        # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        # # self.read_tellfit_output()
        # while process.poll() is None: # check whether process is still running
        #     msg = process.stdout.readline().strip() # read a line from the process output
        #     if msg:
        #         self.redirect_text.write(msg)
        # print('finished')
        # command = ['pypeit_tellfit', '-t', f'{tellfile_out}', f'{std_spec1d_path}']
        # # run_command(command, self.output_text)
        # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,)
        # process.wait()
        telluric_correct(infile_path = std_spec1d_path, star_ra = coord.ra.value, star_dec = coord.dec.value)

        # #the output files would be in the current directory. Need to move them to the output directory
        # files = ['telluric.par',  f'{frame}_standard_spec1d_tellcorr.fits', f'{frame}_standard_spec1d_tellmodel.fits']
        # for f in files:
        #     os.rename(f, f'{self.output}/{f}')

        #load the telluric model
        tellfile = fits.getdata(re.sub('.fits', '_tellcorr.fits', std_spec1d_path))
        tellmodel = tellfile['telluric'] #best-fit telluric model at the airmass of the standard star
        tellmodel = tellmodel**(1. / newhdr['AIRMASS']) #correct the telluric spectrum to be the one at airmass of 1
        self.std['tellmodel'] = tellmodel
        self.insert_text(f"[INFO] Telluric correction on {frame} DONE")


    def read_tellfit_output(self):
        if self.tellfit_process.poll() is None:  # Process is still running
            output_line = self.tellfit_process.stdout.readline()
            if output_line:
                self.redirect_text.write(output_line)
        self.root.after(100, self.read_tellfit_output)  # Check for new output every 100 milliseconds


    def on_click(self, event):
        """
        Mouse button used for general plotting analysis
        """
        # Right click to resume the home window
        if event.button == 3:  
            self.toolbar.home()

    def on_type(self, event):

        """
        Key interaction used for std analysis
        """
        # exclude/include the regions used for standard star analysis
        if event.key == 'e' or event.key == 'i':

            #select the starting point
            if self.std['region_start'] is None:
                self.std['region_start'] = event.xdata

            #select the ending point
            else:
                self.std['region_end'] = event.xdata   

                #flag the region to True to be masked out from fitting
                if event.key == 'e':
                    self.std['flag'][(self.std['wave']>= self.std['region_start']) & (self.std['wave']<= self.std['region_end'])] = 0
                
                #flag the region to False to be included in the fitting
                else:
                    self.std['flag'][(self.std['wave']>= self.std['region_start']) & (self.std['wave']<= self.std['region_end'])] = 1
                    
                self.std['region_start'] = None #reset the starting point for the next input
                self.plot_std(restore_limit = True) #update the std plot

        # add single one continuum data point closest to the mouse location; should be used in the regions with dense sky features
        if event.key == 'a':
            #find the index of the point closest to the mouse location
            # idx = np.argmin((self.std['wave'] - event.xdata)**2 + (self.std['counts'] * self.std['invsens_model_drp'] - event.ydata)**2) 
            idx = np.argmin(np.abs(self.std['wave'] - event.xdata))
            self.std['flag'][idx] = 2
            # print(self.std['use_ind'])
            self.plot_std(restore_limit = True)
        
        # delete single one continuum data point closest to the mouse location; should be used in the regions with dense sky features
        if event.key == 'd':
            if np.sum(self.std['flag'] ==2) > 0:
                # idx = np.argmin((self.std['wave'][self.std['use_ind']] - event.xdata)**2 + ((self.std['counts'] * self.std['invsens_model_drp'])[self.std['use_ind']] - event.ydata)**2) 
                ind_point = np.where(self.std['flag'] ==2 )[0]
                idx = np.argmin(np.abs(self.std['wave'][ind_point] - event.xdata))
                self.std['flag'][ind_point[idx]] = 0
                # self.std['use_ind'] = np.delete(self.std['use_ind'], idx)
                self.plot_std(restore_limit = True)

        if event.key == 'f':
            use = self.std['flag'] > 0
            # use[self.std['use_ind']] = True
            # print(len(use), self.std['wave'])
            self.std['invsens_model'] = self.fit_bspline(self.std['wave'], self.std['invsens_data'], self.std['bspline_bkpt'], self.std['bspline_polyorder'], use)
            
            self.plot_std(restore_limit = True)

        if event.key == 't':
            self.run_telluric()
            self.plot_std()



            # self.std['invsens_model'] = bspline()

    def fit_bspline(self, x_full, y_full, bkpt, k, use):
        """
        A simple BSpline-fit wrapper for flux calibration

        Args:
            x (1D arr): x data
            y (1D arr): y data
            bkpt (int): breakpoints (end point of each interval)
            use (boolen arry, same shape as x and y): True for pixels used for fitting
        """
        x = x_full[use]
        y = y_full[use]
        # x_idx = np.arange(x.size)
        # n_knots = int(x.size // bkpt)        
        # knots_idx = np.linspace(np.min(x[use])+ 0.1, np.max(x[use])-0.1, n_knots)
        # print(n_knots, knots)
        # print(x[use][0], x[use][-1])
        knots_idx = np.arange(bkpt, x.size, bkpt)
        # knots = x[knots_idx]

        # #insert the first and the last index to the knots
        # if x[knots_idx][-1] -1 < x[-1]:
        #     knots = np.append(knots, x[-1] - 0.5)
        # # print
        # knots = np.insert(knots, 0, x[0]+0.5)
        bspline = splrep(x, y, k = k, task=-1, t=x[knots_idx])

        return splev(x_full, bspline)


    # def set_zap_skyseg(self):
    #     if not self.use_multi_skyseg.get():



    
class RedirectText:
    """
    A class redirecting the system output to the text widget in the GUI
    """
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)

    def flush(self):
        pass


def update_text_box(text_box: tk.Text, proc: subprocess.Popen):
    """Create a closure to capture proc and text_box which will be called to update the text box"""

    # closure to capture proc and text_box
    def _update_text_box():
        for line in proc.stdout:
            text_box.insert(tk.END, line.decode())
            text_box.yview(tk.END)

    return _update_text_box


def run_command(command: list[str], text_box: tk.Text):
    """Run a command and redirect output to a text box in real-time"""
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    # Start a separate thread to update the text box with output
    thread = threading.Thread(target=update_text_box(text_box, proc))
    thread.start()
    proc.wait()


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1200x800")
    root.resizable(True, True)
    app = KCWIViewerApp(root)
    root.bind("<FocusIn>", app.focus_in)
    root.mainloop()