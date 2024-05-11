import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from glob import glob
from astropy.io import fits
import os 
import kcwi_tools 
import pyregion



"""
Author: Zhuyun Zhuang, modified from chatGPT (free version!)
Date: 05/2024
Usage: python kcwi_zap_viewer.py

Description: Load data from a directory and allow user to interactively
perform sky subtraction of KCWI/KCRM data with ZAP, and flux calibration. 
"""

#only set it for test purporse.
initial_dir = '/scr/zzhuang/keck_obs/kcwi/2023sep23/red/redux/'

#pick this region to generate the white-lighted image because sky lines are much stronger elsewhere. 
# TODO: Can also make it as an input or variable parameter in the GUI
wlimg_wave_range = [6380, 7200] 
class KCWIViewerApp:
    def __init__(self, root):
        """
        Initialize the GUI window
        """
        self.root = root
        self.root.title("KCWI Data Viewer")

        # Input directory
        self.input_dir_label = tk.Label(root, text="Input Directory:")
        self.input_dir_label.grid(row=0, column=0)
        self.input_dir_entry = tk.Entry(root)
        self.input_dir_entry.grid(row=0, column=1)
        self.browse_button = tk.Button(root, text="Browse", command=self.browse_input_directory)
        self.browse_button.grid(row=0, column=2)

        #Output directory
        self.output_dir_label = tk.Label(root, text="Output Directory:")
        self.output_dir_label.grid(row=0, column=3)
        self.output_dir_entry = tk.Entry(root)
        self.output_dir_entry.grid(row=0, column=4)
        self.browse_button = tk.Button(root, text="Browse", command=self.browse_output_directory)
        self.browse_button.grid(row=0, column=5)

        # Science index input
        self.index_label = tk.Label(root, text="Science Frame No.:")
        self.index_label.grid(row=1, column=0)
        self.index_entry = tk.Entry(root)
        self.index_entry.grid(row=1, column=1)
        self.index_entry.bind("<Return>", self.update_index)

        # Buttons for increasing and decreasing index for science frame
        self.increase_button = tk.Button(root, text="Previous", command=self.decrease_index)
        self.increase_button.grid(row=1, column=2)
        self.decrease_button = tk.Button(root, text="Next", command=self.increase_index)
        self.decrease_button.grid(row=1, column=3)

        #Check box to select whether an off-field sky is used
        self.use_index2_var = tk.BooleanVar()
        self.use_index2_checkbox = tk.Checkbutton(root, text="Use Off-field Sky Frame No.", variable=self.use_index2_var, command=self.toggle_index2_entry)
        self.use_index2_checkbox.grid(row=1, column=4)

        # sky index input
        self.index2_entry = tk.Entry(root)
        self.index2_entry.grid(row=1, column=5)
        self.index2_entry.bind("<Return>", self.update_index2)

        #Input the indices for the regions used for spectrum extraction, and specify the redshift
        # self.redshift = tk.

        # Naming structure selection - used for loading the data cube
        self.structure_label = tk.Label(root, text="Data Product type:")
        self.structure_label.grid(row=2, column=0)
        self.structure_var = tk.StringVar()
        self.structure_options = ["icubed", "icubes", "icube"]
        self.structure_var.set('icubed')
        self.structure_menu = tk.OptionMenu(root, self.structure_var, *self.structure_options)
        self.structure_menu.grid(row=2, column=1)

        #  Load raw data (DRP-reduced cube) button
        self.load_button = tk.Button(root, text="Load Raw Cube", command=lambda: self.load_data('raw'))
        self.load_button.grid(row=2, column=2)

        #  Save cropped data (good wavelength region) button
        self.load_button = tk.Button(root, text="Save Cropped Cube", command=self.save_cropped_data)
        self.load_button.grid(row=2, column=3)

        # Load cropped data button
        self.load_crop_button = tk.Button(root, text="Load Cropped Cube", command=lambda: self.load_data('cropped'))
        self.load_crop_button.grid(row=2, column=4)

        # Input the mask frame No. (can be different from the science one since we usally take the at least three frames)
        # This function is only used for convenience, so that for each pair we only need to create one 
        self.mask_index_label = tk.Label(root, text="Mask Frame No.:")
        self.mask_index_label.grid(row=3, column=0)
        self.mask_entry = tk.Entry(root)
        self.mask_entry.grid(row=3, column=1)
        # self.mask_entry.bind("<Return>", self.update_index)
        self.update_mask_button = tk.Button(root, text = 'Update ZAP mask', command=self.update_zap_mask)
        self.update_mask_button.grid(row =3, column = 2)

        #Setup the continuum 


        # Plotting area
        self.figure = Figure(figsize = (12, 5))
        self.ax = self.figure.add_subplot(111)
        # self.figure.tight_layout()
        # Global font settings
        plt.rc('font', family='serif', size=12)

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=4, column=0, columnspan=6)

        # Add navigation toolbar
        self.toolbarFrame = tk.Frame(master=root)
        self.toolbarFrame.grid(row=5,column=0, columnspan = 6)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)

        # Initialize indices
        self.index = -1 #science frame number
        self.index2 = -1 #sky frame number; if specified
        self.index2_entry.config(state=tk.DISABLED)
        self.update_index_entries()

        # Output text widget
        self.output_text = tk.Text(root, height=5, width=100, font=("Arial", 12))
        self.output_text.grid(row=6, column=0, columnspan=6)

        root.columnconfigure(1, minsize=100, weight=0)


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
        self.index_entry.delete(0, tk.END)
        self.index_entry.insert(tk.END, str(self.index))
        self.index2_entry.delete(0, tk.END)
        self.index2_entry.insert(tk.END, str(self.index2))
        if self.index > 0:
            self.output_text.insert(tk.END, f"[INFO] Set the science frame: {self.prefix}_{self.index:05d}\n")
            if self.index2 > 0:
                self.output_text.insert(tk.END, f"[INFO] Set the sky frame for science frame {self.prefix}_{self.index:05d}: {self.prefix}_{self.index2:05d} \n")


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

    def browse_input_directory(self):
        """Select the input directory"""

        #specify the input directory
        directory = filedialog.askdirectory(initialdir=initial_dir)
        self.input_dir_entry.delete(0, tk.END)
        self.input_dir_entry.insert(tk.END, directory)
        self.base = self.input_dir_entry.get()
        self.output_text.insert(tk.END, f"[INFO] Loading the data from {self.base}\n")

        #find the common string name for a given date
        self.prefix = os.path.basename(glob(f'{directory}/kr*.fits')[0])[:8]

    def browse_output_directory(self):
        """Select the output directory"""
        directory = filedialog.askdirectory()
        self.output_dir_entry.delete(0, tk.END)
        self.output_dir_entry.insert(tk.END, directory)
        self.output = self.output_dir_entry.get()
        # print(f'[INFO] Set the output directory to {self.output}')
        self.output_text.insert(tk.END, f"[INFO] Set the output directory to {self.output}\n")

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
            self.output_text.insert(tk.END, f'[ERROR] datatype not recognized! Need to be "raw" or "cropped" \n')

        # load the science frame
        if self.index > 0:
            self.scihdu = fits.open(f'{base}/{self.prefix}_{self.index:05d}_{ctype}.fits')
            self.objname = self.scihdu[0].header['OBJECT']
            if datatype == 'raw':
                self.scihdu = self.crop_cube(self.scihdu) #crop the data cube to good wavelength region
                self.output_text.insert(tk.END, f"[INFO] Loading the DRP-reduced science frame {self.prefix}_{self.index:05d}\n")
            elif datatype == 'cropped':
                self.output_text.insert(tk.END, f"[INFO] Loading the cropped DRP-reduced science frame {self.prefix}_{self.index:05d}\n")

            self.scihdr = self.scihdu[0].header
            self.obswave = (np.arange(self.scihdr['NAXIS3']) + 1 - self.scihdr['CRPIX3']) * self.scihdr['CD3_3'] + self.scihdr['CRVAL3']

        else:
            self.output_text.insert(tk.END, f"[ERROR] Wrong science frame! Need to set it to a positive integer. Check the KCWI log for the frame number!\n")

        # self.z = 0
        if self.index2 > 0:
            self.skyhdu = fits.open(f'{base}/{self.prefix}_{self.index2:05d}_{ctype}.fits')
            if datatype == 'raw':
                self.skyhdu = self.crop_cube(self.skyhdu) #crop the data cube to good wavelength region
                self.output_text.insert(tk.END, f"[INFO] Loading the DRP-reduced sky frame {self.prefix}_{self.index2:05d}\n")
            else:
                self.output_text.insert(tk.END, f"[INFO] Loading the cropped DRP-reduced sky frame {self.prefix}_{self.index:05d} for {self.prefix}_{self.index:05d}\n")
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

                self.output_text.insert(tk.END, f"[INFO] Saving the cropped datacube for {self.prefix}_{index:05d}\n")

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

        #get the mask frame number
        if not self.mask_entry.get():
            self.output_text.insert(tk.END, f"[ERROR] Need to specify the frame used for generating the object mask!\n")
        
        mindex = int(self.mask_entry.get())

        mhdu = fits.open(f'{self.output}/{self.prefix}_{self.index:05d}_wlimg.fits')
        edgemask = mhdu[1].data > 1

        region_path = f'{self.output}/{self.prefix}_{mindex:05d}.reg'
        if not os.path.exists(region_path):
            self.output_text.insert(tk.END, f"[ERROR] Region file {self.prefix}_{mindex:05d}.reg not exists!\n")
        else:
            self.output_text.insert(tk.END, f"[INFO] Reading region file {self.prefix}_{mindex:05d}.reg for {self.prefix}_{self.index:05d}\n")
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
            self.output_text.insert(tk.END, f"[INFO] Saving ZAP mask to {self.prefix}_{self.index:05d}_zap_mask.fits\n")
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


    def plot_spectrum(self, hdu, z = 0, xrange = [12, 22], yrange = [43, 53], hdu_sky = None):
        """
        plot the spectrum of a given region
        
        Args:
            hdu (fits.HDUList) of a given data cube
            z (float): the redshift of the object
            xrange (list): the start and ending of the x-axis from which the spec is extracted
            yrange (list): the start and ending of the x-axis from which the spec is extracted
        """
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
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1200x800")
    root.resizable(True, True)
    app = KCWIViewerApp(root)
    root.mainloop()