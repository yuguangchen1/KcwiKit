#!/usr/bin/env python3
import sys
import glob
import os
from astropy.io import fits
import numpy as np

# Import the deepCR module. Adjust the import as needed if the API is different.
try:
    from deepCR import CosmicRayDetector
except ImportError:
    sys.exit("Error: deepCR module not found. Please install deepCR before running this script.")

def process_file(filename, detector):
    """Read a FITS file, run deepCR detection, and save the cosmic ray mask."""
    print(f"Processing file: {filename}")
    
    # Open the FITS file and extract the image data (assumed to be in the primary HDU)
    with fits.open(filename) as hdul:
        image_data = hdul[0].data
    
    # Run cosmic ray detection using deepCR. This returns a boolean mask.
    mask = detector.detect(image_data)
    
    # Optionally, you can convert the boolean mask to an integer type (0 and 1) for saving.
    mask_int = mask.astype(np.int16)
    
    # Create an output filename by appending '_crmask' before the .fits extension.
    base, ext = os.path.splitext(filename)
    output_filename = f"{base}_crmask{ext}"
    
    # Write the cosmic ray mask to a new FITS file.
    hdu = fits.PrimaryHDU(mask_int)
    hdu.writeto(output_filename, overwrite=True)
    print(f"Saved cosmic ray mask to: {output_filename}")

def main():
    # Check if at least one file pattern was provided as a command-line argument.
    if len(sys.argv) < 2:
        print("Usage: python deepcr_script.py <file_or_pattern>")
        sys.exit(1)

    # Expand file patterns (e.g. 'kr*.fits') to a list of filenames.
    file_list = []
    for arg in sys.argv[1:]:
        file_list.extend(glob.glob(arg))
    
    if not file_list:
        sys.exit("No files found matching the given pattern(s).")
    
    # Create a deepCR detector instance once to reuse for all images.
    detector = CosmicRayDetector()

    # Process each file.
    for filename in file_list:
        process_file(filename, detector)

if __name__ == '__main__':
    main()