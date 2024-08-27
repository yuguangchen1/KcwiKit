import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
# from scipy.signal import convolve

dir = '/Volumes/Data/Documents/Chris/KCRM_commissioning/VV340a/small_slicer'

pos = {
'cent': range(150,157)
}

msky = np.median([fits.open(f'{dir}/redux/kr230716_{i:05d}_int.fits')[0].data for i in range(161, 165)],axis=0)

msci = np.median([fits.open(f'{dir}/redux/kr230716_{i:05d}_int.fits')[0].data for i in pos['cent']],axis=0)




for i in pos['cent']:
    print(i)
    img = fits.open(f'{dir}/redux/kr230716_{i:05d}_int.fits')[0].data

    crs = img - msci
    crmsk = sigma_clip(crs, sigma = 3, maxiters=1, grow=2).mask.astype(float)

    crmsk_med = crmsk * msci
    crmsk_med[crmsk_med < 1e-4] = 0

    primary = fits.PrimaryHDU(crmsk)
    med = fits.ImageHDU(crmsk_med, name='MEDSCI')
    hdul = fits.HDUList([primary, med])

    hdul.writeto(f'{dir}/redux/kr230716_{i:05d}_crmsk.fits', overwrite=True)
