# CR Final Rejection
# data0 looks like (#, wave, y,x) = (15, 2596, 100, 100)
from astropy.io import fits
import multiprocessing
from tqdm import tqdm
import numpy as np
import sys

# check args
narg=len(sys.argv)

# should be six (including routine name)
if narg != 6:
    print("Usage: python kcwi_crr.py <crrfn> <crrvfn> <crrmfn> <crrefn> [threshold {100}]")
    exit()

# read arg values
crrfn=sys.argv[1]
crrvfn=sys.argv[2]
crrmfn=sys.argv[3]
crrefn=sys.argv[4]

if sys.argv[5]:
    threshold = sys.argv[5]
else:
    threshold = 100

icube = fits.open(crrfn)[0].data
vcube = fits.open(crrvfn)[0].data
mcube = fits.open(crrmfn)[0].data
ecube = fits.open(crrefn)[0].data

waverange = range(icube.shape[1]) #range(310,2310) #range(data.shape[1]) #range(915,930)

def onepix(coords):
    x,y = coords

    if np.sum(icube[:,:,x,y]) == 0:
        return
    for wave in waverange:
        noise = np.sqrt(np.median(vcube[:,wave, x,y]))

        for exp in range(icube.shape[0]):
            if icube[exp,wave,y,x] > threshold*noise:
                return (exp,wave,y,x, 1)


if __name__ == '__main__':

    pool = multiprocessing.Pool() # Create a multiprocessing Pool
    coords = [(x, y) for x in range(icube.shape[3]) for y in range(icube.shape[2])]
    results = list(tqdm(pool.imap(onepix, coords), total=len(coords)))

    for out in results:
        if out is not None:
            icube[out[0],out[1],out[2],out[3]] = np.nan
            vcube[out[0],out[1],out[2],out[3]] = np.nan
            mcube[out[0],out[1],out[2],out[3]] = 4 # CR
            ecube[out[0],out[1],out[2],out[3]] = 0

    fits.PrimaryHDU(icube).writeto(crrfn, overwrite=True)
    fits.PrimaryHDU(vcube).writeto(crrvfn, overwrite=True)
    fits.PrimaryHDU(mcube).writeto(crrmfn, overwrite=True)
    fits.PrimaryHDU(ecube).writeto(crrefn, overwrite=True)
