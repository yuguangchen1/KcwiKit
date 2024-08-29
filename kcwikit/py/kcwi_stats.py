import numpy as np
import pdb



def iter_polyfit(x,y,deg,max_iter=5,nsig=2.5):
    """

    """
    
    y1=y.copy()
    x1=x.copy()

    index=(np.isfinite(y1) & (y1!=0))
    if np.sum(index)==0:
        poly_fit=np.poly1d([0])
        return poly_fit
    y1=y1[index]
    x1=x1[index]

    for i in range(max_iter):
        if deg>=1:
            param=np.polyfit(x1,y1,deg)
            poly_fit=np.poly1d(param)
            y_fit=poly_fit(x1)

            residual=y1-y_fit

            rms=np.sqrt(np.mean(residual**2))
            index=(residual < nsig*rms)
            
            if np.sum(~index)==0:
                break
            else:
                x1=x1[index]
                y1=y1[index]

        else:
            med=np.median(y1)
            poly_fit=np.poly1d([np.mean(y1)])
            residual=y1-med

            rms=np.sqrt(np.mean(residual**2))
            index=(residual < nsig*rms)

            if np.sum(~index)==0:
                break
            else:
                x1=x1[index]
                y1=y1[index]
    
    return poly_fit


            





