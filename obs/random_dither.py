import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy import table


def parser_init():
    description = 'Generate random offset points for dithering.'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--xsig', dest='xsig',
            type=int, default=4, help='pixel sigma in x direction')
    parser.add_argument('--ysig', dest='ysig', 
            type=int, default=5, help='pixel sigma in y direction')
    parser.add_argument('--npoint', dest='npoint', 
            type=int, default=1, help='number of pointings')

    args = parser.parse_args()
    return args


def main(xsig=3, ysig=4, npoint=1):

    fig, ax = plt.subplots()
    ax.set_aspect(1)

    # FoV
    ax.plot([210, 300, 300, 210, 210], [288, 288, 234, 234, 288], color='blue', ls='--')
    ax.plot([233, 277, 277, 233, 233], [288, 288, 234, 234, 288], color='magenta', ls='--')
    ax.plot([243, 266, 266, 243, 243], [288, 288, 234, 234, 288], color='green', ls='--')

    # center
    ax.scatter(255, 261, color='black', marker='+')

    # random points
    xgs = np.random.normal(loc=255, scale=xsig, size=npoint)
    ygs = np.random.normal(loc=261, scale=ysig, size=npoint)

    for xg, yg in zip(xgs, ygs):
        ax.scatter(xg, yg, alpha=0.5)

    
    # print
    np.set_printoptions(threshold=np.inf)
    tab = table.Table([np.arange(npoint), xgs, ygs], names=['#', 'X', 'Y'])
    tab.pprint_all()



    plt.show()


    return


if __name__=='__main__':

    args = parser_init()

    main(**vars(args))
