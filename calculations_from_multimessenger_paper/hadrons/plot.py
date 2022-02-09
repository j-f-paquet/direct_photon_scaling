import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib import colors, ticker, cm
#import scipy.interpolate
#from matplotlib.mlab import bivariate_normal
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mtick
from scipy.stats import linregress
import pandas as pd
import os.path
import glob
import re


font = {'family' : 'URW Gothic',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)

#for init_cond, source_list in system_dict.items():
for init_cond in ['ipg','ipg_kompost']:

    plt.figure()
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0,60)
    plt.ylim(30,1000)

    plt.xlabel(r'Centrality')
    plt.ylabel(r'$dN_{ch}/d\eta$')

    cent_data, dN_data, dN_err_data = np.transpose(np.loadtxt("./data/dNchdeta_phenix.dat"))
    cent_calc, dN_calc, dN_stat_calc = np.transpose(np.loadtxt(os.path.join(init_cond,"AuAu200.dat"))[:,(0,1,2)])

    plt.errorbar(cent_data, dN_data, yerr=dN_err_data, fmt="D",  color='r', label="PHENIX")
    #plt.plot(cent_data, dN_data, "D", color='r', label="PHENIX")
    plt.plot(cent_calc, dN_calc, "-", color='b', label="Calculation ["+init_cond+"]")
    plt.fill_between(cent_calc, dN_calc-dN_stat_calc, dN_calc+dN_stat_calc) #, where=y2 >= y1, facecolor='green', interpolate=True)

    plt.legend(loc='upper right',fontsize=14)
    plt.tight_layout()
    plt.savefig("log_dNch_vs_cent_"+init_cond+".pdf")
    plt.show()

