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

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.xlim(10,1000)
plt.ylim(1e-3,10)

plt.xlabel(r'$dN_{ch}/d\eta$')
plt.ylabel(r'$dN_{\gamma}/dy(p_T^\gamma>1 GeV)$')

colour_list=['blue','green','black','red','purple','orange','cyan']
#symbol_list=['-','--',':','-.','-','--',':']
symbol_list=['D','^','x','o','.','d','D']

filename="../../from_phenix/Ngamma_vs_Nhadron_phenix.dat"
dNch, dNch_err, dNg, dNg_stat, dNg_sys = np.loadtxt(filename).T
plt.errorbar(dNch,dNg,yerr=np.sqrt(dNg_stat**2+dNg_sys**2),fmt="r^",capsize=4, label="PHENIX") #, barsabove=True)

#filename="data/star_photons_hadron_mult.dat"
#dNch, dNg_pT1, dNg_stat_pT1, dNg_sys_pT1, dNg_pT15, dNg_stat_pT15, dNg_sys_pT15 = np.loadtxt(filename).T
#plt.errorbar(dNch,dNg_pT1,yerr=np.sqrt(dNg_stat_pT1**2+dNg_sys_pT1**2),fmt="gD",capsize=4, label="STAR") #, barsabove=True)

filename="./calcs/photons/ipg_kompost_cent0060_20pctbins/dNh_vs_dNg_AuAu200_source=prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem1_pTcut=1.0.dat"
dNch, dNg = np.transpose(np.loadtxt(filename))
plt.plot(dNch, dNg,"-", color='black', label="Prompt+pre.eq+thermal [slope=1.28]")

#filename="calcs/photons/ipg_kompost/dNh_vs_dNg_source=kompost_plus_thermal_reg1_pTcut=1.0.dat"
#dNch, dNg = np.transpose(np.loadtxt(filename))
#plt.plot(dNch, dNg,"o", color='green', label="Thermal")

filename="./calcs/photons/ipg_kompost_cent0060_20pctbins/dNh_vs_dNg_AuAu200_source=prompt_pTcut=1.0.dat"
dNch, dNg = np.transpose(np.loadtxt(filename))
plt.plot(dNch, dNg,"--", color='blue', label="Prompt [slope=1.11]")

plt.legend(loc='upper left',fontsize=10)
plt.tight_layout()
plt.savefig("vs_data_dNg_vs_dNch.pdf")
plt.show()

