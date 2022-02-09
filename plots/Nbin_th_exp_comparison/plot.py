import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mtick
import os
import scipy
from scipy import interpolate

font = {'family' : 'URW Gothic',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)

plt.figure()
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel(r'$N_{bin} (PHENIX)$')
plt.ylabel(r'$N_{bin} (arXiv:2106.11216)$')

dummy, cent_mid_th, Ncoll_th, Ncoll_th_stat_uncert = np.transpose(np.genfromtxt(os.path.join("../../calculations_from_multimessenger_paper/Ncoll/AuAu200.dat")))
cent_mid_exp, Ncoll_exp, Ncoll_exp_stat_uncert = np.transpose(np.loadtxt(os.path.join("../../from_phenix/Ncoll_phenix.dat")))

#Ncoll_th_fct_linear=scipy.interpolate.interp1d(cent_mid_th, Ncoll_th, kind='linear', bounds_error=True)
#Ncoll_th_fct_quad=scipy.interpolate.interp1d(cent_mid_th, Ncoll_th, kind='quadratic', bounds_error=True)
Ncoll_th_fct_cubic=scipy.interpolate.interp1d(cent_mid_th, Ncoll_th, kind='cubic', bounds_error=True)

#print(cent_mid_th)
#print(cent_mid_exp)
#print(Ncoll_th_fct(cent_mid_exp)-cent_mid_exp)

#plt.plot(Ncoll_exp, Ncoll_th_fct_linear(cent_mid_exp), "D", color="red")
#plt.plot(Ncoll_exp, Ncoll_th_fct_quad(cent_mid_exp), "D", color="green")
plt.plot(Ncoll_exp, Ncoll_th_fct_cubic(cent_mid_exp), "D", color="blue")
plt.plot(Ncoll_exp, Ncoll_exp, ":", color="grey")

#plt.legend(loc='upper right',fontsize=16)
plt.tight_layout()
plt.savefig("Nbin_comparison.pdf")
plt.show()

