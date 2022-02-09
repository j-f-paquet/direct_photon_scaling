import numpy as np
import matplotlib as mpl
mpl.use('Agg')
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

fig = plt.figure()
plt.xscale('linear')
plt.yscale('linear')
plt.xlim(0,9)
plt.ylim(0.8,1.7)
#plt.yticks([])
#plt.axes().yaxis.set_minor_formatter(NullFormatter())
#plt.axes().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
#plt.yticks(np.arange(.1, .450, .05))
plt.xlabel(r'$p_T^{min} (GeV)$')
plt.ylabel(r'$\alpha$')


#pTcuts, alpha, d_alpha = np.transpose(np.loadtxt("data/phenix.dat")[:,(0,1,2)])
#plt.errorbar(pTcuts, alpha, yerr=d_alpha, fmt="D",  color='darkgreen', label="PHENIX")

pTmin, slope, intercept = np.transpose(np.loadtxt(os.path.join("./calcs/photons/ipg_kompost_cent0060_20pctbins/beta_vs_pT_cut_AuAu200_phenix_Ncoll_source=prompt.dat")))
plt.plot(pTmin, slope, ":", color="grey", label="Prompt")

pTmin, slope, intercept = np.transpose(np.loadtxt(os.path.join("./calcs/photons/ipg_kompost_cent0060_varied_cents/beta_vs_pT_cut_AuAu200_phenix_Ncoll_source=prompt.dat")))
mask=pTmin>3.5
plt.plot(pTmin[mask], slope[mask], ":", color="grey", label="")


pTmin, slope, intercept = np.transpose(np.loadtxt(os.path.join("./calcs/photons/ipg_kompost_cent0060_20pctbins/beta_vs_pT_cut_AuAu200_phenix_Ncoll_source=prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem1.dat")))
plt.plot(pTmin, slope, "-", color="Red", label="Prompt+thermal")

pTmin, slope, intercept = np.transpose(np.loadtxt(os.path.join("./calcs/photons/ipg_kompost_cent0060_varied_cents/beta_vs_pT_cut_AuAu200_phenix_Ncoll_source=prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem1.dat")))
mask=pTmin>3.5
plt.plot(pTmin[mask], slope[mask], "-", color="Red", label="")

#plt.plot(pT_mins, tmp_data[:,index+1],symbol_list[n],color=colour_list[n], label=label)
#plt.plot(x_range, np.exp(regress_dict[source]['slope']*np.log(x_range)+regress_dict[source]['intercept']),"--",color=colour_list[n])
plt.text(0.02, 0.95, r"$N_{coll}$ and $dN_{ch}/d\eta$""\n""from PHENIX", horizontalalignment='left', verticalalignment='top', transform=fig.axes[0].transAxes, fontsize=12)

plt.arrow(0.3, 0.25, -0.26, 0, linewidth=2, transform=fig.axes[0].transAxes)
plt.text(0.02, 0.22, r"$p_T^{max}=5$ GeV", horizontalalignment='left', verticalalignment='top', transform=fig.axes[0].transAxes, fontsize=12)
plt.text(0.02, 0.15, "Centralities:\n0-20%, 20-40%,\n40-60%", horizontalalignment='left', verticalalignment='top', transform=fig.axes[0].transAxes, fontsize=11)

plt.arrow(0.4, 0.25, 0.55, 0, linewidth=2, transform=fig.axes[0].transAxes)
plt.text(0.41, 0.22, r"$p_T^{max}=14$ GeV", horizontalalignment='left', verticalalignment='top', transform=fig.axes[0].transAxes, fontsize=12)
plt.text(0.41, 0.15, "Centralities:\n0-5%, 5-10%, 10-20%, 20-30%,\n30-40%, 40-50%, 50-60%", horizontalalignment='left', verticalalignment='top', transform=fig.axes[0].transAxes, fontsize=11)

plt.legend(loc='upper right',fontsize=16)
plt.tight_layout()
plt.savefig("alpha_vs_pTmin_phenix_scaling.pdf")
plt.show()

