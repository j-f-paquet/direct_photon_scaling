import numpy as np
import scipy.interpolate
import scipy.integrate
import os.path
import glob
from scipy.stats import linregress

####################################################################
################ Information about the calculations ################ 
####################################################################

system_choice='AuAu200'

calc_dict={
    "AuAu200":{
        #  0-5%, 5-10%, 10-15%, 15-20%, 20-30%, 30-40%, 40-50%, 50-60%
        'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55'},
        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55','C60-70':'65','C70-80':'75','C80-90':'85','C90-100':'95'},
        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35'},
        #'centralities':{'C0-20': '10', 'C20-40':'30', 'C40-60':'50'},
        'hadron_file':"../../../../../calculations_from_multimessenger_paper/hadrons/ipg_kompost/AuAu200.dat",
        'photon_file_dir':"../../../../../calculations_from_multimessenger_paper/photons/",
        'photon_file_name':"average_rms.dat",
        },
#    "PbPb2760":{
#        'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55'},
#        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55','C60-70':'65','C70-80':'75','C80-90':'85','C90-100':'95'},
#        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35'},
#        'hadron_file':"../../hadrons/ipg_kompost/PbPb2760.dat",
#        'photon_file_dir':"../../../../../calcs/kompost_plus_hydro/",
#        'photon_file_name':"average_rms.dat",
#        },
#    "PbPb5020":{
#        'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55'},
#        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35'},
#        #'centralities':{'C0-5': '2.5', 'C5-10':'7.5', 'C10-20':'15','C20-30':'25','C30-40':'35','C40-50':'45','C50-60':'55','C60-70':'65','C70-80':'75','C80-90':'85','C90-100':'95'},
#        'hadron_file':"../../hadrons/ipg_kompost/PbPb5020.dat",
#        'photon_file_dir':"../../../../../calcs/kompost_plus_hydro/",
#        'photon_file_name':"average_rms.dat",
#        },
}

# Photon sources
#photon_sources=['thermal_Tfr105','kompost_and_thermal_Tfr105','kompost_with_suppression_and_thermal_Tfr105','prompt_kompost_and_thermal_Tfr105','prompt_kompost_with_suppression_and_thermal_Tfr105','prompt']
#photon_sources=['kompost_plus_thermal_plus_prompt','prompt','kompost','thermal_reg1','thermal_reg05','kompost_plus_thermal_reg1_plus_prompt']
photon_sources=['prompt', "prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem1"]
#photon_sources=['prompt', "hydro_with_visc_T105_reg1", "kompost_plus_hydro_thermal_with_visc_T105_reg1", "prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1","prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem1","prompt_plus_kompost_plus_hydro_thermal_with_visc_T105_reg1_suppr_tauChem15"]

# p_T cut for the photon momentum
photon_low_pT_cut_list=np.arange(1.0,10.0,1.0)
photon_high_pT_cut=14.0

#######################################################
################ Load the calculations ################ 
#######################################################

curr_calc_dict=calc_dict[system_choice]

hadron_dNdy_file=curr_calc_dict['hadron_file']

centrality_dict=curr_calc_dict['centralities']

################# Load the hadron multiplicities #################

# Load hadron dN/dy
hadron_dNdy_list=np.loadtxt(hadron_dNdy_file)[:,(0,1)]

# Make dictionary out of hadron dN/dy
hadron_dNdy_dict={cent: hadron_dNdy for (cent, hadron_dNdy) in hadron_dNdy_list if float(cent) in [ float(elem) for elem in list(centrality_dict.values()) ]}

################# Loading photon results and calculating their dN/dy #################

# Location of photon results, with @cent@ standing for centrality in 0010, 1020, ... format
# and @source@ standing for the photon source
photon_file_location_pre=os.path.join(curr_calc_dict['photon_file_dir'],"@source@",system_choice,"@cent@/",curr_calc_dict['photon_file_name'])

def get_dNch_and_dNgamma(source,photon_low_pT_cut):

    dNh=[]
    dNg=[]

    # Loop over centralities
    for cent_long, cent_short in centrality_dict.items():

        dNh.append(hadron_dNdy_dict[float(cent_short)])

        # Get file name
        tmp_file_name=photon_file_location_pre
        tmp_file_name=tmp_file_name.replace('@source@',source)
        tmp_file_name=tmp_file_name.replace('@cent@',cent_long)

        photon_file_list=glob.glob(tmp_file_name)
        if (len(photon_file_list)!=1):
            print("Can't find photon calculation in ",photon_file_location_pre)
            exit(1)
        tmp_file_name=photon_file_list[0]

        # Load photon results
        tmp_photon_res=np.loadtxt(tmp_file_name)

        # Store only p_T and dN/dp_T (column 0 and 1)
        # Skip the first two p_T point
        tmp_photon_pT, tmp_photon_Ed3Nd3p=np.transpose(tmp_photon_res[2:,(0,1)])
        # What I must integrate is 2*pi*p_T*(1/2pip_T dN/dp_T)
        # Compute dN/dp_T first
        tmp_dNdpT=2*np.pi*np.multiply(tmp_photon_pT,tmp_photon_Ed3Nd3p)
        # More precise integration with linear interpolation of the log of dN/dpT??
        tmp_interpolator=scipy.interpolate.interp1d(tmp_photon_pT,np.log(tmp_dNdpT),kind='linear',bounds_error=True,assume_sorted=True)
        tmp_integral=scipy.integrate.quad(lambda x : np.exp(tmp_interpolator(x)),photon_low_pT_cut,photon_high_pT_cut,limit=300)
        # Drop the error term 
        tmp_integral=tmp_integral[0]

        dNg.append(tmp_integral)

        # Somewhat less precise integration, with quadratic interpolation of dN/dpT?
        #tmp_interpolator2=scipy.interpolate.interp1d(tmp_photon_pT,tmp_dNdpT,kind='quadratic',bounds_error=True,assume_sorted=True)
        #tmp_integral=scipy.integrate.quad(lambda x : tmp_interpolator2(x),photon_low_pT_cut,photon_high_pT_cut)
        #tmp_integral=tmp_integral[0]

    return dNh, dNg


def get_slope_and_intercept(photon_low_pT_cut):

    tmp_dict={}

    # Loop over photon sources
    for source in photon_sources:

        ####################### For each source, find the slope #######################

        dNh, dNg = get_dNch_and_dNgamma(source,photon_low_pT_cut)

        # Save dN_g vs dN_h
        np.savetxt("dNh_vs_dNg_"+system_choice+"_source="+source+"_pTcut="+str(photon_low_pT_cut)+".dat",np.transpose((dNh,dNg)))

        # Find the slope
        slope, intercept, r_value, p_value, stderr=linregress(np.log(np.array(dNh)),np.log(np.array(dNg)))

        tmp_dict[source]={'slope':slope, 'intercept':intercept}

    return tmp_dict


###############################################################
################ Calculate and save the slopes ################ 
###############################################################

## Name of the file where results are saved
##result_file_name='summary_'+str(photon_low_pT_cut)+'_'+str(photon_high_pT_cut)+'_T105_intLog.dat'
#result_file_name='beta_'+system_choice+'_vs_pT_cut.dat'
#
## Open file where results are saved
#result_file = open(result_file_name, 'w')
#
## Write header
#header="#p_T-cut "
#for source in photon_sources:
#    header+="beta_"+source+"_photon "
#header+="\n"
#result_file.write(header)
#
## Loop over p_T-min cuts
#for photon_low_pT_cut in photon_low_pT_cut_list:
#
#    #tmp_dict={source:{'slope':slope, 'intercept':intercept}}
#    tmp_dict=get_slope_and_intercept(photon_low_pT_cut)
#
#    line_to_write=str(photon_low_pT_cut)+" "
#
#    # Loop over photon sources
#    for source, tmp_dict in tmp_dict.items():
#
#        slope=tmp_dict['slope']
#
#        line_to_write+=str(slope)+" "
#
#    line_to_write+="\n"
#
#    result_file.write(line_to_write)
#
#            
#result_file.close()


# Write header
#header="#p_T-cut "
#for source in photon_sources:
#    header+="beta_"+source+"_photon "
#header+="\n"
#result_file.write(header)

#for source, tmp_dict in tmp_dict.items():

# Loop over photon sources
for source in photon_sources:

    # Name of the file where results are saved
    #result_file_name='summary_'+str(photon_low_pT_cut)+'_'+str(photon_high_pT_cut)+'_T105_intLog.dat'
    result_file_name='beta_vs_pT_cut_'+system_choice+'_source='+source+'.dat'

    # Open file where results are saved
    result_file = open(result_file_name, 'w')

    # Loop over p_T-min cuts
    for photon_low_pT_cut in photon_low_pT_cut_list:

        #tmp_dict={source:{'slope':slope, 'intercept':intercept}}
        tmp_dict=get_slope_and_intercept(photon_low_pT_cut)

        line_to_write=str(photon_low_pT_cut)+" "

        slope=tmp_dict[source]['slope']
        intercept=tmp_dict[source]['intercept']

        line_to_write+=str(slope)+" "+str(intercept)+"\n"

        result_file.write(line_to_write)
            
    result_file.close()
