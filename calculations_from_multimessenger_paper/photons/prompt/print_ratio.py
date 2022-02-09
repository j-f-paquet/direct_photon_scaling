import numpy as np

cent_list=["C0-5","C5-10","C10-20","C20-30","C30-40","C40-50","C50-60"]

for cent in cent_list:
   
    dN_phenix=np.transpose(np.loadtxt("AuAu200_phenix_Ncoll/"+cent+"/average_rms.dat"))[1]
    dN_calc=np.transpose(np.loadtxt("AuAu200/"+cent+"/average_rms.dat"))[1]

    print(cent, dN_phenix/dN_calc)
