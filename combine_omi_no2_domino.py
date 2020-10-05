import numpy as np
import fileinput
import h5py
import glob
import sys
import fileinput
import read_omi_no2_domino


year=2017




meses=np.arange(12,dtype=int)+1
for 	mes in meses:
#for 	mes in [11,12]:
#for 	mes in [1,2,3,4,5,6,7,8,9,10,11]:
	datapath='/home/STG_04/NO2_OMI/L2_V003/OMI-Aura_L2-OMNO2_%im%02i*.he5'% (year,mes)
	#datapath='%s-data/iasi_CO_LATMOS_ULB_metop%s_%i%02i*' % (metop,metop.lower(),year,mes)
	print datapath	
	h5filename='MEXICO/mexico_omi_no2-%i%02i.h5'% (year,mes)
	

	lista =sorted(glob.glob(datapath))
	for f in lista:
	    print f
	    print h5filename
	    read_omi_no2_domino.onefile(f,h5filename)


