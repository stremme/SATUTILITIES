import numpy as np
import fileinput
import h5py
import glob
import sys
import fileinput
import readiasicoh5B


year=2014




meses=np.arange(12,dtype=int)+1
for 	year in [2013,2014,2015,2016,2017]:
	for 	mes in meses:
		metop='B'
		datapath='txt/iasi_CO_LATMOS_ULB_metop%s_%i%02i*' % (metop.lower(),year,mes)
		#datapath='%s-data/iasi_CO_LATMOS_ULB_metop%s_%i%02i*' % (metop,metop.lower(),year,mes)
		print datapath	
		h5filename='MEXICO/mexico_iasi_%s-%i%02i.h5'% (metop,year,mes)
	

		lista =sorted(glob.glob(datapath))
		for f in lista:
		    print f
		    print h5filename
		    readiasicoh5B.onefile(f,h5filename,metop,year,mes)
