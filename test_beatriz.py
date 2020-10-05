import numpy as np
import matplotlib.pyplot as plt
import h5py 
import datetime
import glob
import sys
import pickle
sys.path.append('/home/STG_04/WORK/SATUTILITIES/')
import satmaptools

#path='/home/STG_04/CO_IASI_A/FORLI/MEXICO/*.h5'
path='/home/STG_04/NH3_IASI_B/MEXICO/*.h5'   ##first version pole
#path='/home/D1_SATELITE/IASI/NH3_IASI_A/MEXICO/v2_1/*.h5' ## v2_1 espri
#path='/home/D1_SATELITE/IASI/NH3_IASI_A/MEXICO/*.h5'  ##v2_2 espri
lista=glob.glob(path)
print lista
nx=50
ny=50

latmin=24.0
latmax=26.0
lonmin=-101.0
lonmax=-99.0

latmin=18.4
latmax=20.2
lonmin=-99.8
lonmax=-98.0


latmin=17
latmax=22
lonmin=-99.8
lonmax=-98
'''
latmin=17.0
lonmin=-118.0
latmax=34.0
lonmax=-73.0
'''

lats=latmin+(latmax-latmin)*np.arange(ny)/ny
lons=lonmin+(lonmax-lonmin)*np.arange(nx)/nx
#footprintcircle=6.0/110.0
footprintcircle=6.0/110.0

#matrix= satmaptools.oversampeling(lats,lons,lista,'co_tot',footprintcircle)
#matrix= satmaptools.oversampeling(lats,lons,lista,'column',footprintcircle)
#matrix= satmaptools.oversampelingwithconditions(lats,lons,lista,'column',footprintcircle) #funciona con el compiledatawithcondition sin modificar
#matrix= satmaptools.oversampelingwithconditions(lats,lons,lista,'column',footprintcircle,condition1=('time','0','1')) #'column'
#matrix= satmaptools.makematrixfromcompileddata(lats,lons,lista,'column',footprintcircle)
#matrix= satmaptools.oversampelingremovetime(lats,lons,lista,'co_tot',footprintcircle,npoly=3,nanual=4,nbienal=1,flag='')
#satmaptools.mapplot(lons,lats,matrix,cblabel='NH3_tot [molec/cm2]',cbmin=0.0,cbmax=1.5E16,cbinter=0.2E16) 	
#plt.contourf(lons,lats,matrix)
#plt.colorbar()
#plt.show()
#pickle.dump(matrix, open("matrizCO_A_am.p","wb"))

###ajuste de zuleica


data=satmaptools.compiledata(lats,lons,lista,'column',footprintcircle)
plt.plot(data['column'],'bo')
index=np.where(data['column'] > 0.0)[0]
data=data[index]
plt.plot(data['column'],'rx')
plt.show()

matrix= satmaptools.makematrixfromcompileddata(lats,lons,data,'column',footprintcircle)
#matrix= satmaptools.oversampeling(lats,lons,lista,'TroposphericVerticalColumn',footprintcircle)
satmaptools.mapplot(lons,lats,matrix,cblabel='NH3_tot [molec/cm2]',cbmin=0.0,cbmax=1.5E16,cbinter=0.2E16) 	
#satmaptools.mapplot(lons,lats,matrix,cblabel='NO2_tot [molec/cm2]',cbmin=0.0,cbmax=2E18,cbinter=1.0E17) 
#plt.contourf(lons,lats,matrix)
#plt.colorbar()
plt.show()

'''
def compiledatawithcondition(lats,lons,filelist,name,footprintradio,condition=True):
	#wolf:
	satdatatype=np.dtype([(name,float),('lat',float),('lon',float),('tepoch',int)])

	minlat=(min(lats)-footprintradio/110.0)
	maxlat=(max(lats)+footprintradio/110.0)
	minlon=(min(lons)-footprintradio/110.0)
	maxlon=(max(lons)+footprintradio/110.0)
	#print minlat,maxlat,minlon,maxlon
	ilats=np.arange(int(maxlat)-int(minlat)+1,dtype=int)+int(minlat)
	ilons=np.arange(int(maxlon)-int(minlon)+1,dtype=int)+int(minlon)
	

	for ifile,filename in enumerate(filelist):
		fh5=h5py.File(filename,'r')
		#print ifile
		for lat in ilats:
			for lon in ilons:
				#print lat,lon
				try:
					dset=fh5['%iN%iW' % (int(lat),int(lon))]
					conditionlat=np.logical_and(dset['lat']>=minlat,dset['lat']<=maxlat)
					conditionlon=np.logical_and(dset['lon']>=minlon,dset['lon']<=maxlon)
					vecarr=dset['time']
					#print vecarr
					vec=np.array([ v[0] for v in vecarr])
					#print vec						
					value='1'
					extracondition=(vec==value)
					conditions=np.logical_and(np.logical_and(conditionlat,conditionlon),extracondition)
					index=np.where(conditions)[0]
					newvec=np.zeros((len(index)),dtype=satdatatype)
					newvec['lat']=dset['lat'][:][index]
					newvec['lon']=dset['lon'][:][index]
					newvec['tepoch']=dset['tepoch'][:][index]
					newvec[name]=dset[name][index]
					try:
						data=np.append(data,newvec)
					except:
						#print 'exeption'
						data=np.array(newvec,copy=True)
				except:
					print 'exception dataset no exist'
				#print 'data:',len(data)
				#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

		fh5.close()
	#plt.plot(data['lon'],data['lat'],'.')
	#plt.show()
	return data

'''
