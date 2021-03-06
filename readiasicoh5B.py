import numpy as np
import fileinput
import h5py
import glob
import sys
import fileinput
import datetime
filename='iasi_CO_LATMOS_ULB_metopb_20131030_v20100815.txt'
#data= np.loadtxt(filename)
latmin=13.0
lonmin=-119
latmax=34.0
lonmax=-73


to=datetime.datetime.utcfromtimestamp(0.0)


fields=[('tepoch',int),\
	('lat',float),\
	('lon',float),\
	('date',str,8),\
	('time',str,6),\
	('sza',float),\
	('iFOV',int),\
	('qflag',int,10),\
	('cloud_coverage',float),\
	('dof',float),\
	('rms',float),\
	('fitbias',float),\
	('co_tot',float),\
	('co_err',float),\
	('xa',float,19),\
	('AK_tot',float,19)
	]
iasicotype=np.dtype(fields)
print iasicotype.names
data=np.zeros((1),dtype=iasicotype)
one_measurement=data[0]

def readoneiasi(linearr):
	print linearr
	one_measurement=np.zeros((1),dtype=iasicotype)[0]
	icounter=0
	for name in iasicotype.names[1:]:
		print name
		nlen=1	
		try:
			nlen=len(one_measurement[name])
		except: 
			nlen=1

		print name,linearr[icounter]
		if nlen == 1 or one_measurement[name].dtype == 'str':
			nlen=1
			if one_measurement[name].dtype == 'float':
				val=float(linearr[icounter])				
			if one_measurement[name].dtype == 'int':
				val=int(linearr[icounter])
			if one_measurement[name].dtype == 'str':
				val= linearr[icounter]
				if name == 'time':
					val='%06i' % int(linearr[icounter])
			one_measurement[name]=val
		else:
			print icounter,nlen,name
			subarr=np.array(linearr[icounter:icounter+nlen],dtype=one_measurement[name].dtype )
			one_measurement[name][:]=subarr[:]

		icounter=icounter+nlen
	try:
		one_measurement['tepoch']=(datetime.datetime.strptime(one_measurement['date']+one_measurement['time']+" UTC",'%Y%m%D%H%M%S %Z')-to).total_seconds()
	except: 
		print 'problems'
		year=int(one_measurement['date'][0:4])
		month=int(one_measurement['date'][4:6])
		day=int(one_measurement['date'][6:8])
		hour=int(one_measurement['time'][0:2])
		minute=int(one_measurement['time'][2:4])
		second=int(one_measurement['time'][4:6])
		if second == 60:
			minute=minute+1
			second=0
		if minute == 60:
			minute=0
			hour=hour+1
		one_measurement['tepoch']=(datetime.datetime(year,month,day,hour,minute,second)-to).total_seconds()

	print one_measurement['lat'],one_measurement['lon']
	
	return one_measurement



def addonemeasurement(onemeasurement,fh5):
	print 'addmeasurements'
	lat=onemeasurement['lat']
	lon=onemeasurement['lon']
	print int(lat),int(lon)



	try:
		dset=fh5['%iN%iW' % (int(lat),int(lon))]
		#print 'dset exist'
		#print dset
		#print dset.shape
		n=dset.shape[0]
		#print n
		dset.resize((n+1,))
		dset[n]=onemeasurement
		fh5.flush()
	except:
		#print 'exception'
		datalatlon=[onemeasurement]
		try:
			dset=fh5.create_dataset('%iN%iW' % (int(lat),int(lon)),data=datalatlon,maxshape=(None,))

			fh5.flush()
		except:
			print 'nada works'



def onefile(filename,h5filename,metop,year,mes):
    counter=0		
    try:
	fh5=h5py.File(h5filename,'r+')
	print 'file exist'
    except:
	fh5=h5py.File(h5filename,'w')
	print 'new file'
    datos=np.genfromtxt(filename,dtype=str,usecols=np.arange(60),skip_header=1,skip_footer=1,autostrip=True)		
    for line in datos:
	#print line[0]
	counter=counter+1
	if np.mod(counter,10000) ==0:
	    print counter#,len(datos)
		#break
	if latmin < float(line[0]) < latmax and lonmin < float(line[1]) < lonmax:
	    print float(line[0]),float(line[0])
	    onemeasurement=readoneiasi(line)
	    print onemeasurement['lat'],onemeasurement['lon']
	    addonemeasurement(onemeasurement,fh5)
    zlevels=np.arange(20)
    fh5.attrs.create('zlevels', zlevels, dtype=zlevels.dtype )
    fh5.close()




