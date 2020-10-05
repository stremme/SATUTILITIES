import numpy as np
import fileinput
import h5py
import glob
import sys
import fileinput
import datetime
filename='OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2016m0819t183819.he5'
#data= np.loadtxt(filename)
latmin=13.0
lonmin=-119
latmax=34.0
lonmax=-73
year=2004

to=datetime.datetime.utcfromtimestamp(0.0)
totai93=datetime.datetime(1993,1,1,0,0,0)
dt_epochmenostai93=(to-totai93).total_seconds()

fields=[('tepoch',int),\
	('lat',float),\
	('lon',float),\
	('Time',float),\
	('SolarZenithAngle',float),\
	('SolarAzimuthAngle',float),\
	('ViewingAzimuthAngle',float),\
	('ViewingZenithAngle',float),\
	('XTrackQualityFlags',int),\
	('CloudFraction',int),\
	('ColumnAmountNO2Trop',float),\
	('ColumnAmountNO2TropStd',float),\
	('ColumnAmountNO2',float),\
	('ColumnAmountNO2Std',float),\
	('ColumnAmountNO2Strat',float),\
	('ColumnAmountNO2StratStd',float),\
	('AmfTrop',float),\
	('AmfTropStd',float),\
	('AmfStrat',float),\
	('AmfStratStd',float),\
	('AMFQualityFlags',int),\
	('GroundPixelQualityFlags',int),\
	('FoV75CornerLatitude',float,4),\
	('FoV75CornerLongitude',float,4)
	]
omino2type=np.dtype(fields)
print omino2type.names
data=np.zeros((1),dtype=omino2type)
one_measurement=data[0]

def readoneomino2swath(fhdf,itime):
	one_measurements=np.zeros((60),dtype=omino2type)
	for name in omino2type.names[1:]:
		if name == 'lat':
			dset=fhdf['HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'+'Latitude']
		elif name =='lon':
			dset=fhdf['HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'+'Longitude']
		else:
			try:
				dset=fhdf['HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'+name]
			except:
				dset=fhdf['HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'+name]
		#print 'name',name
		#print dset.shape
		#print one_measurements[name].shape
		if name=='Time':
			one_measurements[name][:]=dset[itime]
			one_measurements['tepoch'][:]=int(dset[itime]-dt_epochmenostai93)
		elif name=='FoV75CornerLatitude' or name=='FoV75CornerLongitude':
			one_measurements[name][:,:]=dset[:,itime,:].T
		else:
			one_measurements[name][:]=dset[itime,:]	
	return one_measurements

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



def onefile(filenamein,h5filenameout):
    counter=0		
    try:
	fh5=h5py.File(h5filenameout,'r+')
	print 'file exist'
    except:
	fh5=h5py.File(h5filenameout,'w')
	print 'new file'
        
    fhdfin=h5py.File(filenamein)
    dsettime=fhdfin['HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time']
    for itime in range(len(dsettime)):
	print itime
	onemeasurements=readoneomino2swath(fhdfin,itime)
	for irow in range(60):	
		onemeasurement=onemeasurements[irow]
		if latmin < onemeasurement['lat'] < latmax and lonmin < onemeasurement['lon'] < lonmax:
		    print onemeasurement['lat'],onemeasurement['lon']
		    if onemeasurement['ColumnAmountNO2Trop'] > -1.0E30:
			addonemeasurement(onemeasurement,fh5)
		
    fhdfin.close()	
    fh5.close()
