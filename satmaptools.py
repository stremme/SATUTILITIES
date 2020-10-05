import numpy as np
import matplotlib.pyplot as plt
import h5py 
import datetime
import glob
try:
	from mpl_toolkits.basemap import Basemap
except:
	print '	mpl_toolkits.basemap import Basemap not working'

try:
	import shapefile
except:
	print 'shapefile no working'

def oversampelingmethod(lats,lons,filelist,name,footprintradio):
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
		
	for ifile,filename in enumerate(filelist):
		fh5=h5py.File(filename,'r')
		for ilat,lat in enumerate(lats):
			for ilon,lon in enumerate(lons):		
				dset=fh5['%iN%iW' % (int(lat),int(lon))]
				distances2=np.array((np.array(dset['lat'])-lat)**2+(np.array(dset['lon'])-lon)**2)
				conditions=(footprintradio**2 > distances2)
				index=np.where(conditions)[0]
				vec=np.array(dset[name][index])
				matrix[ilon,ilat]=matrix[ilon,ilat]+np.sum(vec)
				countermatrix[ilon,ilat]=countermatrix[ilon,ilat]+len(vec)
		fh5.close()
	matrix=matrix/countermatrix
	return matrix

def oversampeling(lats,lons,filelist,name,footprintradio,flag=''):
	data=compiledata(lats,lons,filelist,name,footprintradio)
	#print 'data:'
	#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
	stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	for ilat,lat in enumerate(lats):
		for ilon,lon in enumerate(lons):
			distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
			conditions=(footprintradio**2 > distances2)
			index=np.where(conditions)[0]
			if len(index) > 0:
				vec=data[name][index]
				matrix[ilon,ilat]=np.average(vec)
				stdmatrix[ilon,ilat]=np.std(vec)
				countermatrix[ilon,ilat]=len(vec)
			else:
				pass				
				#plt.plot(data['lon'],data['lat'],'b.')
				#plt.plot([lon],[lat],'ro')
				#plt.show()
	if flag=='full':
		errmatrix=stdmatrix/np.sqrt(countermatrix)
		return matrix.T,stdmatrix.T,errmatrix.T
	else:
		return matrix.T

def makematrixfromcompileddata(lats,lons,data,name,footprintradio,flag=''):
	
#print 'data:'
	#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
	stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	for ilat,lat in enumerate(lats):
		for ilon,lon in enumerate(lons):
			distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
			conditions=(footprintradio**2 > distances2)
			index=np.where(conditions)[0]
			if len(index) > 0:
				vec=data[name][index]
				matrix[ilon,ilat]=np.average(vec)
				stdmatrix[ilon,ilat]=np.std(vec)
				countermatrix[ilon,ilat]=len(vec)
			else:
				pass				
				#plt.plot(data['lon'],data['lat'],'b.')
				#plt.plot([lon],[lat],'ro')
				#plt.show()
	errmatrix=stdmatrix/np.sqrt(countermatrix)
	if flag=='full':
		return matrix.T,stdmatrix.T,errmatrix.T
	else:
		return matrix.T

def makematrixfromcompileddatacartesian(lats,lons,data,name,flag=''):
	dlat=abs(lats[1]-lats[0])	
	dlon=abs(lons[0]-lons[1])
#print 'data:'
	#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	if flag=='full':
		countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
		stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	for ilat,lat in enumerate(lats):
		conditionlat=(np.abs(lat-data['lat'])<dlat)
		for ilon,lon in enumerate(lons):
			conditionlon=(np.abs(lon-data['lon'])<dlon)
			conditions=(np.logical_and(conditionlon,conditionlat))
			index=np.where(conditions)[0]
			if len(index) > 0:
				vec=data[name][index]
				matrix[ilon,ilat]=np.average(vec)
				if flag=='full':
					stdmatrix[ilon,ilat]=np.std(vec)
					countermatrix[ilon,ilat]=len(vec)
			else:
				pass				
				#plt.plot(data['lon'],data['lat'],'b.')
				#plt.plot([lon],[lat],'ro')
				#plt.show()
	if flag=='full':
	
		errmatrix=stdmatrix/np.sqrt(countermatrix)	
		return matrix.T,stdmatrix.T,errmatrix.T
	else:
		return matrix.T


def oversampelingwithconditions(lats,lons,filelist,name,footprintradio,flag='',condition1=True, condition2=True):
#def oversampelingwithcondition(lats,lons,filelist,name,footprintradio,flag='',condition=True):
	#data=compiledata(lats,lons,filelist,name,footprintradio)
	#data=compiledatawithcondition(lats,lons,filelist,name,footprintradio,condition=True)
	data=compiledatawithconditions(lats,lons,filelist,name,footprintradio,condition1=True,condition2=True)
	#print 'data:'
	#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
	stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
		
	for ilat,lat in enumerate(lats):
		for ilon,lon in enumerate(lons):
			distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
			conditions=(footprintradio**2 > distances2)
			index=np.where(conditions)[0]
			if len(index) > 0:
				vec=data[name][index]
				matrix[ilon,ilat]=np.average(vec)
				stdmatrix[ilon,ilat]=np.std(vec)
				countermatrix[ilon,ilat]=len(vec)
			else:
				pass				
				#plt.plot(data['lon'],data['lat'],'b.')
				#plt.plot([lon],[lat],'ro')
				#plt.show()
	errmatrix=stdmatrix/np.sqrt(countermatrix)
	if flag=='full':
		return matrix.T,stdmatrix.T,errmatrix.T
	else:
		return matrix.T

def gettimedependence(data,name,npoly,nanual,nbienal,flag=''):
	ndata=len(data)
	nparam=npoly+2*nanual+2*nbienal
	omega1=2*np.pi/(24*365*3600)
	omega2=2*np.pi/(2*24*365*3600)
	kmat=np.zeros((ndata,nparam))
	tm=np.average(data['tepoch'])
	t=data['tepoch']-tm

	for i in range(npoly):
		kmat[:,i]=t[:]**i

	for i in range(nanual):
		kmat[:,npoly+2*i]=np.sin((i+1)*omega1*t)
		kmat[:,npoly+2*i+1]=np.cos((i+1)*omega1*t)

	for i in range(nbienal):
		kmat[:,npoly+2*nanual+2*i]=np.sin((i+1)*omega2*t)
		kmat[:,npoly+2*nanual+2*i+1]=np.cos((i+1)*omega2*t)

	KTK=np.dot(kmat.T,kmat)
	print KTK
	Gain=np.dot(np.linalg.inv(KTK),kmat.T)
 	x=np.dot(Gain,data[name])
	yfit=np.dot(kmat,x)

	if flag=='':
		return x,yfit
	elif flag=='full':

		t=np.array([datetime.datetime.utcfromtimestamp(tt) for tt in data['tepoch']])
		errormatrix=np.dot(Gain,Gain.T)*np.std((data[name]-yfit))**2
		return x,errormatrix,t,data[name],yfit,kmat


def removetimedependence(data,name,npoly,nanual,nbienal):

	x,yfit=gettimedependence(data,name,npoly,nanual,nbienal)
	#wolf remove time dependence
	
	plt.plot(data['tepoch'],data[name],'bo')
	data[name]=data[name]-yfit
	plt.plot(data['tepoch'],yfit,'ro')
	plt.show()
	return data

def oversampelingremovetime(lats,lons,filelist,name,footprintradio,npoly=1,nanual=0,nbienal=0,flag='',data=[]):
	#wolf
	if data==[]:
		data=compiledata(lats,lons,filelist,name,footprintradio)
	data=removetimedependence(data,name,npoly,nanual,nbienal)
	#print 'data:'
	#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
	matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
	countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
	stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
		
	for ilat,lat in enumerate(lats):
		for ilon,lon in enumerate(lons):
			distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
			conditions=(footprintradio**2 > distances2)
			index=np.where(conditions)[0]
			if len(index) > 0:
				vec=data[name][index]
				matrix[ilon,ilat]=np.average(vec)
				stdmatrix[ilon,ilat]=np.std(vec)
				countermatrix[ilon,ilat]=len(vec)
			else:
				#plt.plot(data['lon'],data['lat'],'b.')
				#plt.plot([lon],[lat],'ro')
				#plt.show()
				pass	
	errmatrix=stdmatrix/np.sqrt(countermatrix)
	if flag=='full':
		return matrix.T,stdmatrix.T,errmatrix.T,countermatrix.T
	else:
		return matrix.T


def compiledatacartesian(ilats,ilons,filelist,name):
	#wolf:
	satdatatype=np.dtype([(name,float),('lat',float),('lon',float),('tepoch',int)])

	for ifile,filename in enumerate(filelist):
		fh5=h5py.File(filename,'r')
		#print ifile
		for lat in ilats:
			for lon in ilons:
				#print lat,lon
				try:
					dset=fh5['%iN%iW' % (lat,lon)]
					#conditionlat=np.logical_and(dset['lat']>=minlat,dset['lat']<=maxlat)
					#conditionlon=np.logical_and(dset['lon']>=minlon,dset['lon']<=maxlon)
					#conditions=np.logical_and(conditionlat,conditionlon)
					conditions=True
					index=np.where(dset[name] > 0 )[0]
					newvec=np.zeros((len(index)),dtype=satdatatype)
					newvec['lat']=dset['lat'][index]
					newvec['lon']=dset['lon'][index]
					newvec['tepoch']=dset['tepoch'][index]
					newvec[name]=dset[name][index]
					try:
						data=np.append(data,newvec)
					except:
						#print 'exeption'
						data=np.array(newvec,copy=True)
				except:
					print 'exception dataset no exist  %iN%iW in %s' % (lat,lon,filename)

					print fh5.keys()

	
				#print 'data:',len(data)
				#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

		fh5.close()
	#plt.plot(data['lon'],data['lat'],'.')
	#plt.show()
	return data


def compiledatacartesianlog(ilats,ilons,filelist,name):
	#wolf:
	satdatatype=np.dtype([(name,float),('log',float),('lat',float),('lon',float),('tepoch',int)])

	for ifile,filename in enumerate(filelist):
		fh5=h5py.File(filename,'r')
		#print ifile
		for lat in ilats:
			for lon in ilons:
				#print lat,lon
				try:
					dset=fh5['%iN%iW' % (lat,lon)]
					#conditionlat=np.logical_and(dset['lat']>=minlat,dset['lat']<=maxlat)
					#conditionlon=np.logical_and(dset['lon']>=minlon,dset['lon']<=maxlon)
					#conditions=np.logical_and(conditionlat,conditionlon)
					conditions=True
					index=np.where(dset[name] > 0 )[0]
					newvec=np.zeros((len(index)),dtype=satdatatype)
					newvec['lat']=dset['lat'][index]
					newvec['lon']=dset['lon'][index]
					newvec['tepoch']=dset['tepoch'][index]
					newvec[name]=dset[name][index]
					newvec['log']=np.log(newvec[name])
					try:
						data=np.append(data,newvec)
					except:
						#print 'exeption'
						data=np.array(newvec,copy=True)
				except:
					print 'exception dataset no exist  %iN%iW in %s' % (lat,lon,filename)

					print fh5.keys()

	
				#print 'data:',len(data)
				#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

		fh5.close()
	#plt.plot(data['lon'],data['lat'],'.')
	#plt.show()
	return data




def compiledata(lats,lons,filelist,name,footprintradio,name2=''):
	#wolf:
	if name2=='':
		satdatatype=np.dtype([(name,float),('lat',float),('lon',float),('tepoch',int)])

	else:
		satdatatype=np.dtype([(name,float),(name2,float),('lat',float),('lon',float),('tepoch',int)])
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
					conditions=np.logical_and(conditionlat,conditionlon)
					index=np.where(conditions)[0]
					newvec=np.zeros((len(index)),dtype=satdatatype)
					newvec['lat']=dset['lat'][:][index]
					newvec['lon']=dset['lon'][:][index]
					newvec['tepoch']=dset['tepoch'][:][index]
					newvec[name]=dset[name][index]
					if name2=='':
						pass
					else:
						newvec[name2]=dset[name2][index]
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



def compilefulldata(ilat,ilon,filelist):
	#wolf:
	for ifile,filename in enumerate(filelist):
		fh5=h5py.File(filename,'r')
		print filename
		#print fh5.keys()
		try:
			dset=fh5['%iN%iW' % (ilat,ilon)]
			satdtype=dset[0].dtype
			newvec=np.zeros((len(dset)),dtype=satdtype)
			for name in satdtype.names: 
				newvec[name]=dset[name]


			try:
				data=np.append(data,newvec)
			except:
				#print 'exeption'
				data=np.array(newvec,copy=True)
		except:
			print '%iN%iW' % (ilat,ilon)
			
			print 'exception dataset no exist'
				#print 'data:',len(data)
				#print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

		fh5.close()
	#plt.plot(data['lon'],data['lat'],'.')
	#plt.show()
	return data




def compiledatawithconditions(lats,lons,filelist,name,footprintradio,condition1,condition2=True):
	#wolf, beatriz:
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
					#vecarr=dset[condition1[0]]
					vec=np.array([ v[0] for v in vecarr])
					#print vec
					#vec=np.array([v[condition1[0]]])						
					value='1' #0 para pm y 1 para am
					#value=condition1[2]
					extracondition=(vec==value)
					conditions=np.logical_and(np.logical_and(conditionlat,conditionlon),(extracondition))
					#vecarr2=dset['date']
					#vec2=np.array([v[2:4] for v in vecarr2]) #v[0:2]
					#value2='16'
					#extracondition2=(vec2==value2)
					#conditions=np.logical_and(np.logical_and(conditionlat,conditionlon),(extracondition,extracondition2))
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



def map_plot(lons,lats,matrix,cblabel='',fontsize = 10,resolution='h',alpha=0.5,show=True,contourflag='',figsize=(20,10),cbmin=0.0,cbmax=0.0,cbinter=0,suptitle='',nparalel=0.5,nmeridian=0.5,plotfig=True):
	if plotfig:
		 plt.figure(figsize=figsize)
	if suptitle=='':
		pass
	else:
		plt.title(suptitle)
	#zuleica:.................................
	latmin=min(lats)
	latmax=max(lats)
	lonmin=min(lons)
	lonmax=max(lons)
	m = Basemap(projection='cyl', resolution=resolution,
		llcrnrlat=latmin, urcrnrlat = latmax,
		llcrnrlon=lonmin, urcrnrlon = lonmax)
	m.drawcoastlines(linewidth=0.5)
	m.drawparallels(np.arange(-90., 90., nparalel), labels=[1, 0, 0, 0])
	m.drawmeridians(np.arange(-180, 180., nmeridian), labels=[0, 0, 0, 1])
	#m.drawstates(linewidth=0.5)
	m.drawcountries(linewidth=0.5)
	m.readshapefile('/home/STG4_SATELITE2/WORK/ZULEICA/ESTADOS', 'ESTADOS')
	m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
	if cbinter> 0:
		values=np.arange(cbmin,cbmax,cbinter)
		if contourflag=='contour':
			plt.contour(lons,lats,matrix, values,alpha=alpha)
		else:
			plt.contourf(lons,lats,matrix,values, alpha=alpha, cmap=plt.cm.jet)		
	else:
		
		if contourflag=='contour':
			plt.contour(lons,lats,matrix, alpha=alpha)
		else:
			plt.contourf(lons,lats,matrix, alpha=alpha, cmap=plt.cm.jet)				
	cbar=m.colorbar(location='bottom',pad="5%")
	cbar.ax.set_xlabel(cblabel, fontsize=fontsize, weight="bold")
	#plt.colorbar()
	if show:
		plt.show()



# #other colors for colormap
	#Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, 
	#Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, 
	#Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r,
	#RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, 
	#Vega10, Vega10_r, Vega20, Vega20_r, Vega20b, Vega20b_r, Vega20c, Vega20c_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr,
	#YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, 
	#coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, 
	#gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r,
	# gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, 
	#pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, 
	#summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r





def compiledatamopitt(lats,lons,filelist,name,footprintradio):
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
					conditions=np.logical_and(conditionlat,conditionlon)
					index=np.where(conditions)[0]
					newvec=np.zeros((len(index)),dtype=satdatatype)
					newvec['lat']=dset['lat'][:][index]
					newvec['lon']=dset['lon'][:][index]
					newvec['tepoch']=dset['tepoch'][:][index]
					#newvec[name]=dset[name][:][index]
					#print np.array(dset[name][:][index])[:,0].shape
					newvec[name]=np.array(dset[name][:][index])[:,0]
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



def mapplot(lons,lats,matrix,footprint, cblabel='',fontsize = 10,resolution='h',alpha=0.5,show=True,contourflag='',figsize=(8,8),cbmin=0.0,cbmax=0.0,cbinter=0,suptitle='',nparalel=0.5,nmeridian=0.5):
	plt.figure(figsize=figsize)
	plt.suptitle(suptitle)
	#zuleica:.................................
	#_________________________________
	#SCATTER
	lonsMD = [-99.176, -99.1492,  -99.24, -99.199]
	latsMD = [19.326, 19.4837 , 19.483, 19.71 ]
	
	xs, ys = np.meshgrid(lonsMD, latsMD)
	dato_unam=3.578E16
	dato_vall=3.15E16
	dato_acat=3.709E16
	dato_cuat=1.97E16
	vcdunam=dato_unam#[:].astype(np.float64)
	vcdcuat=dato_cuat#[:].astype(np.float64)
	vcdvall=dato_vall#[:].astype(np.float64)
	vcdacat=dato_acat#[:].astype(np.float64)
	data = [vcdunam, vcdvall, vcdacat, vcdcuat]
	vu=vcdunam
	vv=vcdvall
	va=vcdacat
	vc=vcdcuat

	#_____________________________________
	latmin=min(lats)
	latmax=max(lats)
	lonmin=min(lons)
	lonmax=max(lons)
	m = Basemap(projection='cyl', resolution=resolution,
		llcrnrlat=latmin, urcrnrlat = latmax,
		llcrnrlon=lonmin, urcrnrlon = lonmax)
	m.drawcoastlines(linewidth=0.5)
	m.drawparallels(np.arange(-90., 90., nparalel), labels=[1, 0, 0, 0])
	m.drawmeridians(np.arange(-180, 180., nmeridian), labels=[0, 0, 0, 1])
	#m.drawstates(linewidth=0.5)
	m.drawcountries(linewidth=0.5)
	m.readshapefile('/home/STG4_SATELITE2/WORK/ZULEICA/ESTADOS', 'ESTADOS')
	m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
	
	if cbinter> 0:
		values=np.arange(cbmin,cbmax,cbinter)
		if contourflag=='contour':
			plt.contour(lons,lats,matrix, values,alpha=alpha)
		else:
			plt.contourf(lons,lats,matrix,values, alpha=alpha, cmap=plt.cm.jet)		
	else:
		
		if contourflag=='contour':
			plt.contour(lons,lats,matrix, alpha=alpha)
		else:
			plt.contourf(lons,lats,matrix, alpha=alpha, cmap=plt.cm.jet)				
	cbar=m.colorbar(location='right',pad="5%")
	cbar.ax.set_xlabel(cblabel, fontsize=fontsize, weight="bold")
	plt.title('NO2 OMI-DOMINO '+footprint+'')
	#plt.colorbar()
	
	#___________________
	import matplotlib.cm as cm
	cm=plt.cm.get_cmap('jet') # autumn and winter
	sc=m.scatter(lonsMD,latsMD, c=data, vmin=0, vmax=3.72E16, s=1000, cmap=cm)
	cbar=m.colorbar(sc, location='bottom',pad="5%")
	#cbar.set_label(units)
	x,y = m(lonsMD,latsMD)
	#labels de cada punto de localizacion de M-DOAS
	labels = [('UNAM',vcdunam), ('VALL', vcdvall),('ACAT', vcdacat),('CUAT',vcdcuat)]
	x_offsets = [-0.15, -0.05, 0, -0.1]
	y_offsets = [-0.1, 0.1, 0.1, 0.1]
	for label, xpt, ypt, x_offset, y_offset in zip(labels,x,y, x_offsets, y_offsets):
		plt.text(xpt+x_offset,ypt+y_offset,label)
	#___________________________
	if show:
		plt.show()



