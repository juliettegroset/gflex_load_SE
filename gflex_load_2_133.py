
def gflex(PATH,NAME,DATA,t,tau_time,gFlex,syst,point_gps,NAME_GPS,grid_inp,addload,uniform_grid,p,rhoq,NB) : 

### Library, etc.
	import gflex
	import numpy as np
	from matplotlib import pyplot as plt
	from matplotlib import colors 
	from matplotlib import cm
	import pandas as pd
	from scipy import linalg
	from pyproj import Transformer

	flex = gflex.F2D()
	flex.Quiet = False
	flex.Debug=True
	flex.Verbose=False

	## FLAGS ########################################

	PLOT = False
	plot_with_gflex = False

	# Output files
	W_file_grd = True
	Epsilon_sigma_file = True

	# type de calcul
	tau_type = 'Time'       # Relaxation visqueuse definie par temps caracteristique ou viscosite 
	rake_input = False

	## INPUT ###########################################

	if NB == 1 :
		DIR_OUT = PATH+'output/'+NAME+'/' + NAME + '_Te' + str(int(gFlex[6]*0.001)) + 'km_tau' + str(tau_time) +'yr/'  # Repertoire de sortie
	else : 
		DIR_OUT = PATH+'output/RANDOM/'+NAME+'/' + NAME + '_Te' + str(int(gFlex[6]*0.001)) + 'km_tau' + str(tau_time) +'yr/'  # Repertoire de sortie

	# Temps caractéristique
	if tau_type == 'Time' : 
		tau = tau_time #yr

	if tau_type == 'Visc':
		mu = 1e21       # viscosité Pa.s
		l = 150000      # longueur d'onde de la charge (m)
		tau = (4*np.pi*mu)/(flex.rho_m*9.81*l)

	# Conversions between WGS84 and UTM32
	if syst == 'UTM32' : 
		wgs2utm = Transformer.from_crs("epsg:4326","epsg:32632")
		utm2wgs = Transformer.from_crs("epsg:32632","epsg:4326")
	elif syst == 'UTM31' :
		wgs2utm = Transformer.from_crs("epsg:4326","epsg:32631")
		utm2wgs = Transformer.from_crs("epsg:32631","epsg:4326")
	elif syst == 'UTM19' :
		wgs2utm = Transformer.from_crs("epsg:4326","epsg:32619")
		utm2wgs = Transformer.from_crs("epsg:32619","epsg:4326")
	else : 
		print('Système de coordonnées non valide')
		

	if rake_input :
		ri = 0

	## GFLEX INPUT ###########################################

	print(' \n ')

	flex.Method = gFlex[0] # Solution method: * FD (finite difference)
                    #                  * SAS (superposition of analytical solutions)
                    #                  * SAS_NG (ungridded SAS)
	flex.PlateSolutionType = 'vWC1994' # van Wees and Cloetingh (1994)
                                   # The other option is 'G2009': Govers et al. (2009)
	flex.Solver = 'direct' # direct or iterative
#convergence = 1E-3 # convergence between iterations, if an iterative solution method is chosen

	flex.g = gFlex[1]          # acceleration due to gravity
	flex.E = gFlex[2]         # Young's Modulus
	flex.nu = gFlex[3]        # Poisson's Ratio
	flex.rho_m = gFlex[4]    # MantleDensity
	flex.rho_fill = gFlex[5] # InfiillMaterialDensity

	flex.Te = gFlex[6]       # Elastic Thickness (m)

	##### Pour la différence finis 
	# Boundary conditions can be:
	# (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
	# For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
	flex.BC_W = '' # west boundary condition
	flex.BC_E = '' # east boundary condition
	flex.BC_S = '' # south boundary condition
	flex.BC_N = '' # nor	th boundary condition

	z1 = flex.Te/2

	## TRAITEMENT PRE-GLFLEX ###########################################
	## NB : gflex computation results are in (x,y,z) cartesian coordinates
	## with x = East, y = North, z = Up (flexure downward is negative)
	####################################################################

	# Chargement du fichier d'entrée
	L = pd.read_table(DATA, delim_whitespace=True, header=6)	
	L = L * gFlex[7]
	
	ASC = open(DATA,'r')
	line = ASC.readlines()
	
	if uniform_grid:
		flex.dx = float(line[4].split()[1]) # grid cell size, x-oriented [m]	
		flex.dy = float(line[4].split()[1]) # grid cell size, x-oriented [m]	
	else:
		flex.dx = float(line[4].split()[1]) # grid cell size, x-oriented [m]	
		flex.dy = float(line[5].split()[1]) # grid cell size, x-oriented [m]	

	# definition de la grille total (UTM32, km)
	if grid_inp[0] == 'NONE' or grid_inp[1] == 'NONE' or grid_inp[2] == 'NONE' or grid_inp[3] == 'NONE' :
		xminmax = [-100,1200]
		yminmax = [4500,5600]
		print('Le calcul est réalisé sur une grille de coordonnées {} {} {}'.format(syst,xminmax,yminmax))
	else :
		xmin, ymin = wgs2utm.transform(grid_inp[0],grid_inp[2])
		xmax, ymax = wgs2utm.transform(grid_inp[1],grid_inp[3])
		xminmax = [int(xmin*0.001),int(xmax*0.001)]
		yminmax = [int(ymin*0.001),int(ymax*0.001)]

	# coordonnées UTM de la position et taille du fichier d'entrée (cf .asc file) en km
	xllcorner = int(np.round((float((line[2].split())[1]))*0.001))
	yllcorner = int(np.round((float((line[3].split())[1]))*0.001))

	ncol = int(line[0].split()[1])
	
	nrow = int(line[1].split()[1])-1 
	
	# Paramètre de lamé
	mu_l = flex.E/(2*(1+flex.nu))
	lamb_l = (flex.E*flex.nu)/((1-2*flex.nu)*(1+flex.nu))

	# Initialisation
	load = np.zeros((nrow,ncol))
	l2 = np.array((L)) #conversion en matrice

	for i in range(1,nrow) : 
		load[i,:]=l2[nrow-i,:] #inversion de la matrice selon y

	# Recalage de la charge sur le modèle globale
	ncellx = ((xminmax[1]-xminmax[0])/(flex.dx/1000))+1		
	ncelly =  ((yminmax[1]-yminmax[0])/(flex.dy/1000))+1

	xpres = (xllcorner-xminmax[0])/(flex.dx/1000)
	ypres = (yllcorner-yminmax[0])/(flex.dy/1000)

	print(' \n ')


	# calcul du vecteur temps
	coeftime = 3.1536e7

	## GLFLEX ###########################################
	#initialisation 
	flex.qs = np.zeros((int(ncelly),int(ncellx)))

	Exx1 = np.zeros((int(ncelly),int(ncellx)))
	Eyy1 = np.zeros((int(ncelly),int(ncellx)))
	Exy1 = np.zeros((int(ncelly),int(ncellx)))
	Exx2 = np.zeros((int(ncelly),int(ncellx)))
	Eyy2 = np.zeros((int(ncelly),int(ncellx)))
	Exy2 = np.zeros((int(ncelly),int(ncellx)))

	psave = p # pas d'enregistrement dans le fichier de sortie (si psave = 40, un point enregistré tous les 40 points)
	j = 0
	D = [10000,10000]

	p = load * rhoq * flex.g
	flex.qs[int(ypres):int(ypres)+nrow,int(xpres):int(xpres)+ncol] = p 
	
	# latitude/longitude solutions are exact for SAS, approximate otherwise
	#latlon = # true/false: flag to enable lat/lon input. Defaults False.
	#PlanetaryRadius = # radius of planet [m], for lat/lon solutions
		
	# run gflex 
	flex.initialize()
	flex.run()              ## Calcul de la flexure 
	flex.finalize()
	
	# calcul de la flexure en fonction du temps d'après Geodynamics (Turcotte et Schubert)
	w = flex.w * np.exp(-t/tau)
	v = (-(w/tau)) * 1e3

	############ écriture du fichier de sortie + calcul des tenseurs de contrainte pour gagner du temps de calcul
	
	# Ouverture des fichiers de sortie
	print('écriture des fichiers de sorties ... \n \n')
	if W_file_grd : 
		newfile = open(DIR_OUT+'gflex_U_'+str(t)+'yr.txt','w+')
		newfile.write('Lon Lat w(m) Ux(m) Uy(m) Vup(mm/a) Vx(mm/a) Vy(mm/a) \n')
		
	if Epsilon_sigma_file :
		FILE1 = open(DIR_OUT +'gflex_E_S_'+str(t)+'yr.txt','w+')
		FILE3 = open(DIR_OUT +'/gflex_S-vpropre_'+str(t)+'yr.txt','w+')
		FILE4 = open(DIR_OUT +'gflex_Epoint_'+str(t)+'yr.txt','w+')
		FILE1.write('Lon Lat Exx Exy Eyy Sxx(MPa) Sxy(MPa) Syy(MPa) Szz(MPa)\n')
		FILE3.write('Lon lat S1(MPa) Az1(deg) S2(MPa) Az2(deg)\n')
		FILE4.write('Lon Lat exx(yr-1) exy(yr-1) eyy(yr-1) eh1(yr-1) eh2(yr-1) Az1(deg) Az2(deg) Emax(yr-1) \n')
	
	# calcul des déplacements et des vitesses horizontales
	Dy, Dx = np.gradient(w, flex.dx, flex.dy) 
	Ux, Uy = -z1*Dx , -z1*Dy
	Vx, Vy = (-Ux/tau)*1e3, (-Uy/tau) * 1e3

	# Calcul de la déformation et contrainte
	# Strain convention : e > 0  <=>  extension (length increase in positive x / y direction)
	# Stress convention : s > 0  <=>  tension (not meca convetion !)
	Dxy, Dxx = np.gradient(Dx, flex.dx, flex.dy)
	Dyy, Dyx = np.gradient(Dy, flex.dx, flex.dy)
	
	Exx1 = -z1 * Dxx
	Exy1 = -z1 * 1/2 * (Dxy + Dyx)
	Eyy1 = -z1 * Dyy

	# Bending stress from gflex	
	Sxx = -z1 * flex.E/(1-flex.nu**2) * (Dxx + flex.nu*Dyy) * 1e-6
	Syy = -z1 * flex.E/(1-flex.nu**2) * (Dyy + flex.nu*Dxx) * 1e-6
	Sxy = -z1 * flex.E/(1+flex.nu**2) * Dxy * 1e-6
	Szz = np.zeros(np.shape(Sxx))
	# Add stress from load (assuming confined medium)
	# warning : stress convention s > 0 = tension, signs must be reversed to add load effect
	if addload :
		Szz =  -flex.qs * 1e-6
		Sxx = Sxx + flex.nu/(1-flex.nu) * Szz
		Syy = Syy + flex.nu/(1-flex.nu) * Szz

	x,y = np.zeros(len(w[0,:])*len(w[:,0])), np.zeros(len(w[0,:])*len(w[:,0]))
	step = 0

			
	# calcul pour chaque point de la grille 
	rate_temp = 0
	for h in range(len(w[:,0])):	
		for k in range(len(w[0,:])):
			rate = int(h/len(w[:,0])*100)	
			if (rate != rate_temp) :
				print(str(rate) + '%')
			rate_temp = np.copy(rate)
					
			E1 = np.array(([Exx1[h,k],Exy1[h,k],0],[Exy1[h,k],Eyy1[h,k],0],[0,0,0]))
			S1 = np.array(([Sxx[h,k], Sxy[h,k], 0], [Sxy[h,k], Syy[h,k], 0], [0, 0, Szz[h,k]]))
		#	Sh = np.array(([Syy[h,k], Sxy[h,k]], [Sxy[h,k], Sxx[h,k]]))
			Sh = np.array(([Sxx[h,k], Sxy[h,k]], [Sxy[h,k], Syy[h,k]]))
			
			# vecteurs propres horizontaux
			Sw, Sv = linalg.eig(Sh)
			Sh1 = Sw.max()
			Sh2 = Sw.min()
		#	Az1 = np.rad2deg(np.arctan2(float(Sv[1,Sw.argmax()]), float(Sv[0,Sw.argmax()])))
			Az1 = np.rad2deg(np.arctan2(float(Sv[0,Sw.argmax()]), float(Sv[1,Sw.argmax()])))
			Az2 = 90 + Az1
				
			# taux de déformation horizontal
			e = -(1/tau) * np.array(([Exx1[h,k],Exy1[h,k]],[Exy1[h,k],Eyy1[h,k]]))
		#	Swe, Sve = linalg.eig(np.array(([e[1,1],e[0,1]],[e[0,1],e[0,0]])))
			Swe, Sve = linalg.eig(np.array(([e[0,0],e[0,1]],[e[0,1],e[1,1]])))
			eh1, eh2 = np.max(Swe), np.min(Swe)
		#	eaz1 = np.rad2deg(np.arctan2(Sve[1,Swe.argmax()], Sve[0,Swe.argmax()]))	
			eaz1 = np.rad2deg(np.arctan2(Sve[0,Swe.argmax()], Sve[1,Swe.argmax()]))	
			eaz2 = eaz1 + 90
				
			xj = xminmax[0]*1e3 + k*flex.dx
			yj = yminmax[0]*1e3 + h*flex.dy
			
			x[step],y[step] = xj,yj
			
			# Enregistrement des fichiers de sorties
			if j == psave :
				lat,lon = utm2wgs.transform(xj,yj)
				if W_file_grd :
					newfile.write(str(np.round(lon,3)) + ' ' + str(np.round(lat,3)) + ' ' + str(w[h,k]) + ' ' + str(np.round(Ux[h,k],3)) + ' ' + str(np.round(Uy[h,k],3)) + ' ' + str(v[h,k])+ ' ' + str(Vx[h,k])+ ' ' + str(Vy[h,k]) + '\n')

				if Epsilon_sigma_file :
					FILE1.write(str(np.round(lon,3)) + ' ' + str(np.round(lat,3)) + ' ' + str(Exx1[h,k]) + ' ' + str(Exy1[h,k]) + ' ' + str(Eyy1[h,k]) + ' ' + str(S1[0,0]) + ' ' + str(S1[0,1]) + ' ' + str(S1[1,1])  + '	' + str(S1[2,2]) +'\n')
					FILE3.write(str(np.round(lon,3)) + ' ' + str(np.round(lat,3)) + ' ' + str(Sh1.real) + ' ' + str(Az1) + ' ' + str(Sh2.real) + ' ' + str(Az2) + '\n')
					FILE4.write(str(np.round(lon,3)) + ' ' + str(np.round(lat,3)) + ' ' + str(e[0,0]) + ' ' + str(e[0,1]) + ' ' + str(e[1,1]) + ' ' + str(eh1.real) + ' ' + str(eh2.real) + ' ' + str(eaz1) + ' ' + str(eaz2) + ' ' + str(np.max(np.abs(eh1) and np.abs(eh2) and np.abs(eh1+eh2))) + '\n')
				j = 0
			
			else : 
				j = j+1
			
			step = step+1

	wt = np.array(np.reshape(w, (len(w[:,0])*len(w[0,:]),1)))
	vt = np.array(np.reshape(v, (len(v[:,0])*len(v[0,:]),1)))
	vxt = np.array(np.reshape(Vx, (len(Vx[:,0])*len(Vx[0,:]),1)))
	vyt = np.array(np.reshape(Vy, (len(Vy[:,0])*len(Vy[0,:]),1)))
	Sxxt = np.reshape(Sxx, (len(v[:,0])*len(v[0,:]),1))
	Sxyt = np.reshape(Sxy, (len(v[:,0])*len(v[0,:]),1))
	Syyt = np.reshape(Syy, (len(v[:,0])*len(v[0,:]),1))
	Szzt = np.reshape(Szz, (len(v[:,0])*len(v[0,:]),1))

	# calcul aux points gps
	if point_gps : 
		print("\ncalcul des vitesses aux points gps ... \n")
		gps = np.array(pd.read_table(NAME_GPS,header=None,delim_whitespace=True))
		model = np.stack((x,y),axis=1)
		xgps,ygps = wgs2utm.transform(gps[:,1],gps[:,2])
		FILE5 = open(DIR_OUT + 'gflex_V_pointgps_'+str(t)+'yr.txt','w+')
		FILE5.write('NAME	LAT	LON	VY	VX	VZ	0	0	0 \n')
		
		rate_temp=0
		for g in range(len(gps[:,0])) :
			rate = int((g/len(gps[:,0]))*100)
			if (rate != rate_temp) :
				print(str(rate) + '%')
			rate_temp = np.copy(rate)
			pointgps = np.array([xgps[g],ygps[g]])
			id = np.array([np.sqrt(xg**2+yg**2) for (xg,yg) in (model[:,[0,1]] - pointgps)]).argmin() 
			latg, long = utm2wgs.transform(x[id],y[id])
			FILE5.write(str(gps[g,0])+' '+str(latg)+' '+str(long)+' '+str(float(vyt[id]))+' '+str(float(vxt[id]))+' '+str(float(vt[id]))+' 0 0 0 \n')


	## Plot figure if PLOT = TRUE || OBSOLETE
	if PLOT :
		print('réalisation de la figure pour t= '+str(t) + 'yrs \n')
		X = np.arange(xminmax[0],xminmax[1],flex.dx/1000)
		Y = np.arange(yminmax[0],yminmax[1],flex.dy/1000)
	
		X,Y = np.meshgrid(X,Y)

		fig = plt.figure()
		plt.subplot(2,1,1)
		norm = colors.Normalize(vmin=np.min(w),vmax=np.max(w))
		plt.scatter(X,Y,c=w,cmap='RdBu',s=5)
		plt.scatter(xfault/1000,yfault/1000,c='black',s=15)
		plt.scatter((xfault/1000)+D[0],(yfault/1000)+D[1],c='yellow',s=10, label = 'position de la faille')
		clb1 = plt.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'))
		clb1.set_label('m')
		plt.title('Flexure pour t= '+str(t)+' yr en m')
		plt.legend(loc=0)

		plt.subplot(2,1,2)
		norm = colors.Normalize(vmin=np.min(v),vmax=np.max(v))
		plt.scatter(X,Y,c=v,cmap='RdBu',s=5)
		plt.scatter(xfault/1000,yfault/1000,c='black',s=15)
		plt.scatter((xfault/1000)+D[0],(yfault/1000)+D[1],c='yellow',s=10, label = 'position de la faille')
		clb2 = plt.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'))
		clb2.set_label('mm/an')
		plt.title('Vitesses verticales pour t= '+str(t)+' yr en mm/an')
		plt.tight_layout()
		plt.legend(loc=0)
		fig.savefig('flexure_'+str(t)+'yr')
		
	print('\n \nFlexure max = '+ str(np.max(np.abs(w))) + ' m')
	print('DONE for t = '+str(t)+ 'yrs \n ############################################ \n ')
	
	## OUTPUT ###########################################

	## OUPUT GFLEX 
	if plot_with_gflex :
	# If you want to plot the output
		flex.plotChoice='both'
	# An output file for deflections could also be defined here
	# flex.wOutFile = 
		flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plotChoice and/or flex.wOutFile have been set
	# TO OBTAIN OUTPUT DIRECTLY IN PYTHON, you can assign the internal variable,
	# flex.w, to another variable -- or as an element in a list if you are looping
	# over many runs of gFlex:
	#deflection = flex.w
	#print('Flexure max = '+str(np.max(np.abs(flex.w)))+' m \n \n ')

	
def FAULTS (PATH,NAME,tau_time,Te,syst,fault,t) :
	import numpy as np
	import pandas as pd
	from pyproj import Transformer
	
	dip = np.deg2rad(fault[1])
	
	Angle_friction = False
	A_friction = 30 # si angle de friction = True

	if Angle_friction :
		CFE = np.tan(np.deg2rad(fault[3]))
	else :
		CFE = fault[3]

	print('CFE = ' + str(CFE))
	
	DIR_OUT = PATH+'output/'+NAME+'/' + NAME + '_Te' + str(Te) + 'km_tau' + str(tau_time) +'yr/'  # Repertoire de sortie

	Sfile = pd.read_table(DIR_OUT +'gflex_E_S_'+str(t)+'yr.txt',delim_whitespace=True,header=0)
	Wfile = pd.read_table(DIR_OUT+'gflex_U_'+str(t)+'yr.txt',delim_whitespace=True,header=0)
	lon,lat,w,v = np.array(Wfile['Lon']),np.array(Wfile['Lat']),np.array(Wfile['w(m)']),np.array(Wfile['Vup(mm/a)'])
	Sxx,Sxy,Syy,Szz = np.array(Sfile['Sxx(MPa)']),np.array(Sfile['Sxy(MPa)']),np.array(Sfile['Syy(MPa)']),np.array(Sfile['Szz(MPa)'])

	if syst == 'UTM32' : 
		wgs2utm = Transformer.from_crs("epsg:4326","epsg:32632")
		utm2wgs = Transformer.from_crs("epsg:32632","epsg:4326")
	elif syst == 'UTM31' :
		wgs2utm = Transformer.from_crs("epsg:4326","epsg:32631")
		utm2wgs = Transformer.from_crs("epsg:32631","epsg:4326")
	
	x,y = wgs2utm.transform(lat,lon)
	
	if fault[0] != 'NONE' : 
		
		coord = np.array(pd.read_table(PATH+'failles/coord_'+fault[0],header=None,delim_whitespace=True))
		fault_lat = coord[:,1]
		fault_lon = coord[:,0]
		xf,yf = wgs2utm.transform(fault_lat,fault_lon)


		xfault,yfault,az,Dcsf,r = np.zeros(len(xf)-1), np.zeros(len(yf)-1), np.zeros(len(yf)-1), np.zeros(len(yf)-1), np.zeros(len(yf)-1)
		U = np.array(np.stack((x,y),axis=1))
		FILE6 = open(DIR_OUT+'/'+fault[0]+'/'+NAME+'_faille'+fault[0]+'_dip'+str(fault[1])+'_DN-'+str(fault[2])+'_frict'+str(fault[3])+'.txt','w+')		
		FILE6.write('lon lat W(m) Vup(mm/a) azimuth(deg) CFS(MPa) rake(deg) \n')
		for i in range(len(xfault)):
			xfault[i],yfault[i] = (np.abs(xf[i+1]-xf[i])/2)+np.min([xf[i],xf[i+1]]), (np.abs(yf[i+1]-yf[i])/2)+np.min([yf[i],yf[i+1]])
			az[i] = -(np.arctan2((yf[i+1]-yf[i]),(xf[i+1]-xf[i])))+(np.pi/2)
			if az[i] > np.pi :
				az[i] = az[i] - np.pi
			if az[i] < 0 :
				az[i] = az[i] + np.pi
			if fault[2] :
				az[i] = az[i] + np.pi
			
			latfault,lonfault = utm2wgs.transform(xfault[i],yfault[i])
			print('Point {:6.2f}, {:6.2f}'.format(latfault, lonfault))
			print('  fault azimuth = {:4.1f} deg.N'.format(np.rad2deg(az[i])))
			pointi = [xfault[i],yfault[i]]
			id = np.array([np.sqrt(xi**2+yi**2) for (xi,yi) in (U[:,[0,1]] - pointi)]).argmin()
			print('  w = {:4.1f} m'.format(float(w[id])))
			fd  = ([np.cos(az[i])*np.cos(dip), -np.sin(az[i])*np.cos(dip), -np.sin(dip)])
			fs  = ([np.sin(az[i]), np.cos(az[i]), 0])
			fn  = ([np.cos(az[i])*np.sin(dip), -np.sin(az[i])*np.sin(dip), np.cos(dip)])

			S = np.array(([Sxx[id], Sxy[id], 0], [Sxy[id], Syy[id], 0], [0, 0, Szz[id]]))
			print('  stress tensor : ', S)
			V1  = np.dot(S,fn)
			Sn1 = np.dot(V1,fn)
			Td1 = np.dot(V1,fd)
			Ts1 = np.dot(V1,fs)
			Tt1 = np.sqrt(Ts1**2 + Td1**2)
		#	print('  Td, Ts, Tt : ', Td1, Ts1, Tt1)
			r[i]  = np.rad2deg(np.arctan2(-Td1, Ts1))
			Dcsf[i] = Tt1 - CFE * (-1*Sn1) 
			print('  CFS = {:7.2f} MPa, rake = {:6.1f} deg'.format(Dcsf[i],r[i]))

			FILE6.write(str(lonfault)+ ' ' +str(latfault)+ ' ' + str(np.round(float(w[id]),3))+ ' ' +str(np.round(float(v[id]),3))+ ' ' +str(np.round(np.rad2deg(az[i]),2))+ ' ' +str(np.round(Dcsf[i],3))+ ' ' +str(np.round(r[i],2))+ '\n')
