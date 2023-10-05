from gflex_load_2_133 import gflex,FAULTS
import os
import numpy as np
import pandas as pd

FILE = 'GIA1_2-148'

PATH = '/Users/juliettegrosset/work/GIA_failles/' 


NAME_ini = FILE

DATA_ini = PATH + 'data/mey_3km.asc'      # Input data


gps_point = True
GPS = '/Users/juliettegrosset/work/GPS/Alpes/comp_vel/final/Blewitt_final.txt'


NAME_fault = 'NONE' #NONE if unused
dip_left = [True]
dip = [60,70,80]
friction = [0.1,0.6]

ADDLOAD = False

grid_coord = ['NONE','NONE','NONE','NONE']
#grid_coord = [38,52,-6,20]
#grid_coord = [47,51,-72,-66]  # NONE ou coordonnées de la grille en WGS84 (latmin, latmax, lonmin,lonmax)
#Si ValueError: could not broadcast input array from shape = augmenter la grille (grid_coord)
systeme_coord  = 'UTM32'
#systeme_coord  = 'UTM31'
uniform_grid=True # True if (dx = dy)


## gflex input 

tfin = 15000   # t end(yr)           
tau =  [3000,3500,4000,4500,5000,5500,6000,6500,7000] #relaxation time (yr)

method = 'SAS'
g = 9.81        #gravité
E = 1E11        #Module d'Young
nu = 0.25       #coef de poisson
rho_m = 3300    #densité du manteau
rho_q = 1000   #densité de la charge
rho_f = 1000    #densité du matériel de remplissage
Te = [5,10,15,20,25,30,35,40]      #épaisseur élastique
coef = 1	 #1	#-1E-2 #-0.28E-2	#Si besoin, coefficient de conversion pour le fichier d'entrée


psave = 1 # pas d'enregistrement dans le fichier de sortie (si psave = 40, un point enregistré tous les 40 points)


#
nb = 1
nb_exist = 0

RUN_GFLEX = True





###################################################################################################
###    Loop over tau and Te   -  No need to modify anything below                               ###
###################################################################################################


for f in np.arange(nb_exist,nb,1) :
	if nb == 1 :
		#DATA = PATH + PATHFILE + FILE + '.asc'
		DATA = DATA_ini
		NAME = NAME_ini
	else : 
		DATA = PATH + PATHFILE + 'random/' + FILE + '_' + str(f) + '.asc'
		NAME = NAME_ini + '_' + str(f)
		print('tirage = ' + str(f))
	
	
	if RUN_GFLEX :		
		for i in tau :
			print('tau =' + str(i))
			for j in Te :
				print('TE =' + str(j))
				NAME_model = NAME + '_Te' + str(j) + 'km_tau' + str(i) +'yr'
			
				if nb == 1 :
					os.makedirs(PATH+'output/'+NAME+'/'+NAME_model,exist_ok=True)
					INPUT = open(PATH+'output/'+NAME+'/'+NAME_model+'/input_'+NAME_model,'w+')
				
				else :
					os.makedirs(PATH+'output/RANDOM/'+NAME+'/'+NAME_model,exist_ok=True)
					INPUT = open(PATH+'output/RANDOM/'+NAME+'/'+NAME_model+'/input_'+NAME_model,'w+')
			
				INPUT.write('Model name = {} \nFault name = {} \nDip = {}     Dip north = {} \nFriction = {} \n'.format(NAME_model,NAME_fault,dip,dip_left,friction))
				INPUT.write('Tau = {} \nTfin = {} \ng = {} \nE = {} \nnu = {} \nrho_m = {} \nrho_fill = {} \nTe = {}'.format(tau,tfin,g,E,nu,rho_m,rho_f,Te))

				gFlex_param = [method,g,E,nu,rho_m,rho_f,j*1000,coef]
		
			
				#RUN#
				gflex(PATH,NAME,DATA,tfin,i,gFlex_param,systeme_coord,gps_point,GPS,grid_coord,ADDLOAD,uniform_grid,psave,rho_q,nb)


if NAME_fault != 'NONE' : 
	for i in tau :
		for j in Te :
		
			print('\n \n \n \n \n \n \n \n \n TE =' + str(j) + '\n ')
			
			NAME_model = NAME + '_Te' + str(j) + 'km_tau' + str(i) +'yr'
			os.makedirs(PATH+'output/'+NAME+'/'+NAME_model+'/'+NAME_fault,exist_ok=True)
			
			for d in dip :
				for dn in dip_left :
					for f in friction :
	
						fault = [NAME_fault,d,dn,f] #latitude, longitude, azimuth, pendage, friction
						FAULTS(PATH,NAME,i,j,systeme_coord,fault,tfin)
							
