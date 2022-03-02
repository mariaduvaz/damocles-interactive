#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 14:58:36 2021

@author: maria
"""


import os
import damocleslib as model
import numpy as np
import matplotlib.pyplot as plt
import sys
import fileinput
import FUNCTIONS as fn 
from mpl_toolkits import mplot3d
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button,CheckButtons,TextBox
from matplotlib.artist import Artist


#get directory root
path = os.path.dirname(os.path.realpath(__file__))

#specify redshift of object here
z_red=0.034

#read in datafile, first column as wavelength second column as flux values
obsfile ="iPTF14hls_2016-11-08_14-31-56_FTN_FLOYDS-N_iPTF.ascii"
obswav,obsflux= fn.datafile_2_array(obsfile,isint=False,zipped=True)
obsflux = [i-4.2e-17 for i in obsflux]

#trim_lims determine the limits of velocity space you want to view your line profile at
trim_lims = (-5000,7000)
#choose an area of continuum in the spectra where there is no emission features in velocity space. the standard dev of flux in this region gives observational uncertainty which is used in the chi sq calculation
bg_lims = (-50000,-25000)
snip_regs = () #(-315,470)



#obsvels,obsflux = fn.trim_wav_flux(obsvels,obsflux,trim_lim_1,trim_lim_2)
#for chi sq calc
#bg1 = 
#bg2 =

##########################################                            #############################################
##########################################  PARAMETERS WE'RE VARYING  #############################################
##########################################                            #############################################





#INITIAL PARAMETERS THAT MAKE OUR GRID OF GAS

v_max_init=4130 #maximum velocity at edge of shell
Rrat_init = 0.27 #radius of inner shell boundary divided by ratio of outer shell. Determines how 'thick' your gas shell is
rho_index_init=1.32   #density is proportional to radius^(-rho_index)
mdust_init=4.6e-05        #mass of dust in solar masses
grain_size_init=0.12    #size of dust grain in microns



#put in the wavelength of the spectral line transition you want to create a model of in NANOMETRES
#set this to true if you're modelling a doublet (two close together lines that have blended into each other like O2+ 4959,5007)
doublet= "false"
wavelength_peak_1= 656.3 #729.15
#set this to 0 if you don't have a doublet
wavelength_peak_2= 732.39




#calculate observational uncert on input spectrum using background regions provided

lam_o= (1+z_red)*(wavelength_peak_1*10.0)
obsvels = fn.convert_wav_to_vel(obswav,lam_o,wavelength_peak_1*10)


bg_vels,bg_flux = fn.trim_wav_flux(obsvels,obsflux,bg_lims[0],bg_lims[1])
#plt.plot(bg_vels,bg_flux)
#plt.plot(obsvels,obsflux)
#plt.show()

obs_err = np.std(bg_flux)
#write processed spec to line.out, which is used to do the chi calcn in damocles! 

obsvels,obsflux = fn.trim_wav_flux(obsvels,obsflux,trim_lims[0],trim_lims[1])
obsflux = fn.snip_spect(obsvels,obsflux,*snip_regs)
#plt.plot(obsvels,obsflux)
#plt.show()

#scale which models line profile peaks to roughly the same size, may need to be manually adjusted depending on SNR
obs_scale = np.amax(obsflux)*60

filey = open(path+"/input/line.in",'w')
filey.write(str(len(obsflux))+ " " + str(obs_err) + "\n")
for j in range(len(obsflux)):
    
                    filey.write(str(obsvels[j]) + ' ' + str(obsflux[j]) + "\n")

filey.close()



#age of supernova remnant at the time that observational data was taken, in days.
#this is important because the older the supernova is, the larger its radius would be as its had more time to expand. A larger supernova means you need more dust to get the same amount of light being absorbed

age_d = 778 #1170                                 #INSERT AGE IN DAYS HERE!!!
##no of grid cells in x,y,z direction. higher grid number means higher resolution. INCREASE THIS WHEN YOU HAVE A GOOD ESTIMATE OF WHAT YOUR PARAMETERS SHOULD BE
grid_divs = 20   
#no of photon packets w. which to run simulation. more packets = more time for each simulation to run and higher SNR model line profile
photno = 70000


###HERE WE'RE PUTTING IN THE PARAMETERS WE SET ABOVE INTO THE CODE, FEEDING THEM INTO THE INPUT FILES
#Read in input files as arrays
input_file = "input/input.in"
dust_file = "input/dust.in"
gas_file = "input/gas.in"
spec_file = "input/species.in"
inlines = fn.datafile_2_array(input_file,isint=False,zipped=False)
dustlines = fn.datafile_2_array(dust_file,isint=False,zipped=False)
gaslines = fn.datafile_2_array(gas_file,isint=False,zipped=False)
speclines = fn.datafile_2_array(spec_file,isint=False,zipped=False)



#Replace values in files in input.in fortran file with values defined in this script


fi = fileinput.FileInput(files=(input_file,dust_file,gas_file,spec_file),inplace=True)
for line in fi:
   if 'day' in line:
       line=fn.replace_str(age_d,0,line)
   if  'number of photons' in line:
       line=fn.replace_str(photno,0,line)
   if  'total flux of line to scale model' in line:
       line=fn.replace_str(obs_scale,0,line)
   if 'doublet?' in line:
       line=fn.replace_str(doublet,0,line)
   if 'Msun' in line:
       line=fn.replace_str(mdust_init,0,line)
   if "first doublet component" in line:
       line=fn.replace_str(wavelength_peak_1,0,line)
   if "second doublet component" in line:
       line=fn.replace_str(wavelength_peak_2,0,line)
   if  'dustData' in line:
       line=fn.replace_str(grain_size_init,3,line)
       line=fn.replace_str(grain_size_init,4,line)
  
   sys.stdout.write(line)  

fi.close()




#creating a gridfile of the supernova here
#where we have 4 1d arrays; a list of x,y,z and density points
#this allows us to plot what the model looks like 




fig2 = plt.figure(2,figsize=(8,8))
plt.xlabel("Velocity km/s",fontsize=20)
plt.ylabel("Brightness",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=20)

#button to clear line emission plotting window
resetax = plt.axes([0.25, 0.0, 0.65, 0.04])

button = Button(resetax, 'Clear')

#button to add in factor by which to rescale model line profile with
scaleax = plt.axes([0.5, 0.1, 0.2, 0.04])
scale_input = TextBox(scaleax, 'Scale model line profile by:', initial='1.0', color='.95', hovercolor='1', label_pad=0.01)

plt.subplots_adjust(left=0.25, bottom=0.25)




ax = fig2.add_subplot(111)


ax.set_xlim([trim_lims[0],trim_lims[1]])
#ax.set_ylim([0,3])

#set axis labels
fig = plt.figure(1,figsize=(10,10))
#fig = plt.subplot(1,2,2)
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.25)



x,y,z,d = fn.make_Grid(v_max_init,Rrat_init,rho_index_init,age_d,grid_divs)
g_s= np.amax(x)
#function contains axis orientation, limits and labels
def setax():
	ax.view_init(elev=30, azim=50)
	ax.set_xlabel('X axis (cm)')
	ax.set_ylabel('Y axis (cm)')
	ax.set_zlabel('Z axis (cm)')
	ax.set_xlim([-2*g_s,2*g_s])
	ax.set_ylim([-2*g_s,2*g_s])
	ax.set_zlim([-2*g_s,2*g_s])
	ax.set_title("Model of gas distribution in a Supernova")



setax()


#legend tells you params
#param_label = "vmax: "+str(v_max_init)+","+"Rrat: "+str(Rrat)+","+"rho_index: "+str(rho_index) dont need label as we have sliders
l = ax.scatter(x,y,z,c=d,cmap="nipy_spectral")
cbar = fig.colorbar(l)
cbar.set_label('density', rotation=270,size=18,labelpad=20)
cbar.ax.tick_params(labelsize=13)




ax_vmax = plt.axes([0.25, 0.0, 0.65, 0.02])
ax_r = plt.axes([0.25, 0.05, 0.65, 0.02])
ax_rho = plt.axes([0.25, 0.1, 0.65, 0.02])
ax_dm = plt.axes([0.25, 0.15, 0.65, 0.02])
ax_gs = plt.axes([0.25, 0.2, 0.65, 0.02])

s_vmax = Slider(ax_vmax, 'V$_{max}$', 1000, 15000,valinit=v_max_init)
s_vmax.label.set_size(15)
s_r = Slider(ax_r, 'R$_{in}$/R$_{out}$', 0.001, 1,valinit=Rrat_init)
s_r.label.set_size(15)
s_rho = Slider(ax_rho, r'$\beta$', 0.001, 3,valinit=rho_index_init)
s_rho.label.set_size(15)
s_dm = Slider(ax_dm, 'Dust mass (M$\odot$)', 0.0, 0.005,valfmt='%9.7f',valinit=mdust_init)
s_dm.label.set_size(15)
s_gs = Slider(ax_gs, 'Grain radius ($\mu$m)', 0.005, 0.5,valinit=grain_size_init)
s_gs.label.set_size(15)

amc_ax= plt.axes([0.15, 0.26, 0.08, 0.07])
#IF AMC BOX IS TICKED THEN 100% AMC USED, IF NOT THEN 100 SIL IS USED
amc_button=CheckButtons(amc_ax, ['AmC?'],actives=[False])


clump_ax= plt.axes([0.21, 0.26, 0.1, 0.07])
#IF AMC BOX IS TICKED THEN 100% AMC USED, IF NOT THEN 100 SIL IS USED
clump_button=CheckButtons(clump_ax, ['Clumped?'],actives=[False])




def update(val):
  #plt.cla()
  

  newv = s_vmax.val
  newr = s_r.val
  newrho = s_rho.val
  newdm = s_dm.val
  newgs	= s_gs.val
  
  
  #using the updated parameters, remake the grid points and save this grid to dust_grid.in file 
  x,y,z,d = fn.make_Grid(newv,newr,newrho,age_d,grid_divs)
  
  
  
  plt.figure(1)
  ax.clear()
  
  
  
  setax()
  l= ax.scatter(x,y,z,c=d,cmap="nipy_spectral") # update the plot
  
  #writing values we've updated to shell version of damocles as something is going wrong in the code w. arbitrary setting
  fi2 = fileinput.FileInput(files=(dust_file),inplace=True)
  for lineo in fi2:	
    if 'max dust velocity' in lineo:
     lineo=fn.replace_str(newv,0,lineo)
    if 'Rin/Rout' in lineo:
     lineo=fn.replace_str(newr,0,lineo)
    if 'rho~r^-q' in lineo:
     lineo=fn.replace_str(newrho,0,lineo)
    if 'Total dust mass' in lineo:
     lineo=fn.replace_str(newdm,0,lineo)
    sys.stdout.write(lineo)

  fi2.close()

  fi3 = fileinput.FileInput(files=(spec_file),inplace=True)
  for lineo in fi3:	
       if  'dustData' in lineo:
          lineo=fn.replace_str(newgs,3,lineo)
          lineo=fn.replace_str(newgs,4,lineo)
       sys.stdout.write(lineo)  

  fi3.close()

  mod = model.run_damocles_wrap()    

  outfile = path + "/output/integrated_line_profile_binned.out"
  vel,flux,flux_e= fn.datafile_2_array(outfile,isint=False,zipped=True)
  
  ####read in tau value from model_properties.out file
  mod_prop = path + "/output/model_properties.out"
  props = fn.datafile_2_array(mod_prop,isint=False,zipped=False)

   #scaling flux by additional factor provided in Fig 2
  flux = [i * float(scale_input.text) for i in flux]
  
  
  def reset(event):
    plt.cla()
    plt.tick_params(axis='both', which='major',labelsize=20)
	
    plt.xlabel("Velocity km/s",fontsize=20)
    plt.ylabel("Brightness",fontsize=20)
  button.on_clicked(reset)

  
  def set_amc(event):
      
      fi3 = fileinput.FileInput(files=(spec_file),inplace=True)
      if amc_button.get_status() == [True]:
       
          
          for lineo in fi3:	
             
            if  'dustData' in lineo:
                 
                 lineo=fn.replace_str('\'dustData/amC-zb1.nk\'',1,lineo)
            sys.stdout.write(lineo)     
            
      else:
         for lineo in fi3:	
             if  'dustData' in lineo:
                 lineo=fn.replace_str('\'dustData/sil-dlee.nk\'',1,lineo)
         
             sys.stdout.write(lineo)
      fi3.close()
  
      #get status of button
  amc_button.on_clicked(set_amc)
  
  
  def set_clump(event):
      
      fi2 = fileinput.FileInput(files=(dust_file),inplace=True)
      if clump_button.get_status() == [True]:
       
          
          for lineo in fi2:	
             
            if  'fraction' in lineo:
                 
                 lineo=fn.replace_str('1.0',0,lineo)
            sys.stdout.write(lineo)     
            
      else:
         for lineo in fi2:	
             if  'fraction' in lineo:
                 lineo=fn.replace_str('0.0',0,lineo)
         
             sys.stdout.write(lineo)
      fi2.close()
  
  clump_button.on_clicked(set_clump)



  obs_file = path +"/" + obsfile
  obswav,obsflux= fn.datafile_2_array(obs_file,isint=False,zipped=True)
    

  
  
  lam_o= (1.0+z_red)*wavelength_peak_1*10.0
  obsvels = fn.convert_wav_to_vel(obswav,lam_o,wavelength_peak_1*10.0)


  

  obsvels,obsflux = fn.trim_wav_flux(obsvels,obsflux,trim_lims[0],trim_lims[1])
  obsflux = [i-4.2e-17 for i in obsflux]
  #obsflux = fn.snip_spect(obsvels,obsflux,-188,362,601,1286)
  obsflux = fn.snip_spect(obsvels,obsflux,*snip_regs)
  
  
  #extracting values from model_properties.out file by word matching. We then present the reduced chi sq value, as this is the most useful for determining whether a model is a good fit to the data
  for i in props:
       
        if 'optical' in i and 'depth' in i and 'wavelength' in i and 'cell' not in i and '(absn)' not in i:
           tau = i[-1]
        if 'chi' in i:   
             
           chi_sqq = float(i[-1])/len(flux)
  
  
  
  plt.figure(2)
  
  plt.plot(vel,flux)
  plt.plot(obsvels,obsflux)
      
  plt.tick_params(axis='both', which='major',labelsize=20)
	
  plt.xlabel("Velocity km/s",fontsize=20)
  plt.ylabel("Brightness",fontsize=20)
  
  frame1=plt.text(5000,np.amax(obsflux),"tau=" + str(tau),fontsize=20)
  frame2=plt.text(3500,np.amax(obsflux)-np.amax(obsflux)/10,"chi sq=" + str(chi_sqq),fontsize=20)        
  
  
  fig2.canvas.draw()
  Artist.remove(frame1)
  Artist.remove(frame2)	
  
 
  
  

s_vmax.on_changed(update)
s_r.on_changed(update)
s_rho.on_changed(update)
s_dm.on_changed(update)
s_gs.on_changed(update)




plt.show()
#'''
