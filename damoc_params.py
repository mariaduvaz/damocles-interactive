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
import ipywidgets as widgets


#################################################################################################################################################
##########################################  PARAMETERS TO BE SPECIFIED BEFORE RUNNING CODE          #############################################
##########################################  These parameters are parsed to the damocles input files #############################################
#################################################################################################################################################


#specify redshift of object here
z_red=0.034



#trim_lims determine the limits of wavelength space where your line profile is at. 
trim_lims = (6650,6950)
#choose an area of continuum in the spectra where there is no emission features in velocity space. the standard dev of flux in this region gives observational uncertainty which is used in the chi sq calculation
bg_lims = (5809,6300)
#specify a list here of regions in the line profile where there could be contaminating features also in wavelength space which you could remove
snip_regs = () 


#put in the wavelength of the spectral line transition you want to create a model of in NANOMETRES
#set this to true if you're modelling a doublet (two close together lines that have blended into each other like O2+ 4959,5007)
is_doublet= "false"
wavelength_peak_1= 656.3 
wavelength_peak_2= 732.39 
doublet_ratio = 3.13

#These are our default model parameters values for the sliders
-5000,7000
v_max_init=4130  #maximum velocity at edge of shell
Rrat_init = 0.27 #radius of inner shell boundary divided by ratio of outer shell. Determines how 'thick' your gas shell is
rho_index_init=1.32   #density is proportional to radius^(-rho_index)
mdust_init=4.6e-05        #mass of dust in solar masses
grain_size_init=0.12    #size of dust grain in microns



#age of supernova remnant at the time that observational data was taken, in days.
age_d = 778  
##no of grid cells in x,y,z direction used to make the spherical shell model in damocles
grid_divs = 20   
#no of photon packets w. which to run simulation. more packets = more time for each simulation to run and higher SNR model line profile
phot_no = 70000



##########################################
########  READING IN DATAFILES    ########
##########################################


path = os.path.dirname(os.path.realpath(__file__))

obsfile ="iPTF14hls_2016-11-08_14-31-56_FTN_FLOYDS-N_iPTF.ascii"
input_file = "input/input.in"
dust_file = "input/dust.in"
gas_file = "input/gas.in"
spec_file = "input/species.in"

obswav,obsflux= fn.datafile_2_array(obsfile,isint=False,zipped=True)
inlines = fn.datafile_2_array(input_file,isint=False,zipped=False)
dustlines = fn.datafile_2_array(dust_file,isint=False,zipped=False)
gaslines = fn.datafile_2_array(gas_file,isint=False,zipped=False)
speclines = fn.datafile_2_array(spec_file,isint=False,zipped=False)


###########################################################################################
######################## initial processing of observed spectrum  #########################
###########################################################################################

#line needs to be removed, this is just for particular example spec of ipt14hls as continuum needs to be moved down
obsflux = [i-4.2e-17 for i in obsflux]


#calculate observational uncert on input spectrum using background regions provided
bg_vels,bg_flux = fn.trim_wav_flux(obswav,obsflux,bg_lims[0],bg_lims[1])
obs_err = np.std(bg_flux)


####trim spectra to line you want
obswav,obsflux = fn.trim_wav_flux(obswav,obsflux,trim_lims[0],trim_lims[1])
#snip out narrow line or contaminating region
obsflux = fn.snip_spect(obswav,obsflux,*snip_regs)

#scale which models line profile peaks to roughly the same size, may need to be manually adjusted depending on SNR
#####NEED TO CHECK WHETHER THIS VALUE CHANGES GIVEN CHANGING PHOTON NO 
obs_scale = np.amax(obsflux)*60


#convert from wavelength space to velocity space and correct for redshift of supernova's host galaxy
lam_o= (1+z_red)*(wavelength_peak_1*10.0)
obsvels = fn.convert_wav_to_vel(obswav,lam_o,wavelength_peak_1*10)



#write observed line to line.out file, which is used to do the chi sq calcn in damocles! 
filey = open(path+"/input/line.in",'w')
filey.write(str(len(obsflux))+ " " + str(obs_err) + "\n")
for j in range(len(obsflux)):    
                    filey.write(str(obsvels[j]) + ' ' + str(obsflux[j]) + "\n")
filey.close()





#Replace values in files in input.in fortran file with initial slider values defined in script

fi = fileinput.FileInput(files=(input_file,dust_file,gas_file,spec_file),inplace=True)
for line in fi:
   if 'day' in line:
       line=fn.replace_str(age_d,0,line)
   if  'number of photons' in line:
       line=fn.replace_str(phot_no,0,line)
   if  'total flux of line to scale model' in line:
       line=fn.replace_str(obs_scale,0,line)
   if 'doublet?' in line:
       line=fn.replace_str(is_doublet,0,line)
   if 'Msun' in line:
       line=fn.replace_str(mdust_init,0,line)
   if "first doublet component" in line:
       line=fn.replace_str(wavelength_peak_1,0,line)
   if "second doublet component" in line:
       line=fn.replace_str(wavelength_peak_2,0,line)
   #unless otherwise specified via the interactive button, the default dust distribition is smooth
   if  'fraction in clumps' in line:                
                 line=fn.replace_str('0.0',0,line)
   if  'dustData' in line:
       line=fn.replace_str(grain_size_init,3,line)
       line=fn.replace_str(grain_size_init,4,line)
  
   sys.stdout.write(line)  

fi.close()




########################################################################################################################################################
##########setting environment parameters and creating widgets for the observed spectral line and damocles model comparison plotting window #############
########################################################################################################################################################

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

trim_lims_vel = fn.convert_wav_to_vel(trim_lims,lam_o,wavelength_peak_1*10.0)
ax.set_xlim([trim_lims_vel[0],trim_lims_vel[1]])

plt.plot(obsvels,obsflux)

########################################################################################################################################################
##########setting environment parameters and creating widgets for the 3D grid and sliders plotting window #############
########################################################################################################################################################



fig = plt.figure(1,figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.25)


x,y,z,d = fn.make_Grid(v_max_init,Rrat_init,rho_index_init,age_d,grid_divs)
grid_lim= np.amax(x)
fn.setax(ax,grid_lim)


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

#IF AMC BOX IS TICKED THEN 100% AMC USED, IF NOT THEN 100 SIL IS USED
amc_ax= plt.axes([0.15, 0.26, 0.08, 0.07])
amc_button=CheckButtons(amc_ax, ['AmC?'],actives=[False])

#IF CLUMP BOX IS TICKED THEN CLUMPED DISTRIBUTION IS USED, IF NOT THEN A SMOOTH ONE IS
clump_ax= plt.axes([0.21, 0.26, 0.1, 0.07])
clump_button=CheckButtons(clump_ax, ['Clumped?'],actives=[False])






###EVERY TIME A SLIDER OR BUTTON IS INTERACTED WITH, THIS FUNCTION IS CALLED####

def update(val):
  #plt.cla()
  
  
    
  def reset(event):
    plt.cla()
    plt.tick_params(axis='both', which='major',labelsize=20)
    
    obs_file = path +"/" + obsfile
    obswav,obsflux= fn.datafile_2_array(obs_file,isint=False,zipped=True)
  
    obswav,obsflux = fn.trim_wav_flux(obswav,obsflux,trim_lims[0],trim_lims[1])
    obsflux = [i-4.2e-17 for i in obsflux]
    obsflux = fn.snip_spect(obswav,obsflux,*snip_regs)
    lam_o= (1.0+z_red)*wavelength_peak_1*10.0
    obsvels = fn.convert_wav_to_vel(obswav,lam_o,wavelength_peak_1*10.0)

    plt.plot(obsvels,obsflux)
    plt.xlabel("Velocity km/s",fontsize=20)
    plt.ylabel("Brightness",fontsize=20)

  
  def set_amc(event):
      
      fi3 = fileinput.FileInput(files=(spec_file),inplace=True)
      if amc_button.get_status() == [True]: 
          for line in fi3:	           
            if  'dustData' in line:                 
                 line=fn.replace_str('\'dustData/amC-zb1.nk\'',1,line)
            sys.stdout.write(line)     
            
      else:
         for line in fi3:	
             if  'dustData' in line:
                 line=fn.replace_str('\'dustData/sil-dlee.nk\'',1,line)        
             sys.stdout.write(line)
      fi3.close()
  

  def refit_scaled(event,arg1,arg2,arg3,arg4,arg5,arg6):
               

        arg1 = [i * float(scale_input.text) for i in arg1]

        chi_sq_rescale = fn.chi_sq(arg2,arg1,arg3,arg4)
        print("chi sq: ",chi_sq_rescale)

        plt.plot(arg5,arg1)
        
   #     frame3=plt.text(3500,np.amax(obsflux)-np.amax(obsflux)/10,"chi sq=" + str(chi_sq_rescale),fontsize=20)  
        
        arg6.canvas.draw()
     #   Artist.remove(frame3)	
  
  def set_clump(event): 
      fi2 = fileinput.FileInput(files=(dust_file),inplace=True)
      if clump_button.get_status() == [True]:        
          for line in fi2:	             
            if  'fraction' in line:                
                 line=fn.replace_str('1.0',0,line)
            sys.stdout.write(line)     
            
      else:
         for line in fi2:	
             if  'fraction' in line:
                 line=fn.replace_str('0.0',0,line)         
             sys.stdout.write(line)
      fi2.close()
  
  
  #if any of the buttons are pressed, this calls the inset functions
  amc_button.on_clicked(set_amc)
  clump_button.on_clicked(set_clump)

  button.on_clicked(reset)
 
  newv = s_vmax.val
  newr = s_r.val
  newrho = s_rho.val
  newdm = s_dm.val
  newgs	= s_gs.val
  
  
  #remaking the 3D grit plot with new slider values
  plt.figure(1)
  ax.clear()
   
  x,y,z,d = fn.make_Grid(newv,newr,newrho,age_d,grid_divs)
  grid_lim= np.amax(x)
  fn.setax(ax,grid_lim)
  
  l= ax.scatter(x,y,z,c=d,cmap="nipy_spectral") # update the plot
  
  
  #writing values we've updated to shell version of damocles 
  fi2 = fileinput.FileInput(files=(dust_file,spec_file),inplace=True)
  for line in fi2:	
    if 'max dust velocity' in line:
     line=fn.replace_str(newv,0,line)
    if 'Rin/Rout' in line:
     line=fn.replace_str(newr,0,line)
    if 'rho~r^-q' in line:
     line=fn.replace_str(newrho,0,line)
    if 'Total dust mass' in line:
     line=fn.replace_str(newdm,0,line)
    if  'dustData' in line:
          line=fn.replace_str(newgs,3,line)
          line=fn.replace_str(newgs,4,line)     
    sys.stdout.write(line)
  fi2.close()



  #This line is calling damocles to run after all input parameter files have been updated 
  model.run_damocles_wrap()    

  #reading model output file after damocles has run
  outfile = path + "/output/integrated_line_profile_binned.out"
  vel,flux,flux_e= fn.datafile_2_array(outfile,isint=False,zipped=True)
  
  #file that tells you tau and chi sq
  mod_prop = path + "/output/model_properties.out"
  props = fn.datafile_2_array(mod_prop,isint=False,zipped=False)

  #run chi sq here 
  #scaling model flux by additional factor provided in Fig 2
  
  #have to use the lambda function to pass the refit_scaled def more than one arg, as by default it only accepts one, the event handler


  
  
  #extracting values from model_properties.out file by word matching. We then present the reduced chi sq value, as this is the most useful for determining whether a model is a good fit to the data
  for i in props:      
        if 'optical' in i and 'depth' in i and 'wavelength' in i and 'cell' not in i and '(absn)' not in i:
           tau = i[-1]


  
  scale_input.on_submit(lambda event, args=(flux,obsflux,obs_err,flux_e,vel,fig2): refit_scaled(event,*args))
  flux = [i * float(scale_input.text) for i in flux]
  
  
  chi_sqq = fn.chi_sq(obsflux,flux,obs_err,flux_e) 
        
  
  
  ###plotting model against observed line profile panel
  
  plt.figure(2)
  
  plt.plot(vel,flux)
  
      
  plt.tick_params(axis='both', which='major',labelsize=20)
	
  plt.xlabel("Velocity km/s",fontsize=20)
  plt.ylabel("Brightness",fontsize=20)
  
  frame1=plt.text(5000,np.amax(obsflux),"tau=" + str(tau),fontsize=20)
  frame2=plt.text(3500,np.amax(obsflux)-np.amax(obsflux)/10,"chi sq=" + str(chi_sqq),fontsize=20)        
  

  
  #need to reconfigure this so every time a button is activated the text is not removed
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
