#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 14:58:36 2021

@author: maria
"""
import tkinter as tk
import tkinter.font as TkFont

import os
import damocleslib as model
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import sys
import fileinput
from FUNCTIONS import *
from mpl_toolkits import mplot3d
import numpy as np
from matplotlib.figure import Figure
###


class Slider():
    def __init__(self,frame_a,frame_b,frame_c,frame_d,init_val,range_vals,var_name,step_size):
        
       
 
        self.var_name = var_name
        self.sliderfont = TkFont.Font(family='bitstream charter', size=14)
        self.lab = tk.Label(frame_a,text=self.var_name,font=self.sliderfont)
        self.slider = tk.Scale(frame_a,from_=range_vals[0],to=range_vals[1],orient='horizontal',resolution=step_size,width=17,length=900,font=self.sliderfont)      
        self.slider.set(init_val)
        self.slider.bind("<ButtonRelease-1>", self.slider_command)      
        self.lab.pack(fill='x', padx=1)      
        self.slider.pack()
        self.frame_a=frame_a
        self.frame_b=frame_b
        self.frame_c=frame_c
        self.frame_d=frame_d

    def slider_command(self,event):
       
        slider_vals = []
        for child_widget in self.frame_a.winfo_children():
	          if child_widget.winfo_class() == 'Scale':
    	            slider_vals.append(child_widget.get())
        
        #de-log the dust mass and grain size values
  
        slider_vals[3] = 10**(slider_vals[3])
        slider_vals[4] = 10**(slider_vals[4])
        
        #pass all slider values upon a slider change 
        GasGrid.update_gasgrid(self.frame_b,slider_vals)
        DamoclesInput.update_damocfile_input(slider_vals)
        #run damocles after sliders have been changed and input files have been changed with slider values
        model.run_damocles_wrap() 
        Plotting_window.plot_model(self.frame_c)
        Plotting_window.update_tautextbox(self.frame_d)
        Plotting_window.update_chitextbox(self.frame_d)

        
    
        
class DamoclesInput():
    
    "This class contains all functions and variables which pass input parameters to the DAMOCLES code"
    
    def __init__(self):
        self.spec_file = "input/species.in"
        self.dust_file = "input/dust.in"
        self.buttonfont = TkFont.Font(family='bitstream charter', size=20)
        
        self.obswav_init,self.obsflux_init= datafile_2_array(obsfile,isint=False,zipped=True)
        self.obswav,self.obsflux = trim_wav_flux(self.obswav_init,self.obsflux_init,trim_lims[0],trim_lims[1])
        self.obsflux = snip_spect(self.obswav,self.obsflux,*snip_regs)
        #self.obsflux = [i-0.4e-16 for i in self.obsflux]
        self.obsvels = convert_wav_to_vel(self.obswav,(1+G_red)*(wavelength_peak_1*10.0),wavelength_peak_1*10)

  
        self.obs_err = self.get_obserr()
        self.write_obsfile_out()
        
        
    def get_obserr(self):
            #calculate observational uncert on input spectrum using background regions provided
                bg_vels,bg_flux = trim_wav_flux(self.obswav_init,self.obsflux_init,bg_lims[0],bg_lims[1])
                return np.std(bg_flux)
    
    def write_obsfile_out(self):
    #write observed line to line.out file, which is used to do rebin the damocles model 
        filey = open(path+"/input/line.in",'w')
        filey.write(str(len(self.obsflux))+ " " + str(self.obs_err) + "\n")
        for j in range(len(self.obsflux)):    
                        filey.write(str(self.obsvels[j]) + ' ' + str(self.obsflux[j]) + "\n")
        filey.close()


    
    def make_clump_button(self,frame):
        
        confirmation = tk.BooleanVar()
        clump_button = tk.Checkbutton(frame, text="Clump?",variable=confirmation,onvalue='True',offvalue='False',command=lambda: self.is_Clump(confirmation),bg='red',
                                      activebackground='red',font=self.buttonfont,borderwidth=5,indicatoron=0)
        clump_button.pack(fill='x')
            

           
    def is_Clump(self,conf): 
      fi2 = fileinput.FileInput(files=(self.dust_file),inplace=True)
      
      if conf.get() == True:        
          for line in fi2:	             
            if  'fraction' in line:                
                 line=replace_str('1.0',0,line)
            sys.stdout.write(line)     
            
      else:
         for line in fi2:	
             if  'fraction' in line:
                 line=replace_str('0.0',0,line)         
             sys.stdout.write(line)
      fi2.close()
    
    
    def update_damocfile_input(self,params):
         #writing values we've updated via sliders (where params comes from) to shell version of damocles 
         fi2 = fileinput.FileInput(files=(self.dust_file,self.spec_file),inplace=True)
         for line in fi2:	
             if 'max dust velocity' in line:
                 line=replace_str(params[0],0,line)
             if 'Rin/Rout' in line:
                 line=replace_str(params[1],0,line)
             if 'rho~r^-q' in line:
                 line=replace_str(params[2],0,line)
             if 'Total dust mass' in line:
                 line=replace_str(params[3],0,line)
             if  'dustData' in line:
                 line=replace_str(params[4],3,line)
                 line=replace_str(params[4],4,line)  
             if  'amC' in line:
                  line=replace_str(params[5],2,line)
             if 'sil' in line:
                  line=replace_str(str(float(1-params[5])),2,line)
             sys.stdout.write(line)
         fi2.close()
         
    def initialise_damocfile_input(self):
        fi = fileinput.FileInput(files=(input_file,dust_file,gas_file,spec_file),inplace=True)
        for line in fi:
            if 'day' in line:
                line=replace_str(age_d,0,line)
            if  'number of photons' in line:
                line=replace_str(phot_no,0,line)
            if  'total flux of line to scale model' in line:
                line=replace_str(str(np.amax(self.obsflux)*60),0,line)
            if 'doublet?' in line:
                line=replace_str(is_doublet,0,line)
            if "first doublet component" in line:
                line=replace_str(wavelength_peak_1,0,line)
            if "second doublet component" in line:
                line=replace_str(wavelength_peak_2,0,line)
       #unless otherwise specified via the interactive button, the default dust distribition is smooth
            if  'fraction in clumps' in line:                
                     line=replace_str('0.0',0,line)
                    
      
            sys.stdout.write(line)  

        fi.close()   
        
        
        
     

class GasGrid():
    
    def __init__(self,frame):
        self.fig = Figure(figsize=(8.5, 5.3), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.figure_canvas= FigureCanvasTkAgg(self.fig,frame)
    
    
    
    def make_Grid(self,v_max,Rrat,rho_index,age,divno):
	#creating a gridfile of the supernova here
#where we have 4 1d arrays; a list of x,y,z and density points
#this allows us to plot what the model looks like 

        v_min = v_max * Rrat
        Rout = v_max * age * 8.64e-6 * 1e15
        Rin = v_min* age * 8.64e-6 * 1e15
      
        #grid divisions for a uniform grid
        grid_arr = np.linspace(-Rout,Rout,divno)

        #These values contain every point in a 40*40*40 grid from the limits of -Rout,Rout 
        Y,X,Z = np.meshgrid(grid_arr,grid_arr,grid_arr) 
        #turning these into 1d arrays
        Xf = X.flatten()
        Yf = Y.flatten()  
        Zf = Z.flatten()
        rads= [np.sqrt(Xf[i]**2 + Yf[i]**2 + Zf[i]**2) for i in range(len(Xf))]

        rho = np.zeros((len(Xf)))


        plotdens,plotx,ploty,plotz = [],[],[],[]

	    #looping through radius of every grid point.
	    #if rad is within Rout and Rin, set density using r^-(rho_index) law
        for j in range(len(rads)):
           if abs(rads[j]) <= Rout and abs(rads[j]) >= Rin:
              if Xf[j] > 0 and Yf[j] > 0 and Zf[j] > 0:
                 rho[j] = (abs(rads[j]))**(-rho_index)*1e20   #randomly rescaling the density to make it a reasonable number
              else:
	           
                rho[j] = (abs(rads[j]))**(-rho_index)*1e20
                plotdens.append(rho[j])
                plotx.append(Xf[j])
                ploty.append(Yf[j])
                plotz.append(Zf[j])

  
        return plotx,ploty,plotz,plotdens



    def setax(self,axis,g_s):
       axis.view_init(elev=30, azim=50)
       axis.set_xlabel('X axis (cm)')
       axis.set_ylabel('Y axis (cm)')
       axis.set_zlabel('Z axis (cm)')
       axis.set_title("Model of gas distribution in a Supernova")
       axis.set_xlim([-1.5*g_s,1.5*g_s])
       axis.set_ylim([-1.5*g_s,1.5*g_s])
       axis.set_zlim([-1.5*g_s,1.5*g_s])

    
    
    def initialise_grid_axis(self,frame):
        
        #plotting initial grid made by specified initial parameter values and setting up color-bar scale and axis 
        
        x,y,z,d = self.make_Grid(v_max_init,Rrat_init,rho_index_init,age_d,grid_divs)
        grid_lim= np.amax(x)
        self.setax(self.ax,grid_lim)

        l = self.ax.scatter(x,y,z,c=d,cmap="nipy_spectral")
        cbar = self.fig.colorbar(l)
        cbar.set_label('density', rotation=270,size=18,labelpad=20)
        cbar.ax.tick_params(labelsize=13)
        self.fig.tight_layout()
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill='x', expand=1)
       


    def update_gasgrid(self,frame,params):
        
         try: 
             self.figure_canvas.get_tk_widget().pack_forget()
         except AttributeError: 
            pass 
    
         x,y,z,d = self.make_Grid(params[0],params[1],params[2],age_d,grid_divs)
        
          
         l = self.ax.scatter(x,y,z,c=d,cmap="nipy_spectral")
               
         self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill='x', expand=1)
         self.figure_canvas.draw()
  
    

class Plotting_window(DamoclesInput):
     def __init__(self,frame_a,frame_b,frame_c):
        self.frame_a = frame_a
        self.frame_b = frame_b
        self.frame_c = frame_c
        
        self.DamoclesInput =DamoclesInput()
        self.buttonfont = TkFont.Font(family='bitstream charter', size=16)
        
        self.fig = Figure(figsize=(7.5, 7.7), dpi=100)
        self.figure_canvas= FigureCanvasTkAgg(self.fig,self.frame_a)
        self.toolbar = NavigationToolbar2Tk(self.figure_canvas, self.frame_a)
        self.toolbar.update()
        self.outfile = path + "/output/integrated_line_profile_binned.out"
        self.mod_prop_file = path + "/output/model_properties.out" 
        
        tk.Label(self.frame_b,text = 'Optical depth (\u03C4) ',font=self.buttonfont).pack(anchor='w',side=tk.LEFT)
        self.tau_text = tk.Text(self.frame_b,height=2,width=20,font=self.buttonfont)
        self.tau_text.pack(anchor='w',side=tk.LEFT)
        
        tk.Label(self.frame_b,text = '\u03A7^2: ',font=self.buttonfont).pack(side=tk.LEFT)
        self.chi_text = tk.Text(self.frame_b,height=2,width=20,font=self.buttonfont)
        self.chi_text.pack(side=tk.LEFT)
       
        
        
         
     def initialise_plotwindow(self,frame):
           
          print("PLOT WINDOW",)
          trim_lims_vel = convert_wav_to_vel(trim_lims,(1+G_red)*(wavelength_peak_1*10.0),wavelength_peak_1*10.0)
          
          self.ax = self.fig.add_subplot(111)
          self.ax.axes.set_xlabel("Velocity (km/s)",fontsize=20)
          self.ax.axes.set_ylabel("Flux ($ergs$  $cm^{-2}$  $s^{-1}$  $\AA^{-1}$)",fontsize=20)
          self.ax.tick_params(axis='both', which='major',labelsize=20)
          self.ax.set_xlim([trim_lims_vel[0],trim_lims_vel[1]])
          self.ax.plot(self.DamoclesInput.obsvels,self.DamoclesInput.obsflux) 
          self.fig.tight_layout()
          self.figure_canvas.get_tk_widget().pack(fill='x', expand=1)
      
    
     def plot_model(self,frame):
         
          modvel,modflux,modflux_e = datafile_2_array(self.outfile,isint=False,zipped=True)
          #modflux = convolve_spectra(res_kms, modvel,modflux)
          modflux= [i * np.amax(DamoclesInput.obsflux)/np.amax(modflux) for i in modflux]
          
          #self.ax.text(np.amax(modvel)/2,Line.text_high1,Line.line_id_lab,fontsize=22)
          #self.ax.text(np.amax(modvel)/2,Line.text_high2,Line.dust_type_lab,fontsize=22)           
          #self.ax.text(np.amax(modvel)/2,Line.text_high1+(Line.text_high1-Line.text_high2),"iPTF14hls d"+str(age_d),fontsize=22)
          self.ax.plot(modvel,modflux) 
          self.figure_canvas.draw()
            
     def make_reset_button(self,frame):
         
        reset_button = tk.Button(frame, text="Reset",command=lambda: self.clear_pane(frame),bg='orange',
                                 activebackground='grey',font=self.buttonfont,borderwidth=5)
        reset_button.pack(fill='x',side=tk.BOTTOM,anchor='s')
     
     def make_model_scalebox(self,frame):
            labelText=tk.StringVar()
            labelText.set("Scale model by:")
            labelDir= tk.Label(frame, textvariable=labelText, height=3,font=self.buttonfont)
            labelDir.pack(fill='x',side=tk.LEFT,padx='20',pady='10')
         
            scale_var = tk.StringVar(value="1.0")
            
            def scalebox_command(var):
                sf = float(scale_var.get())
                modvel,modflux,modflux_e = datafile_2_array(self.outfile,isint=False,zipped=True)
                modflux= [i * np.amax(DamoclesInput.obsflux)/np.amax(modflux) * sf for i in modflux]
               # modflux = convolve_spectra(res_kms, modvel,modflux)
                chi = chi_sq(DamoclesInput.obsflux,modflux,DamoclesInput.obs_err,modflux_e) 
                chi = round(chi,2)
                self.chi_text.delete(1.0,6.0)
                self.chi_text.insert(tk.END,str(chi))
                
                self.ax.plot(modvel,modflux) 
                self.figure_canvas.draw()
            
            scale_var_entry = tk.Entry(self.frame_c, textvariable=scale_var)
            scale_var_entry.bind("<Return>", scalebox_command)   
            scale_var_entry.pack(fill = 'x', side=tk.LEFT, padx='20',pady='10')
  

     def clear_pane(self,frame):
         self.fig.clear() #clear your figure
         self.initialise_plotwindow(frame)
         
     def update_tautextbox(self,frame):
         props = datafile_2_array(self.mod_prop_file,isint=False,zipped=False)   
         for i in props:      
             if 'optical' in i and 'depth' in i and 'wavelength' in i and 'cell' not in i and '(absn)' not in i:
                 tau = i[-1]
         
         self.tau_text.delete(1.0,4.0)
         self.tau_text.insert(tk.END,str(tau))

         
     def update_chitextbox(self,frame):
         
         modvel,modflux,modflux_e= datafile_2_array(self.outfile,isint=False,zipped=True)
         #modflux = convolve_spectra(res_kms, modvel,modflux)
         chi = chi_sq(DamoclesInput.obsflux,modflux,DamoclesInput.obs_err,modflux_e) 
         chi = round(chi,2)
         self.chi_text.delete(1.0,6.0)
         self.chi_text.insert(tk.END,str(chi))
  
        
  
class InputWindow(tk.Tk):
    def __init__(self):
       
       super().__init__()
       self.geometry('500x800')
       
       self.z_var = tk.DoubleVar(value=0.0)
       self.SN_name = tk.StringVar()
       self.Line_name = tk.StringVar()
       self.bg_lims = tk.StringVar()
       self.trim_lims = tk.StringVar()
       self.snip_reg = tk.StringVar()
       self.is_doublet= tk.StringVar(value="false")
       self.wavelength_peak_1= tk.DoubleVar(value=656.3) 
       self.wavelength_peak_2= tk.DoubleVar(value=732.39) 
       self.doublet_ratio = tk.DoubleVar(value=3.13)
       self.resolution= tk.DoubleVar(value=10.0)
       #age of supernova remnant at the time that observational data was taken, in days.
       self.age_d = tk.DoubleVar(value=778)
       ##no of grid cells in x,y,z direction used to make the spherical shell model in damocles
       self.grid_divs = tk.IntVar(value=20)   
       #no of photon packets w. which to run simulation. more packets = more time for each simulation to run and higher SNR model line profile
       self.phot_no = tk.IntVar(value=30000)
       
       

       self.start_vars = {"Host Redshift:": [self.z_var,None]}#,"SN name": ['iptf14hls',None]}
       
       for i in list(self.start_vars.keys()):
           #print(i,self.start_vars.get(i)[0])
           self.make_label_entry(i,self.start_vars.get(i)[0])
       
        
      
      # Button to be clicked which opens up modelling app when fields are complete
       tk.Button(self,
                 text='Open modelling app',
                 command=self.open_window).pack(side=tk.BOTTOM,expand=True)
       
    def make_label_entry(self,labelname,variablename):
       
       labelText=tk.StringVar()
       labelText.set(labelname)
       labelDir= tk.Label(self, textvariable=labelText, height=3)
       labelDir.pack(side=tk.LEFT)
      
       z_entry = tk.Entry(self, textvariable=variablename)
       z_entry.pack(fill = 'x', side=tk.LEFT, padx='20',pady='10')
         

     
    
    def open_window(self,event=None):
      for i in list(self.start_vars.keys()):
          self.start_vars.get(i)[1] = self.start_vars.get(i)[0].get()
      print(self.start_vars)
      window = App(self)
      window.grab_set()
      
    
    
class App(tk.Toplevel,GasGrid,Slider,Plotting_window):
    def __init__(self,parent):
      
        super().__init__(parent)
        
        print("IN APP",InputWindow.start_vars)
        self.geometry("2000x1500")
        self['bg'] = 'blue'
 
        
        frame_2 = tk.Frame(self)
        frame_1 = tk.Frame(self)
        frame_3 = tk.Frame(self)
        frame_4 = tk.Frame(self)
        frame_5 = tk.Frame(self)
        
        frame_2.pack()
        frame_1.pack()
        frame_3.pack()
        frame_4.pack()
        frame_5.pack()
 
        
        
        self.DamoclesInput = DamoclesInput()
        self.GasGrid = GasGrid(frame_2)
        
        self.Plotting_window = Plotting_window(frame_3,frame_4,frame_5)
       
        
        self.DamoclesInput.initialise_damocfile_input()
        self.DamoclesInput.make_clump_button(frame_1)
        
        #create sliders for parameters that are changed by user
        v_slider = Slider(frame_1,frame_2,frame_3,frame_4,v_max_init,(1000, 15000),"Vmax (km/s)",1)
        r_slider = Slider(frame_1,frame_2,frame_3,frame_4,Rrat_init,(0.01, 1),"Rin/Rout",0.0005)
        rho_slider = Slider(frame_1,frame_2,frame_3,frame_4,rho_index_init,(-6, 6),"Density index (\u03B2)",0.01)
        md_slider = Slider(frame_1,frame_2,frame_3,frame_4,mdust_init,(-9, 0.2),"log(Dust mass (M\u2609))",0.001)
        gs_slider = Slider(frame_1,frame_2,frame_3,frame_4,grain_size_init,(-2.3, 0.5),'log(Grain radius (\u03BCm))',0.001)
        amc_frac_slider = Slider(frame_1,frame_2,frame_3,frame_4,0.0,(0.0,1.0),'AmC Fraction',0.01)
        
        
        self.GasGrid.initialise_grid_axis(frame_2)
       
        self.Plotting_window.initialise_plotwindow(frame_3)
        self.Plotting_window.make_reset_button(frame_3)
        self.Plotting_window.make_model_scalebox(frame_5)
        
        frame_1.place(x=950,y=515)
        frame_2.place(x=980,y=0) 
        frame_3.place(x=100,y=10)
        frame_4.place(x=105,y=865)
        frame_5.place(x=235,y=920)
        
    
#################################################################################################################################################
##########################################  PARAMETERS TO BE SPECIFIED BEFORE RUNNING CODE          #############################################
##########################################  These parameters are parsed to the damocles input files #############################################
#################################################################################################################################################

#Line = Model_spect('ha',"Ha 6563$\AA$","100% sil clump","iPTF14hls_2016-11-08_14-31-56_FTN_FLOYDS-N_iPTF.ascii","models/ha-sil-clump-d778",6563,-4.2e-17,1.75,1.9,752,5.6e-16,5.3e-16,-0.5e-17,6.4e-16,(),(-8000,8000))

#specify redshift of object here
G_red=0.034
SN_name = 'iPTF14hls'
Line_name = 'Ha'

#trim_lims determine the limits of wavelength space where your line profile is at. 
trim_lims = (6650,7050)
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

#specify resoln
resolution= 10.0
res_kms = resolution/(wavelength_peak_1*10) * 299792 
print(res_kms)
#These are our default model parameters values for the sliders

v_max_init = 4130  #maximum velocity at edge of shell
Rrat_init = 0.31 #radius of inner shell boundary divided by ratio of outer shell. Determines how 'thick' your gas shell is
rho_index_init=1.330   #density is proportional to radius^(-rho_index)
mdust_init=-3.67       #mass of dust in solar masses
grain_size_init=-1.235    #size of dust grain in microns



#age of supernova remnant at the time that observational data was taken, in days.
age_d = 778  
##no of grid cells in x,y,z direction used to make the spherical shell model in damocles
grid_divs = 20   
#no of photon packets w. which to run simulation. more packets = more time for each simulation to run and higher SNR model line profile
phot_no = 30000



##########################################
########  READING IN DATAFILES    ########
##########################################


path = os.path.dirname(os.path.realpath(__file__))

obsfile ="iPTF14hls_2016-11-08_14-31-56_FTN_FLOYDS-N_iPTF_contsub.ascii"
input_file = "input/input.in"
dust_file = "input/dust.in"
gas_file = "input/gas.in"
spec_file = "input/species.in"


inlines = datafile_2_array(input_file,isint=False,zipped=False)
dustlines = datafile_2_array(dust_file,isint=False,zipped=False)
gaslines = datafile_2_array(gas_file,isint=False,zipped=False)
speclines = datafile_2_array(spec_file,isint=False,zipped=False)


###########################################################################################
######################## initial processing of observed spectrum  #########################
###########################################################################################








if __name__ == '__main__':
  
  InputWindow = InputWindow()
  InputWindow.mainloop()
  
  
'''
  rootwindow = tk.Tk()
  rootwindow.geometry("2000x1500")
  rootwindow['bg'] = 'blue'
  
  
  frame_2 = tk.Frame(height=50,width=50)
  frame_1 = tk.Frame()
  frame_3 = tk.Frame()
  frame_4 = tk.Frame()
  frame_5 = tk.Frame()
  
  DamoclesInput = DamoclesInput()
  DamoclesInput.initialise_damocfile_input()
  DamoclesInput.make_clump_button(frame_1)
  
  #create sliders for parameters that are changed by user
  v_slider = Slider(frame_1,v_max_init,(1000, 15000),"Vmax (km/s)",1)
  r_slider = Slider(frame_1,Rrat_init,(0.01, 1),"Rin/Rout",0.0005)
  rho_slider = Slider(frame_1,rho_index_init,(-6, 6),"Density index (\u03B2)",0.01)
  md_slider = Slider(frame_1,mdust_init,(-9, 0.2),"log(Dust mass (M\u2609))",0.001)
  gs_slider = Slider(frame_1,grain_size_init,(-2.3, 0.5),'log(Grain radius (\u03BCm))',0.001)
  amc_frac_slider = Slider(frame_1,0.0,(0.0,1.0),'AmC Fraction',0.01)
  
  
  GasGrid = GasGrid(frame_2)
  GasGrid.initialise_grid_axis(frame_2)
  
  Plotting_window = Plotting_window(frame_3,frame_4,frame_5)
  Plotting_window.initialise_plotwindow(frame_3)
  Plotting_window.make_reset_button(frame_3)
  Plotting_window.make_model_scalebox(frame_5)

  frame_1.place(x=950,y=515)
  frame_2.place(x=980,y=0) 
  frame_3.place(x=100,y=10)
  frame_4.place(x=105,y=865)
  frame_5.place(x=235,y=920)
  

  rootwindow.mainloop()
'''