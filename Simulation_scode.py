
import Tkinter as tk    # Importing the Tkinter (tool box) library 
from PIL import ImageTk
import numpy as np
import Image as im
import math as ma
import cmath as cma
import pylab as pl
import Pade1




print "\n \n ###################### Start ################################## \n\n"



###### simulation settings start !!!! 
height_im = 100
# z le_pz=100*le_px --> one pixel in z direction = 1um
width_im=1200

cladding_color=[10,10,10]
n2=1.55
core_color=[255,255,255]
n1=1.60

# 1 pixel = 0.1um
le_px=1*ma.pow(10,-7)


#BPM constants
BPM_freq =250*ma.pow(10,12)
BPM_A = 1
# Eworld - eigenvalue for simulation world
error=1*ma.pow(10,-17);
hzMode =0

# Create black empty images
blank_image = np.zeros((height_im,width_im,3), np.uint8)
#blank_image[:,0:0.5*width] = (255,0,0)      # (R, G, B)
#blank_image[:,0.5*width_im:width_im] = (0,255,0)
blank_image[:,:] = (0,0,0)


#Draw a waveguide
waveguide_image = np.zeros((height_im,width_im,3),np.uint8)
#Cladding
waveguide_image[:,:]=cladding_color #-> n=1.6
#Waveguide
zz1=0.05*width_im
zz2=0.4*width_im
x1=0.12*height_im
x2=0.24*height_im
################ upper waveguide
waveguide_image[(x1):(x2),0:zz1]=core_color #

#Draw a right narrow curve

xx1=x1
xx2=x2

while (xx2<0.49*height_im) and (zz1<(width_im-13)):
   
   waveguide_image[xx1:xx2,zz1]=core_color 
   waveguide_image[xx1:xx2,zz1+1]=core_color 
   waveguide_image[xx1:xx2,zz1+2]=core_color 
   zz1=zz1+3
   xx1=xx1+1
   xx2=xx2+1
   
waveguide_image[(xx1-1):(xx2-1),zz1:width_im]=core_color

#Draw a right width curve

##while (xx2<0.51*height_im) and (zz1<(width_im-13)):
##   
##   waveguide_image[xx1:xx2,zz1]=core_color 
##   waveguide_image[xx1:xx2,zz1+1]=core_color 
##   waveguide_image[xx1:xx2,zz1+2]=core_color 
##   waveguide_image[xx1:xx2,zz1+3]=core_color 
##   waveguide_image[xx1:xx2,zz1+4]=core_color 
##   waveguide_image[xx1:xx2,zz1+5]=core_color
##   zz1=zz1+5
##   xx1=xx1+1
##   xx2=xx2+1

##############################   lower waveguide
zz1=0.05*width_im

x3=0.8*height_im
x4=0.92*height_im
waveguide_image[ (x3):(x4),0:zz1]=core_color # -> n=1.5

#Draw a left narrow curve

xx3=x3
xx4=x4

while (xx3>0.52*height_im) and (zz1<(width_im-13)):
   
   waveguide_image[xx3:xx4,zz1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+2]=core_color # -> n=1.5
   zz1=zz1+3
   xx3=xx3-1
   xx4=xx4-1
zz2=zz1+0.0833*width_im
waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

#Draw a right narrow curve
x3=xx3+1
x4=xx4+1


while (x4<0.92*height_im) and (zz2<(width_im-13)):
   
   waveguide_image[x3:x4,zz2]=core_color 
   waveguide_image[x3:x4,zz2+1]=core_color 
   waveguide_image[x3:x4,zz2+2]=core_color 
   waveguide_image[x3:x4,zz2+3]=core_color 
   zz2=zz2+3
   x3=x3+1
   x4=x4+1

xx3=x3-2
xx4=x4-2
zz1=zz2
zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

zz1=zz2

#draw left curve

while (xx3>0.6*height_im) and (zz1<(width_im-13)):
   
   waveguide_image[xx3:xx4,zz1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+2]=core_color # -> n=1.5
   zz1=zz1+3
   xx3=xx3-1
   xx4=xx4-1

zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

x3=xx3+1
x4=xx4+1

#draw right curve
while (x4<0.92*height_im) and (zz2<(width_im-13)):
   
   waveguide_image[x3:x4,zz2]=core_color 
   waveguide_image[x3:x4,zz2+1]=core_color 
   waveguide_image[x3:x4,zz2+2]=core_color 
   waveguide_image[x3:x4,zz2+3]=core_color 
   zz2=zz2+3
   x3=x3+1
   x4=x4+1

xx3=x3-2
xx4=x4-2
zz1=zz2
zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

zz1=zz2

#draw left curve
while (xx3>0.6*height_im) and (zz1<(width_im-13)):
   
   waveguide_image[xx3:xx4,zz1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+2]=core_color # -> n=1.5
   zz1=zz1+3
   xx3=xx3-1
   xx4=xx4-1

zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

x3=xx3+1
x4=xx4+1

#draw right curve
while (x4<0.92*height_im) and (zz2<(width_im-13)):
   
   waveguide_image[x3:x4,zz2]=core_color 
   waveguide_image[x3:x4,zz2+1]=core_color 
   waveguide_image[x3:x4,zz2+2]=core_color 
   waveguide_image[x3:x4,zz2+3]=core_color 
   zz2=zz2+3
   x3=x3+1
   x4=x4+1

xx3=x3-2
xx4=x4-2
zz1=zz2
zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

zz1=zz2

#draw left curve
while (xx3>0.52*height_im) and (zz1<(width_im-13)):
   
   waveguide_image[xx3:xx4,zz1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+1]=core_color # -> n=1.5
   waveguide_image[xx3:xx4,zz1+2]=core_color # -> n=1.5
   zz1=zz1+3
   xx3=xx3-1
   xx4=xx4-1

zz2=zz1+150

waveguide_image[ (xx3+1):(xx4+1),zz1:zz2]=core_color

x3=xx3+1
x4=xx4+1

#draw right curve
while (x4<0.92*height_im) and (zz2<(width_im-13)):
   
   waveguide_image[x3:x4,zz2]=core_color 
   waveguide_image[x3:x4,zz2+1]=core_color 
   waveguide_image[x3:x4,zz2+2]=core_color 
   waveguide_image[x3:x4,zz2+3]=core_color 
   zz2=zz2+3
   x3=x3+1
   x4=x4+1

xx3=x3-2
xx4=x4-2
zz1=zz2
zz2=zz1+10

waveguide_image[ (xx3+1):(xx4+1),zz1:width_im]=core_color





####
####   


######## simulation settings end

#print 0.48*height_im
#print 0.52*width_im

# show simulation world
image=waveguide_image
#pic = im.fromarray(image, 'RGB')
#pic.save('image.ppm')
#pic.save('image.jpg')
#master = tk.Tk()
#master.title("My World")
#canvas = tk.Canvas(master, width= width_im, height=height_im)
#canvas.pack()
#img = tk.PhotoImage(file="image.ppm")
#canvas.create_image(0,0, anchor=tk.NW, image=img)
#tk.mainloop()

    # ############################## Direction #############################

left_to_right = 1
dir = 0

#which direction is wave travelling

dir = left_to_right


if (dir):
   dir=dir
   #print "wave travels from left to right -->"
else:
    print "!!! check direction !!!"


#if ((raw_input())=='y'):
#    print 'Finish'


# ############################## Analyse structure #############################

# 

def function_colordecoder(color):
   if np.array_equal(color,[10, 10, 10]):
      return n2
   if np.array_equal(color,[255, 255, 255]):
      return n1



n_mesh_ar=np.zeros((height_im,width_im,1),np.float64)
image=waveguide_image
colorchange_list=[]
old_value=waveguide_image[0,0]
i=0
#create array with mesh of n for structure
while i<height_im:
   ii=0
   while ii<width_im:
       new_value= waveguide_image[i,ii]
       n_mesh_ar[i,ii]=function_colordecoder(new_value)
       ii=ii+1
   i=i+1


#find start and end of cladding + core
i=0
while i<height_im:
    
   new_value= waveguide_image[i,0]
   n_mesh_ar[i,0]=function_colordecoder(new_value)
   
   if not(np.array_equal(new_value,old_value)):
      colorchange_list.append(i)
   old_value=new_value
   i=i+1
   

#print "1 px=",le_px,"m"
print "1 px \n             ->",le_px*ma.pow(10,6),"um \n\n"

#higher cladding
claddingH=[colorchange_list[1]+1, height_im]

#core
core=[colorchange_list[0], colorchange_list[1]]

#lower cladding
claddingL=[0, colorchange_list[0]-1]
   
# calculate width of core
d_core=(core[1]-core[0])*le_px/2
print "d core of waveguide \n             ->",d_core

#center of core (px)
center_for_BPM =  core[0]+ (core[1]-core[0])/2


# ########################## Mode solving for TE-Mode #############################

# calculate Mode TE 1



#Some constants
c_vacuum =2.998
c_vacuum =c_vacuum*ma.pow(10,8)
#print "c=",c_vacuum

print "\n \n++++++ Start Mode analysis ++++++++++++++++++++++++++++++++++++++++\n \n"

# #######################################Frequency for analysis

#given Frequency
f_given = 10*ma.pow(10,12)
# print "Frequency for analysis  \n             ->",f_given*ma.pow(10,-12),"THz"


# ###################################### Variables for analysis of world
n1=n1 # core refraction index
n2=n2 #cladding refraction index

#give wavelenght
L_vacuum = c_vacuum/f_given
#k0- wavevector
k_vacuum = 2*ma.pi/L_vacuum
# k -vector core
k1= n1*k_vacuum
# k -vector cladding
k2= n2*k_vacuum
# Dn^2 constant for simulation world
#Dnworld= ma.pow(n1,2)-ma.pow(n2,2)

# Vworld -  constant for simulation world
#Vworld= (2*ma.pi*d_core)*(1/L_vacuum)*ma.sqrt(Dnworld)
#print "V constant for simulation world \n             ->", Vworld



# Functions used by mode solver

def function_find_neff (E_error, E_p,  E_n1, E_n2, E_k0,E_d):
   "Calculates the eigenvalue for the waveguide by using bisection"

   #initialise variables
   E_neff_max = E_n1
   E_neff_min =  E_n2
   E_neff_up = E_neff_max
   E_neff_down = E_neff_min
   E_error_apriori = (E_neff_max-E_neff_min)/2
   i=0

   while (E_error<=E_error_apriori):
         
          E_neff_center = E_neff_down+ (E_neff_up - E_neff_down)/2
          
          if (function_neff(E_neff_up,E_p,  E_n1, E_n2, E_k0,E_d)*function_neff(E_neff_center,E_p,  E_n1, E_n2, E_k0,E_d)<=0):
              E_neff_down=E_neff_center     
          elif (function_neff(E_neff_down,E_p,  E_n1, E_n2, E_k0,E_d)*function_neff(E_neff_center,E_p,  E_n1, E_n2, E_k0,E_d)<=0):
              E_neff_up=E_neff_center
              
          E_error_apriori= E_error_apriori/2
          if i>ma.pow(10,4):
                print"break loop"
                print "apriori error ",E_error_apriori
                print "eigenvalue" ,E_neff_center
                break
               
          #control if E_neff_up and E_neff_down inside boundaries
          if E_neff_up>E_neff_max:
             print "++++++++++++++ E_neff_up out of boundaries !   find_neff"
             E_neff_up=E_n1
          if E_neff_down<E_neff_min:
             print "++++++++++++++ E_neff_down out of boundaries !    find_neff"
             E_neff_down=E_n2
             
   return E_neff_center

def function_neff(f_neff ,f_p,  f_n1, f_n2, f_k0,f_d):
   "calculates one value of the eigenvalue function with effective index"
   # atan(inf)=pi/2
   if f_n1==f_neff:
      return ma.pi/2
   
   t = ma.pow(f_neff,2)-ma.pow(f_n2,2)
   t= t/ (ma.pow(f_n1,2)-ma.pow(f_neff,2))
   t = ma.atan(ma.sqrt(t))
   t = t + f_p*ma.pi/2 - f_k0*f_d*ma.sqrt(ma.pow(f_n1,2)-ma.pow(f_neff,2))
   return t

#init

neff_list=[]
freq_list=[]

f_given=BPM_freq
L_vacuum = c_vacuum/f_given
k_vacuum = 2*ma.pi/L_vacuum
k1= n1*k_vacuum
k2= n2*k_vacuum
error=1*ma.pow(10,-15);
mode =hzMode
neffective = function_find_neff(error,mode,n1,n2,k_vacuum,d_core)
BPM_neff=neffective

f_given = 20*ma.pow(10,12)
#Loop to calculate neff vs freq.
while f_given<=1000*ma.pow(10,12):
   L_vacuum = c_vacuum/f_given
   k_vacuum = 2*ma.pi/L_vacuum
   k1= n1*k_vacuum
   k2= n2*k_vacuum
   error=1*ma.pow(10,-15);
   mode =hzMode
   neffective = function_find_neff(error,mode,n1,n2,k_vacuum,d_core)
   neff_list.append(neffective)
   freq_list.append(f_given*1*ma.pow(10,-12))
   f_given=f_given + 20*ma.pow(10,12)
#plot
#pl.figure(211)
#pl.plot(freq_list,neff_list)
#pl.xlabel("f/Thz")
#pl.ylabel("neff")
#plot_title = "mode "+ str(mode)
#pl.title(plot_title)
#pl.draw()   

########################################## BPM Simulation ###################################################
print "\n \n+++++++++ Start BPM sim ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n \n"

print "Frequency for BPM \n                     --> ",BPM_freq*ma.pow(10,-12),"THZ"
print "effective index for BPM \n                     --> ",BPM_neff
print "A (amplitude) for BPM \n                     --> ",BPM_A

print "\n ++++++++ Calc initial mode \n"

L_vacuum = c_vacuum/BPM_freq
k_vacuum = 2*ma.pi/L_vacuum
k1= n1*k_vacuum
k2= n2*k_vacuum
BPM_bheta= k_vacuum*BPM_neff

#Hz(d) =  boundary condition, cos 
BPM_B= BPM_A*ma.cos(d_core*ma.sqrt(ma.pow(k1,2)-ma.pow(k_vacuum*BPM_neff,2)))



######## Functions to calculate envelope at z=0
def function_Hz_core(x,k1,bheta,A,d):
   "calculate the value of Hz in the core abs(x)<=d for z=0 depending on x"
   if abs(x)>d:
      print "++++++++++++++ x out of boundaries ! function_Hz_core" 
      return 0
   Hz_core=A*ma.cos(x*ma.sqrt(ma.pow(k1,2)-ma.pow(bheta,2)))
   return Hz_core

def function_Hz_cladding(x,k2,bheta,B,d):
   "calculate the value of Hz in the cladding abs(x)>d for z=0 depending on x"
   if abs(x)<=d:
      print "++++++++++++++ x out of boundaries ! function_Hz_cladding" 
      return 0
   if x>d:
      Hz_cladding=B*ma.exp(-(x-d)*ma.sqrt(ma.pow(bheta,2)-ma.pow(k2,2)))
      return Hz_cladding
   if x<d:
      Hz_cladding=B*ma.exp((x+d)*ma.sqrt(ma.pow(bheta,2)-ma.pow(k2,2)))
      return Hz_cladding

####### calculate envelope and plot
#zeropoint of calculation --> center of waveguide

#start at the bottom and go up
x_axis=[]
y_axis=np.zeros(height_im)
yAxList=[]
#x_normalized = center of core represents x=0 an one px = a*le
#x_normalized= (x_image - center_for_BPM)*le_px

x_counter=0

while x_counter<height_im:
   x_axis.append(x_counter)

   #normalize x
   x_normalized=(x_counter - center_for_BPM)*le_px

   if abs(x_normalized)>d_core:
      y_axis[x_counter]=function_Hz_cladding(x_normalized,k2,BPM_bheta,BPM_B,d_core)
      yAxList.append(function_Hz_cladding(x_normalized,k2,BPM_bheta,BPM_B,d_core))
   if abs(x_normalized)<=d_core:
      y_axis[x_counter]=function_Hz_core(x_normalized,k1,BPM_bheta,BPM_A,d_core)
      yAxList.append(function_Hz_core(x_normalized,k1,BPM_bheta,BPM_A,d_core))
   x_counter=x_counter+1

pl.figure(700)
pl.plot(x_axis,y_axis)

pl.draw()
pl.show()


# def PadeApproximation(startAr ,nref, freq, nMeshAr,xSize,zSize,xStep,zStep):
t=Pade1.PadeApproximation(y_axis, BPM_neff, BPM_freq, n_mesh_ar, height_im, width_im, le_px,le_px*10)    
print t
