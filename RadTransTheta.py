#  RadTransTheta.py  - @ by  Iván Martí-Vidal and Miguel Pérez-Torres
#
#    To run the script, type "python RadTransTheta.py" on a Terminal window.
#
#    Requirements: python2.7
#
#    Output: Rff_2D.dat  
#  
#    Code to compute the 2D radiative transfer of synchrotron radiation 
#    observed from a counterjet of 1 unit of intrinsic emission 
#    (emission is in normalized units), taking into account free-free absorption from 
#    a torus and an external medium as described in Mattila & Perez-Torres, et 
#    al. It also takes into account Doppler boosting/de-boosting of the counterjet, 
#    as well as light-travel time effects.
#
#    If you use this script, or a modification of it, we'd
#    appreciate if you gave due credit by citing the following work:
  

# Various imports
import pylab as pl
import numpy as np
import sys
import pickle as pk
import scipy.interpolate as spopt

# Make auxiliary plot to check density distribution?
CHECK_PLOT = False

# Number of pixels of simulation window - the larger, the finer and more
# accurate. 
NPIX = 2048

# Size of the simulation window, in pc
TotSize = 3.0  

# SIZE OF TORUS - Given as twice the torus radius, in pc.
TorusSize = 0.22*2. 


# DENSITY STUFF
# 
########################
# DENSITY IN TORUS, in #/cm^3
TorusDens = 3.7e4 

# DENSITY IN EXTERNAL MEDIUM (AT TorusSize), in #/cm^3:
MedDens = 3.7e2 

# RADIAL DENSITY PROFILE OF EXTERNAL MEDIUM, taken as in Mattila, PT et al.
MedPow = -2.5
########################


# TEMPERATURE STUFF
#
########################
# TEMPERATURE OF INNER TORUS:
T1 = 1.e4 # K
# TEMPERATURE OF EXTERNAL MEDIUM:
T0= 3.e3 # K
# SIZE AT WHICH TEMPERATURE DROPS FROM T1 TO T0:
TSize = 0.195*2. # pc
########################


# RANGE OF INTRINSIC VELOCITIES TO SIMULATE:
Betas = np.linspace(0.01,0.99,100)

# AGE OF THE JET AT SIMULATION:
Age = (56601. - 53300.)*86400.  # seconds

# NORMALIZED INTRINSIC FLUX:
I0 = 1.


# RANGE OF VIEWING ANGLES TO SIMULATE:
Theta = np.linspace(0.1,75.,180) # deg.

#cols = ['r',   'g',  'b',  'k',  'm',  'c']

# Frequency of simulation (Hz):
Nu0 = 8.4e9


##################
# SCRIPT STARTS: #
##################


pc2cgs = 3.09e18

# Number of viewing angles:
NRAY = len(Theta)


# PREPARE GRIDDING:
XX = np.zeros((NPIX,NPIX))

X = np.linspace(-TotSize/2.,TotSize/2.,NPIX)
Y = np.copy(X) # Return an array copy of X, so if X is changed, this doesn't alter Y

# make X a column vector by inserting an axis along 2nd dimension
XX[:] = X[:,np.newaxis]  
YY = np.copy(np.transpose(XX))

RR = np.sqrt(XX**2. + YY**2.)



# PREPARE DENSITY AND TEMPERATURE PROFILES:
Dens = np.zeros((NPIX,NPIX))
Temp = np.copy(Dens)


# Set Density of Torus:
Mask = RR<TorusSize/2.
Dens[Mask] = TorusDens
Out = np.logical_not(Mask)
Dens[Out] = MedDens*np.power(2.*RR[Out]/TorusSize,MedPow)


# Set Temperature of torus:
Mask = RR<TSize/2.
Temp[Mask] = T0
Temp[np.logical_not(Mask)] = T1

# Get trajectory of light rays:
rad2deg = 180./np.pi


# This is the outgoing intensity (for an initial I of 1Jy):
AbsCurve = np.zeros((len(Betas),NRAY))


for bb,Beta in enumerate(Betas):


 Coords = [[] for ti in Theta]

 for t in range(NRAY):

  ViewAng = (Theta[t])/rad2deg
 
# Beta_app of the counterjet (notice + sign in denominator):
  BetaApp = Beta*np.sin(ViewAng)/(1. + Beta*np.cos(ViewAng))

# Initial coordinates (Age is the age of the jet in seconds): 
  X0 = 0.0
  Y0 = -(BetaApp*3.e10*Age)/pc2cgs

# TanT (=Tangent of T) is used to determine DeltaY based on DeltaX
# (i.e., this is to determine the trajetory of the light rays)
  TanT = np.tan((90.-Theta[t])/rad2deg)
  Xini = np.argmin(np.abs(X-X0))


  for i,xi in enumerate(X[Xini:]):

# Distance traveled by light in X direction:
   Dx = xi - X[Xini]

# Distance traveled by light in Y direction:
   Dy = Dx*TanT

# Closest pixel coordinate to Y position of light ray:
   yi = np.argmin(np.abs(Y - Y0 - Dy))

   if i==0:
     Coords[t].append([Xini+i,yi])
   else:
     yi1 = Coords[t][-1][1]

# Fill-in missing pixels in Y direction 
# (i.e., if jump in Dy is of more than 1 pixel size, due
# to the value of Theta[t]):
     if yi>yi1:
      for j in range(1,yi - yi1+1):
       Coords[t].append([Xini+i,yi1+j])
     else:
      for j in range(1,yi1 - yi+1):
       Coords[t].append([Xini+i,yi+j])

# This is just for figure title:
 Dist = np.sqrt(X0**2.+Y0**2.)



# Checking plot (i.e., 2D plot of density + light-rays):
 if CHECK_PLOT:
  PlDens = np.copy(Dens)
  

# Set pixels of light rays to nan (just for plotting!)
# Nans are showed as white pixels)
  for t in range(NRAY):
    for coord in Coords[t]:
      PlDens[coord[0],coord[1]] = np.nan

  pl.imshow(np.transpose(np.log10(PlDens)),origin='lower',extent=(-TotSize/2.,TotSize/2.,-TotSize/2.,TotSize/2.),interpolation='nearest')
  cb = pl.colorbar()
  cb.set_label(r'Log. Density (cm$^{-3}$)')
  pl.xlabel('X (pc)'); pl.ylabel('Y (pc)')
  pl.title('Distance: %.3f pc.'%Dist)
  pl.savefig('RTrans.png')
  pl.show()



  raw_input('CONTINUE? (Ctrl+D to abort; ENTER to accept)')


# Gauntt factor:
 Gauntt = 9.8e-3




# Rad. transfer:

 Idec = [np.zeros(len(Coords[x])) for x in range(NRAY)]

 print '\n\n   GOING TO SIMULATE BETA %.3f\n\n'%Beta

 for t in range(NRAY):
  Ifin = float(I0)

# Pixel size, projected in the direction of propagation of the light ray):
  dpix = TotSize/float(NPIX)*pc2cgs/np.max([np.sin(Theta[t]/rad2deg),np.cos(Theta[t]/rad2deg)])


  for i,coord in enumerate(Coords[t]):

   if not i%100:
     sys.stdout.write('\rRAY %i of %i. Integrating pixel %i of %i   '%(t+1,NRAY,i+1,len(Coords[t])))
     sys.stdout.flush()

   dens = Dens[coord[0],coord[1]]
   temp = Temp[coord[0],coord[1]]

# This is the actual radiative transfer:
   kappa = dens**2.*(temp**(-1.5))*(Nu0**(-2.1))*Gauntt*(17.7 + np.log(temp**1.5/Nu0))
   Ifin -= Ifin*kappa*dpix
   Idec[t][i] = Ifin/I0

# If all I is absorbed, just set it to zero and break the loop:
   if Ifin < 0.0:
     Ifin = 0.0
     Idec[t][i] = 0.0
     break


# Update AbsCurve with the results for the current beta:
 AbsCurve[bb,:] = np.array([Idec[t][-1] for t in range(NRAY)])


# Save all result into an external file:
ofile = open('Rff_2D.dat','w')
pk.dump([Theta,Betas,AbsCurve],ofile)
ofile.close()



