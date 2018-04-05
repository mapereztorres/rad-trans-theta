#
# betaThetaPlot-2D.py - @ by  Iván Martí-Vidal and Miguel Pérez-Torres
#
#   Input: Rff_2D.dat  (the output from RadTransTheta.py)
#   Output: BetaTheta_2D.png, .pdf
# 
#   Script to generate the Theta-beta parameter space plot, which constrains the
#   allowed values of beta and theta for an expanding jet. 
#   It uses constraints from radio interferometric observations. Namely the 
#   apparent speed of the jet, as well as the jet-to-counterjet radio, R. 
#   The values in the code apply to the TDE Arp299B-AT1, and follows the
#   observations reported in Mattila, Pérez-Torres et al., but the code is of
#   general use and can thus be used for any other jet.   
# 
#   If you use this script, or a modification of it, we'd
#   appreciate if you gave due credit by citing the following work:
# 
#   Mattila, Perez-Torres et al. (2018)

import pylab as pl
import numpy as np
import pickle as pk
import scipy.interpolate as spopt

### Initializing values

# Gridding in beta direction:
Nval = 1000

# p = k - alpha, with k=2 (continue jet) and k=3 for a discrete jet,
# and alpha the observed spectral index of the jet, defined as (S_nu ~ nu^(-alpha))
# The value below assumes k=3 and alpha = 1
p = 2.0 

# List of values constraining Bapp: 
BetaApps = [0.16,0.34]

# List of values constraining the (observed) brightness ratio (the value below corresponds to 5 sigma!):
Rats = [2.4]

# Plot dashed line at values of beta_intrinsic (modelled) within 1 sigma 
BetaMod = 0.22 
BetaModSigma = 0.02

#### Values for the plot
# Coordinates of labels:
CR1 = [0.49,60.] # X and Y of label for R curve

#CA1 = [0.02,40.] # X and Y of label for lower beta_app
CA1 = [0.30,4.0] # X and Y of label for lower beta_app
CA2 = [0.45,30.] # X and Y of label for higher beta_app

# Text (and axis-labels) fonts:
fs = 20 
fontname = 'Latin Modern Roman'

#################
# SCRIPT STARTS #
#################


# Open aux file (written by RadTransTheta.py) 
# and create interpolation function:

ifile = open('Rff_2D.dat')
Theta,Betas,AbsCurve = pk.load(ifile)
ifile.close()

# Interpolation function. It returns the fraction of remaining intensity
# that passed the free-free screen:
Rextra = spopt.RectBivariateSpline(Betas,Theta,AbsCurve,kx=3,ky=3)

# Array of betas for the plot:
Beta = np.linspace(0.010,0.999,Nval)

rad2deg = 180./np.pi

# This function returns be intrinsic beta, given viewing angle AND 
# brightness ratio. It accounts for the free-free effects as well
# (i.e., it decouples the Rextra from the Doppler boosting/deboosting):
def BetaThetaR(theta,R):

# theta is an array. This function returns an array as well.
  print 'Doing interpolation'
  output = np.zeros(len(theta))
  NITERMAX = 1000
# Loop to get a self-consistent beta:
  for i in range(len(theta)):
    BetaTest = 0.0
    BetaFin = 0.5
    NITER = 0
    while (np.abs(BetaTest-BetaFin)>0.005 and NITER<NITERMAX):

# Get the free-free contribution for BetaTest:
      BetaTest = BetaFin
      R1 = Rextra(BetaTest,theta[i])

# Doppler contribution to the (observed) constraint in brightness ratio:
      Rdeboost = np.max([R*R1,1.0])

# Intrinsic beta given Rdeboost:
      Rpow = np.power(Rdeboost,1./p)
      BetaFin = ((Rpow-1)/(Rpow+1))/np.cos(theta[i]/rad2deg)

      NITER += 1

    if NITER>=NITERMAX:
      print 'WARNING! NO CONVERGENCE FOR %i: %.2f'%(i,BetaTest-BetaFin)

    output[i] = BetaFin

  return output

#def ThetaBetaR(beta,R):
#  Rpow = np.power(R,1./p)
#  return rad2deg*np.arccos(1./beta*((Rpow-1)/(Rpow+1)))

# Given intrinsic beta and apparent beta, it returns the viewing angle.
# It accounts for the maximum allowed beta_app (i.e., no solution of theta 
# out for allowed region):
def ThetaBetaAppP(beta,betaApp):
  result = np.zeros(len(beta))
  mask = betaApp<beta/np.sqrt(1.-beta**2.)
  result[mask] = (betaApp**2. + np.sqrt((beta[mask]**2.-1)*betaApp**2. + beta[mask]**2.))/(beta[mask]*betaApp**2. + beta[mask])
  return rad2deg*np.arccos(result)


# This is the second solution for theta (but not interesting for us,
# since it covers the range theta>90.:

#def ThetaBetaAppM(beta,betaApp):
#  result = np.zeros(len(beta))
#  mask = betaApp<beta/np.sqrt(1.-beta**2.)
#  result[mask] = (betaApp**2. - np.sqrt((beta[mask]**2.-1)*betaApp**2. + beta[mask]**2.))/(beta[mask]*betaApp**2. + beta[mask])
#  return rad2deg*np.arccos(result)

fig = pl.figure(figsize=(8,6))
sub1 = fig.add_subplot(111)
fig.subplots_adjust(left=0.11,right=0.97,top=0.97,bottom=0.11)

# For the brightness curve, we use "theta" as input, and derive beta: 
Thetas = np.linspace(0.,70.,200)
BR1 = BetaThetaR(Thetas,Rats[0])
pl.plot(BR1,Thetas,'--b',linewidth=3.0)

# This is just for the shadowed area. We re-sample the R curve homogeneously in beta,
# using an interpolation function:
ThetaInterp = np.linspace(31.,90.,70)
BRInterp = BetaThetaR(ThetaInterp,Rats[0])
Binterp = spopt.interp1d(BRInterp,ThetaInterp,kind='linear',bounds_error=False)
FillB = Binterp(Beta)

# Curves for beta_app:
BA1 = ThetaBetaAppP(Beta,BetaApps[0]) ; mask1 = BA1>0.0
pl.plot(Beta[mask1],BA1[mask1],'-r',linewidth=3)

BA2 = ThetaBetaAppP(Beta,BetaApps[1]) ; mask2 = BA2>0.0
pl.plot(Beta[mask2],BA2[mask2],'-r',linewidth=3)

pl.plot([BetaMod-3*BetaModSigma,BetaMod-3*BetaModSigma],[0.,80.],'-.k',linewidth=2.5)
pl.plot([BetaMod+3*BetaModSigma,BetaMod+3*BetaModSigma],[0.,80.],'-.k',linewidth=2.5)

pl.xlabel(r'$\beta$',fontsize=fs,family=fontname)
pl.ylabel(r'$\theta$ (deg.)',fontsize=fs,family=fontname)

pl.ylim((0,70.))
#pl.xlim((0.001,1.0))
pl.xlim((0.0,1.0))

pl.setp(sub1.get_xticklabels(),fontsize=fs)
pl.setp(sub1.get_yticklabels(),fontsize=fs)
pl.setp(sub1.get_xticklabels(),family=fontname)
pl.setp(sub1.get_yticklabels(),family=fontname)

# This is for the shadowing:
C1 = np.minimum(FillB,BA2)
C2 = BA1

# Print the values of beta_app and R in the figure:
pl.text(CR1[0],CR1[1],'R = %.2f'%Rats[0],fontsize=fs,family=fontname,color='b')
pl.text(CA1[0],CA1[1],r'$\beta_{app}$ = %.2f'%BetaApps[0],fontsize=fs,family=fontname,color='r')
pl.text(CA2[0],CA2[1],r'$\beta_{app}$ = %.2f'%BetaApps[1],fontsize=fs,family=fontname,color='r')

X0 = np.argmin(np.abs(0.22-Beta))
pl.fill_between(Beta[X0:],C1[X0:],C2[X0:],facecolor='green',alpha=0.25)

# Change linewidth of axes and ticks
for axis in ['top','bottom','left','right']:
  sub1.spines[axis].set_linewidth(2.0)

sub1.tick_params(length=10, width=2)

pl.tight_layout()

pl.savefig('BetaTheta_2D.png')
pl.savefig('BetaTheta_2D.pdf')

pl.show()
