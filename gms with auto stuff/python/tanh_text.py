from matplotlib import rc, rcParams
from pylab import arange, tanh, sqrt, real, imag, plot, show, text, xlabel, ylabel, title, axis, savefig, subplots_adjust, figure, semilogy, loadtxt, grid
import os
#from scipy.optimize import broyden1, anderson2

#understanding the roots of equation 18a from the hopf biffurcation in 1d patterns (GMS)
#grabbed stuff from hopf2.m

# Have to fix it a bit to generate the two figure, for kappa=.65 and kappa=1


#--------------------drawing stuff-------------------

fig_width_pt = 546.0  # Get this from LaTeX using 

inches_per_pt = 1.0/72.0                # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches 
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
		  'axes.labelsize': 24,
		  'title.fontsize': 20,
		  'text.fontsize': 16,
		  'xtick.labelsize': 15,
		  'ytick.labelsize': 15,
		  'text.usetex': True,
		  'figure.figsize': fig_size}

rcParams.update(params)
rc('lines', linewidth=1)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))


#--------------------drawing stuff-------------------

Titre = 'kappa_1'
lam_p = loadtxt('../data/'+Titre+'_lambda_plus.txt')
lam_m = loadtxt('../data/'+Titre+'_lambda_minus.txt')
tau_p = loadtxt('../data/'+Titre+'_tau_plus.txt')
tau_m = loadtxt('../data/'+Titre+'_tau_minus.txt')
LL = loadtxt('../data/'+Titre+'_L.txt')
  
#figure(1)
#plot(LL,lam_p,'r'), plot(LL,lam_m,'b')
#axis([0, 14, -1e-6, 1e-6])
#grid('True')
#text(2,-37e-8,r'\lambda_+',fontsize=24);
#text(4,-13e-8,r'\lambda_-',fontsize=24);
#xlabel(r'L')
#ylabel(r'\lambda_{\pm}')

##--------------------------------------------
####this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
#Title1 = 'lambda_L_'+Titre
#savefig('../images/'+Title1+'.pdf')
#os.system('pdf2ps ../images/'+Title1+'.pdf ../images/'+Title1+'.ps')
#os.system('ps2eps -f ../images/'+Title1+'.ps')
#os.system('rm ../images/'+Title1+'.ps')

figure(2)
semilogy(LL,tau_p,'r'), semilogy(LL,tau_m,'b')
semilogy(LL,-tau_p,'r--'), semilogy(LL,-tau_m,'b--')
text(2,4e3,r'\tau_+',fontsize=24);
text(4,6000,r'\tau_-',fontsize=24);
text(11,56,r'-\tau_+',fontsize=24);
text(13,250,r'-\tau_-',fontsize=24);
xlabel(r'L')
ylabel(r'\tau_{\pm}')
axis([0, 14, 10, 1e8])
grid('True')
show()

##--------------------------------------------
###this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
Title2 = 'tau_L_'+Titre
savefig('../images/'+Title2+'.pdf')
os.system('pdf2ps ../images/'+Title2+'.pdf ../images/'+Title2+'.ps')
os.system('ps2eps -f ../images/'+Title2+'.ps')
os.system('rm ../images/'+Title2+'.ps')

###---------------------------------------------------
###---------------------------------------------------
###---------------------------------------------------

Titre = 'kappa_0.65'
lam_p = loadtxt('../data/'+Titre+'_lambda_plus.txt')
lam_m = loadtxt('../data/'+Titre+'_lambda_minus.txt')
tau_p = loadtxt('../data/'+Titre+'_tau_plus.txt')
tau_m = loadtxt('../data/'+Titre+'_tau_minus.txt')
LL = loadtxt('../data/'+Titre+'_L.txt')
  
#figure(3)
#plot(LL,lam_p,'r'), plot(LL,lam_m,'b')
#axis([0, 14, -1e-6, 1e-6])
#grid('True')
#text(2,-36e-8,r'\lambda_+',fontsize=24);
#text(4,-9e-8,r'\lambda_-',fontsize=24);
#xlabel(r'L')
#ylabel(r'\lambda_{\pm}')

##--------------------------------------------
####this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
#Title1 = 'lambda_L_'+'kappa_065'
#savefig('../images/'+Title1+'.pdf')
#os.system('pdf2ps ../images/'+Title1+'.pdf ../images/'+Title1+'.ps')
#os.system('ps2eps -f ../images/'+Title1+'.ps')
#os.system('rm ../images/'+Title1+'.ps')

#figure(4)
#semilogy(LL,tau_p,'r'), semilogy(LL,tau_m,'b')
#semilogy(LL,-tau_p,'r--'), semilogy(LL,-tau_m,'b--')
#text(2,1e4,r'\tau_+',fontsize=24);
#text(3,24e3,r'\tau_-',fontsize=24);
#text(11,22e1,r'-\tau_-',fontsize=24);
#text(12,5e2,r'-\tau_+',fontsize=24);
#xlabel(r'L')
#ylabel(r'\tau_{\pm}')
#axis([0, 14, 10, 1e8])
#grid('True')
#show()

##--------------------------------------------
###this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
#Title2 = 'tau_L_'+'kappa_065'
#savefig('../images/'+Title2+'.pdf')
#os.system('pdf2ps ../images/'+Title2+'.pdf ../images/'+Title2+'.ps')
#os.system('ps2eps -f ../images/'+Title2+'.ps')
#os.system('rm ../images/'+Title2+'.ps')