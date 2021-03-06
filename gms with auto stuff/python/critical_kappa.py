from matplotlib import rc, rcParams
from pylab import arange, tanh, sqrt, real, imag, plot, show, text, xlabel, ylabel, title, axis, savefig, subplots_adjust, figure, semilogy, loadtxt, grid
import os
#from scipy.optimize import broyden1, anderson2

#understanding the roots of equation 18a from the hopf biffurcation in 1d patterns (GMS)
#grabbed stuff from hopf2.m


#--------------------drawing stuff-------------------

fig_width_pt = 546.0  # Get this from LaTeX using 

inches_per_pt = 1.0/72.0                # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches 
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
		  'axes.labelsize': 18,
		  'title.fontsize': 20,
		  'text.fontsize': 16,
		  'xtick.labelsize': 12,
		  'ytick.labelsize': 12,
		  'text.usetex': True,
		  'figure.figsize': fig_size}

rcParams.update(params)
rc('lines', linewidth=1)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))


#--------------------drawing stuff-------------------

Titre = 'hopf_critical_'
kap_p = loadtxt('../data/'+Titre+'_Kappa_full_p.txt')
kap_m = loadtxt('../data/'+Titre+'_Kappa_full_m.txt')
l_p = loadtxt('../data/'+Titre+'_L_critical_plus.txt')
l_m = loadtxt('../data/'+Titre+'_L_critical_minus.txt')
  
figure(1)
plot(kap_p[2:],l_p[2:],'r'), plot(kap_m[2:],l_m[2:],'b')
axis([0.2, 1.5, 4, 14])
grid('True')
show()
text(0.7,10.1,r'\lambda_+',fontsize=20);
text(0.6,9,r'\lambda_-',fontsize=20);
xlabel(r'\kappa_{critical}')
ylabel(r'L_{critical}')

##--------------------------------------------
####this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
Title2 = 'kappa_'+Titre
savefig('../images/'+Title2+'.pdf')
os.system('pdf2ps ../images/'+Title2+'.pdf ../images/'+Title2+'.ps')
os.system('ps2eps -f ../images/'+Title2+'.ps')
os.system('rm ../images/'+Title2+'.ps')