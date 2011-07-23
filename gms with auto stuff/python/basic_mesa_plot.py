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
		  'axes.labelsize': 20,
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

Titre = 'basic_mesa_'
kap1 = loadtxt('../data/'+Titre+'kappa05.txt')
kap2 = loadtxt('../data/'+Titre+'kappa1.txt')
kap3 = loadtxt('../data/'+Titre+'kappa25.txt')
kap4 = loadtxt('../data/'+Titre+'kappa5.txt')
X = loadtxt('../data/'+Titre+'X.txt')
  
figure(1)
plot(X,kap1,'b'), plot(X,kap2,'b')
plot(X,kap3,'b'), plot(X,kap4,'b')
xlabel(r'x')
ylabel(r'u(x)')
text(0.15,1.8,r'\kappa = 0.25',fontsize=18);
text(0.22,1.2,r'\kappa = 1',fontsize=18);
text(0.34,0.75,r'\kappa = 2.5',fontsize=18);
text(0.5,0.35,r'\kappa = 5',fontsize=18);

show()

##--------------------------------------------
###this is in order to generate an eps image, matplotlib seems to be having trouble with that right now
Title2 = 'profiles_'+Titre
savefig('../images/'+Title2+'.pdf')
os.system('pdf2ps ../images/'+Title2+'.pdf ../images/'+Title2+'.ps')
os.system('ps2eps -f ../images/'+Title2+'.ps')
os.system('rm ../images/'+Title2+'.ps')
