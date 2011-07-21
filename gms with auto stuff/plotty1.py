from pylab import*

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
          'text.fontsize': 12,
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

clf()

A1 = loadtxt('matlab/GMS_rho_0.01_N_750.txt')
A2 = loadtxt('matlab/GMS_rho_0.001_N_750.txt')
A3 = loadtxt('matlab/GMS_rho_0.0001_N_750.txt')
A4 = loadtxt('matlab/GMS_rho_1e-05_N_750.txt')

B1 = loadtxt('branch1')
B2 = loadtxt('branch2')
B3 = loadtxt('branch3')
B4 = loadtxt('branch4')

x1 = B1[:,0] 
x2 = B2[:,0]
x3 = B3[:,0]
x4 = B4[:,0]

y1 = B1[:,3]
y2 = B2[:,3]
y3 = B3[:,3]
y4 = B4[:,3]

x = arange(0,max(x4),0.01)
y = ones(len(x))*0.5603**2

figure(1)
semilogx(A4[:,0],A4[:,4],'0.01',linewidth=3)
semilogx(A3[:,0],A3[:,4],'0.30',linewidth=3)
semilogx(A2[:,0],A2[:,4],'0.60',linewidth=3)
semilogx(A1[:,0],A1[:,4],'0.80',linewidth=3)

leg = legend((r'$\rho=0.00001$',r'$\rho=0.0001$',r'$\rho=0.001$',r'$\rho=0.01$'),loc=6,fancybox=True)
leg.get_frame().set_alpha(0.5)

int1 = find(y1==min(y1))
int2 = find(x1==max(x1))
stable = '-b'
unstable = '--b'

semilogx(B1[0:int1,0],B1[0:int1,3],unstable)
semilogx(B1[int1:int2,0],B1[int1:int2,3],stable)
semilogx(B1[int2:-1,0],B1[int2:-1,3],unstable)
#semilogx(x,y,unstable)

int3 = find(y2==min(y2))
int4 = find(x2==max(x2))

semilogx(B2[0:int3,0],B2[0:int3,3],unstable)
semilogx(B2[int3:int4,0],B2[int3:int4,3],stable)
semilogx(B2[int4:-1,0],B2[int4:-1,3],unstable)

int5 = find(y3==min(y3))
int6 = find(x3==max(x3))

semilogx(B3[int5:-1,0],B3[int5:-1,3],unstable)
semilogx(B3[int6:int5,0],B3[int6:int5,3],stable)
semilogx(B3[0:int6,0],B3[0:int6,3],unstable)

int7 = find(y4==min(y4))
int8 = find(x4==max(x4))

semilogx(B4[int7:-1,0],B4[int7:-1,3],unstable)
semilogx(B4[int8:int7,0],B4[int8:int7,3],stable)
semilogx(B4[0:int8,0],B4[0:int8,3],unstable)

#xlabel(r"$\textrm{L}$")
#ylabel(r"$\textrm{max(V)}$")
xlabel('L')
ylabel('max(V)')
axis([0.25, 20, 0.29, 0.32])

figure(2)
semilogx(B1[0:int1,0],B1[0:int1,3],unstable)
semilogx(B1[int1:int2,0],B1[int1:int2,3],stable)
semilogx(B1[int2:-1,0],B1[int2:-1,3],unstable)
semilogx(x,y,unstable)

int3 = find(y2==min(y2))
int4 = find(x2==max(x2))

semilogx(B2[0:int3,0],B2[0:int3,3],unstable)
semilogx(B2[int3:int4,0],B2[int3:int4,3],stable)
semilogx(B2[int4:-1,0],B2[int4:-1,3],unstable)

int5 = find(y3==min(y3))
int6 = find(x3==max(x3))

semilogx(B3[int5:-1,0],B3[int5:-1,3],unstable)
semilogx(B3[int6:int5,0],B3[int6:int5,3],stable)
semilogx(B3[0:int6,0],B3[0:int6,3],unstable)

int7 = find(y4==min(y4))
int8 = find(x4==max(x4))

semilogx(B4[int7:-1,0],B4[int7:-1,3],unstable)
semilogx(B4[int8:int7,0],B4[int8:int7,3],stable)
semilogx(B4[0:int8,0],B4[0:int8,3],unstable)
semilogx(A4[:,0],A4[:,4],'0.01',linewidth=3)

semilogx(0.3729, y[0],'r.',markersize=12)
semilogx(0.7458, y[0],'r.',markersize=12)
semilogx(1.4917, y[0],'r.',markersize=12)
semilogx(2.9833, y[0],'r.',markersize=12)

xlabel('L')
ylabel('max(V)')
axis([0.25, 20, 0.29, 0.32])

#subplot(2,2,1)
##plot(b1,a1,'rx')
#hist(a[0], bins = n_large/6, normed = 0)
##xlabel(r"$\textrm{eigenvalue}$")
#ylabel(r"$\textrm{frequency}$")
#title(r"$\textrm{B-A graph}$")
#axis([0,25,0,max(a1)])


show()
