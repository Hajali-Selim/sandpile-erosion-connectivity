from scipy.stats import linregress
from inits import *
from itertools import chain
from collections import Counter
import matplotlib.patches as patches
from math import log10, sqrt
from scipy.optimize import curve_fit
from scipy.stats import skew, kurtosis



x1, x2 = 20, 30
G, not_sinks = tgrid(x1,x2)

vP, eP, VP, EP = np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=dict), np.zeros(4, dtype=dict)
for k in range(4):
    vP[k], eP[k], VP[k], EP[k] = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'vegcover.asc', skiprows=6)[1+int(k==3):-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'dem.asc', skiprows=6)[1:-1,1:-1], {}, {}
    for i in range(60):
        for j in range(20):
            VP[k][(j,59-i)], EP[k][(j,59-i)] = vP[k][i,j]/100, eP[k][i,j]

#### 3D MAPS

import matplotlib.animation as animation
from celluloid import Camera
from matplotlib.colors import LightSource, ListedColormap
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)

k = 3
V, E, SC = VP[k], EP[k], pickle.load(open('SCw_p3_40_steps_paths.txt','rb'))

Vmat, Emat, MTmat, SC1mat, SC2mat, FC1mat, FC2mat, FChmat, FCsmat = vP, eP, np.array([eP[k]-np.array([[np.mean(eP[k][i]) for l in range(20)] for i in range(60)]) for k in range(4)]), np.array([nx_to_mat(SC_st[k]) for k in range(4)]), np.array([nx_to_mat(SC_eff[k]) for k in range(4)]), np.array([nx_to_mat(FC_st[k]) for k in range(4)]), np.array([nx_to_mat(FC_eff[k]) for k in range(4)]), DSCHG, SEDTR



Z = Emat[k]
ls = LightSource(270, 45)
fc = ls.shade(FC2mat[k], cmap=plt.cm.terrain_r, vert_exag=0.1, blend_mode='soft')
sc = ls.shade(SC2mat[k], cmap=plt.cm.terrain_r, vert_exag=0.1, blend_mode='soft')
mt = ls.shade(MTmat[k], cmap=plt.cm.rainbow, vert_exag=0.1, blend_mode='soft')
v = ls.shade(Vmat[k], cmap=plt.cm.Greens, vert_exag=0.1, blend_mode='soft')

Z = Emat[k]
ls = LightSource(270, 45)
fc = ls.shade(np.log2(1+FC1mat[k]), cmap=plt.cm.terrain_r, vert_exag=0.1, blend_mode='soft')
sc = ls.shade(np.log2(1+SC1mat[k]), cmap=plt.cm.terrain_r, vert_exag=0.1, blend_mode='soft')
mt = ls.shade(MTmat[k], cmap=plt.cm.rainbow, vert_exag=0.1, blend_mode='soft')
v = ls.shade(Vmat[k], cmap=plt.cm.Greens, vert_exag=0.1, blend_mode='soft')

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect(aspect=(2.5,1,1))
X, Y = np.meshgrid(np.arange(0,20,1), np.arange(0,60,1))

ax.plot_surface(X-60, Y, Z, rstride=1, cstride=1, facecolors=fc, linewidth=0, antialiased=True, shade=False)
ax.plot_surface(X-40, Y, Z, rstride=1, cstride=1, facecolors=sc, linewidth=0, antialiased=True, shade=False)
ax.plot_surface(X-20, Y, Z, rstride=1, cstride=1, facecolors=mt, linewidth=0, antialiased=True, shade=False)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=v, linewidth=0, antialiased=True, shade=False)
ax.xaxis._axinfo["grid"]['color'] = (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] = (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] = (1,1,1,0)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.text(x=18, y=0, z=1.08, s='Vegetation', fontsize=13)
ax.text(x=-2, y=0, z=1.08, s='Microtopography', fontsize=13)
ax.text(x=-21.5, y=0, z=1.08, s='Structural Connectivity', fontsize=13)
ax.text(x=-42, y=0, z=1.08, s='Functional Connectivity', fontsize=13)
ax.plot([19,19], [30, 60], [1, .5], c='tab:red', linewidth=2)
arrow = Arrow3D([19, 19], [59.7, 60.7], [.51, .497], mutation_scale=20, lw=3, arrowstyle="-|>", color="r", shrinkA=0, shrinkB=0)
ax.add_artist(arrow)
plt.show()

import matplotlib as mpl
fig, ax = plt.subplots(figsize=(1,7), tight_layout=True)
cmap = mpl.cm.rainbow
#cmap = mpl.cm.terrain_r.with_extremes(over='blue')
#norm = mpl.colors.BoundaryNorm(list(range(6)), cmap.N, extend='max')
norm = mpl.colors.Normalize(vmin=0, vmax=1)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='vertical', ticks=list(range(6)), extendfrac='auto')
ax.set_yticks([])
plt.show()

0,7,8

V, E = Vs[p][t], Es[p][t]
E, MT = generate_elevation(G, V, Ep, alpha, alpha_std, beta)

Vmat, Emat, MTmat = np.zeros((60,20), int), np.zeros((60,20), float), np.zeros((60,20), float)
for i in range(60):
    for j in range(20):
        Vmat[i,j], Emat[i,j], MTmat[i,j] = V[(j,59-i)], E[(j,59-i)], MT[(j,59-i)]
    for j in range(20):
        MTmat[i,j] = Emat[i,j] - np.mean(Emat[i])

plt.imshow(Vmat, cmap='Greens')
plt.figure()
plt.imshow(MTmat, cmap='rainbow')
plt.show()


### CLUSTERING FIGURES
Vs, Es = np.zeros((3,100), dict), np.zeros((3,100), dict)
for p in range(3):
    pc = .2+p*.3
    for t in range(100):
        Vs[p][t] = tclustering_square(G,.3,pc)
        Es[p][t], _ = generate_elevation(G, Vs[p][t], EP[3], alpha_3, alpha_std_3, beta_3)

Vs, Es = [pickle.load(open('grids/grassland/Vs.txt','rb')), pickle.load(open('grids/shrubland/Vs.txt','rb'))], [pickle.load(open('grids/grassland/Es.txt','rb')), pickle.load(open('grids/shrubland/Es.txt','rb'))]

scp_p4 = pickle.load(open('SCw_p4_40_steps_paths.txt','rb'))
scp_p1 = pickle.load(open('SCw_p1_30_steps_paths.txt','rb'))

plt.imshow(nx_to_mat(scp_p1),cmap='terrain_r')
plt.figure()
plt.imshow(nx_to_mat(scp_p4),cmap='terrain_r')
plt.show()

SCPs = pickle.load(open('time_convergence_data_previous.txt','rb'))[3]


##### PROCESSING GRIDS DATA #####

### SET-UP V_SC_FC (TRY WITH SCE_FCV, LATER, SCP_FCL) + SIZE, and SIZE




#constrained_layout=True
fig, ax = plt.subplots(2,4,figsize=(12,6),tight_layout=True)

fig = plt.figure(figsize=(12,8),constrained_layout=True)
gs = gridspec.GridSpec(2,4, width_ratios=[1,1,1,1.5], figure=fig)
ax = fig.add_subplot(gs[:,0]), fig.add_subplot(gs[:,1]), fig.add_subplot(gs[:,2]), fig.add_subplot(gs[0,3]), fig.add_subplot(gs[1,3])
ax[0].imshow(nx_to_mat(Vs[p][t]), cmap='Greens')
ax[1].imshow(SCmat[p], cmap='terrain_r')
ax[2].imshow(FCmat[p], cmap='terrain_r')
_ = [ax[i].set_xticks([]) for i in range(3)]
_ = [ax[i].set_yticks([]) for i in range(3)]
u1, u2 = np.unique(Siz[p], return_counts=True)
ax[3].scatter(u1,u2,marker='x',s=20,c='k')
ax[3].plot(u1,sum(u2)*plot_weibull3(u1[1:],.999,wb0.alpha,wb0.beta),c='r',lw=.7)
ax[3].loglog([])
ax[3].set_xlabel('avalanche size', size=12)
ax[3].set_ylabel('count', size=12)
ax[4].scatter([i+1 for i in np.unique(Dur[p])], avgSizes[p],marker='x',s=20,c='k')
ax[4].plot([0,max(avgSizes[p])],[0,max(avgSizes[p])], linestyle='--',c='grey',alpha=.7)
ax[4].set_xlabel('duration', size=12)
ax[4].set_ylabel('average size per duration', size=12)
ax[0].set_title(' (a)', loc='left')
ax[1].set_title(' (b)', loc='left')
ax[2].set_title(' (c)', loc='left')
ax[3].set_title(' (d)', loc='left')
ax[4].set_title(' (e)', loc='left')
plt.show()

# LOG BINNING

def drop_zeros(a_list):
    return [i for i in a_list if i>0]

def log_binning(counter_dict,bin_count=35):
    max_x = log10(max(counter_dict.keys()))
    max_y = log10(max(counter_dict.values()))
    max_base = max([max_x,max_y])
    min_x = log10(min(drop_zeros(counter_dict.keys())))
    bins = np.logspace(min_x,max_base,num=bin_count)
    # Based off of: http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
    bin_means_y = (np.histogram(list(counter_dict.keys()),bins,weights=list(counter_dict.values()))[0] / np.histogram(list(counter_dict.keys()),bins)[0])
    bin_means_x = (np.histogram(list(counter_dict.keys()),bins,weights=list(counter_dict.keys()))[0] / np.histogram(list(counter_dict.keys()),bins)[0])
    return bin_means_x,bin_means_y

SS0, SS1 = pickle.load(open('grassland_sizes.txt','rb')), pickle.load(open('shrubland_sizes.txt','rb'))
SSp = np.loadtxt('sizes.txt')


ba_c = dict(Counter(SS0))
ba_x,ba_y = log_binning(ba_c,100)
X, Y = ba_x[~np.isnan(ba_x)], ba_y[~np.isnan(ba_y)]

post_bin0 = np.array([X[i] for i in range(len(X)) for k in range(int(Y[i]))])
post_bin1 = np.array([X[i] for i in range(len(X)) for k in range(int(Y[i]))])

_0b = exponweib.fit(post_bin0)
_1b = exponweib.fit(post_bin1)

_0, _0b
_1, _1b


plt.plot(X, exponweib.pdf(X,*_0b))
plt.loglog([])
plt.show()



skew(post_bin)
kurtosis(post_bin)

skew(SS1)
kurtosis(SS1)


#X, Y = X[1:], Y[1:]
popt, pcov = curve_fit(func3,X,Y)
plt.plot(X,func3(X,*popt),label=str(round(popt[1],2)))
plt.scatter(X,Y)
plt.loglog([])
plt.legend(loc='lower left')
plt.show()

X, Y = ba_x[~np.isnan(ba_x)], ba_y[~np.isnan(ba_y)]
X -= 1.9999

plt.scatter(X,Y, c='k')
popt, pcov = curve_fit(func,X,Y)
mse = round(sqrt(mean_squared_error(Y, func(X,*popt))),2)
plt.plot(X,func(X,*popt),label=str(round(popt[1],2)))
for k in range(10):
    X, Y = X[1:], Y[1:]
    popt, pcov = curve_fit(func,X,Y)
    mse = round(sqrt(mean_squared_error(Y, func(X,*popt))),2)
    plt.plot(X,func(X,*popt),label=str(round(popt[1],2)))

plt.loglog([])
plt.legend(loc='lower left')
plt.show()

post_bin = np.array([X[i] for i in range(len(X)-10) for k in range(int(Y[i]))])
#post_bin = post_bin[np.where(post_bin>2)[0]]
wb2 = Fit_Weibull_3P(post_bin)
plt.figure()
plt.scatter(X,Y,c='k',marker='x',s=30)
Y2 = plot_weibull3(X[:-10], wb2.gamma, wb2.alpha, wb2.beta)
plt.plot(X[:-10],Y2*sum(Y[:-10]),c='r',lw=.7)
plt.loglog([])
plt.show()


def func(x,a,b):
    return a*x**b

def func2(x,a,b):
    return a*np.exp(-x**b)

def func3(x,a,b,loc=1.99999):# a:shape, b:scale, loc:shift=mu
    return (b/a)*(((x-loc)/a)**(b-1))*np.exp(-((x-loc)/a)**b)

def func4(x,c):
    return c*x**(c-1) *np.exp(-x**c)

popt, pcov = curve_fit(func4,X,Y)#,p0=[1/Y[0],1])
plt.plot(X,func4(X,*popt)*sum(X))#,label=str(round(popt[1],2)))
plt.loglog([])
plt.show()

kurtosis(SS0, bias=True)
kurtosis(SS0, bias=False)




plt.scatter(X,Y)
plt.loglog([])
plt.legend(loc='lower left')
plt.show()

plt.plot(SSp, exponweib.pdf(SSp, *exponweib.fit(SSp, 1, 1, scale=2, loc=2)))
plt.loglog([])
plt.show()


plt.plot(post_bin, exponweib.pdf(post_bin, *exponweib.fit(post_bin, 1, 1, scale=2, loc=2)))
plt.loglog([])
plt.show()


plt.plot(SSp, stats.exponweib.pdf(data, *stats.exponweib.fit(data, 1, 1, scale=02, loc=0)))


plt.plot(Y, exponweib.pdf(Y, *exponweib.fit(Y, 2, 1, scale=1, loc=2)))
plt.loglog([])
plt.show()

def wb2LL(p, x): #log-likelihood
    return sum(log(stats.weibull_min.pdf(x, p[1], 0., p[0])))

params1 = weibull_min.fit(post_bin)
params2 = weibull_min.fit(post_bin, floc=1.9999)






X, Y = ba_x[~np.isnan(ba_x)], ba_y[~np.isnan(ba_y)]
alpha, beta = np.polyfit(np.log10(X),np.log10(Y),1)
plt.plot(np.log10(X), [i*alpha+beta for i in np.log10(X)], c='r', lw=2)
plt.scatter(np.log10(X),np.log10(Y),c='k',marker='s',s=50)
plt.loglog([])
plt.show()

slopes, BB = [], dict(Counter(SS1))
for i in range(1,51):
    _ = BB.pop(max(BB.keys()))
    ba_x,ba_y = log_binning(BB,60)
    X, Y = np.log10(ba_x[~np.isnan(ba_x)]), np.log10(ba_y[~np.isnan(ba_y)])
    r = linregress(X,Y)
    slopes.append(r.slope)
    #slopes.append(np.polyfit(X,Y,1)[0])

plt.plot(range(1,51), slopes)
plt.xlabel('removed later sizes')
plt.ylabel('power-law slope')
plt.title('Shrubland, maximal size = 129, number of different sizes = 109')
plt.show()

# for shrubland plot at 0, 25 and 50→keep the binning after deleting last points
BB = dict(Counter(SS1))
ba_x,ba_y = log_binning(BB,50)
plt.scatter(ba_x,ba_y,marker='x',s=40)
X, Y = np.log10(ba_x[~np.isnan(ba_x)]), np.log10(ba_y[~np.isnan(ba_y)])
r = linregress(X,Y)
plt.plot(ba_x, [r.intercept*i**r.slope for i in ba_y], label='0')
for i in range(1,3):
    for k in range(25):
        _ = BB.pop(max(BB.keys()))
    ba_x,ba_y = log_binning(BB,50)
    X, Y = np.log10(ba_x[~np.isnan(ba_x)]), np.log10(ba_y[~np.isnan(ba_y)])
    r = linregress(X,Y)
    plt.plot(ba_x, [r.intercept*i**r.slope for i in Y], label=str(i*25))

plt.legend()
plt.xlabel('sizes')
plt.ylabel('count')
plt.loglog([])
plt.show()



popt, pcov = curve_fit(func2,ba_x[~np.isnan(ba_x)], ba_y[~np.isnan(ba_y)])
plt.scatter(ba_x[~np.isnan(ba_x)],ba_y[~np.isnan(ba_y)],marker='x',s=30)
#plt.plot(ba_x[~np.isnan(ba_x)],func2(ba_x[~np.isnan(ba_x)],*popt))
plt.loglog([])
plt.show()

## CUTTED-OFF POWER-LAW DETECTION
Siz = [sum(Sizes_eff[k],[]) for k in range(2)]
u1, u2 = np.zeros(2, np.ndarray), np.zeros(2, np.ndarray)
u1[0], u2[0] = np.unique(Siz[0], return_counts=True)
u1[1], u2[1] = np.unique(Siz[1], return_counts=True)

# slope stabilization

from scipy.optimize import curve_fit
def func(x,a,b):
    return a*x+b

def func2(x,a,b):
    return a*x**b

popt, pcov = curve_fit(func2,u1[k],u2[k])
plt.scatter(u1[k],u2[k],marker='x',s=30)
plt.plot(u1[k],func2(u1[k],*popt))
plt.loglog([])
plt.show()



for i in range(1,10):
    popt, pcov = curve_fit(func2,u1[k][:-i*20],u2[k][:-i*20])
    plt.figure()
    plt.scatter(u1[k],u2[k],marker='x',s=30)
    plt.plot(u1[k][:-i*20],func2(u1[k][:-i*20],*popt))
    plt.loglog([])
    plt.title('removing last '+str(i*10)+' points, slope='+str(round(popt[1],2))+'\nnumber of points='+str(len(u1[k][:-i*20])))

plt.show()

for i in range(1,10):
    popt, pcov = curve_fit(func2,u1[k][i:],u2[k][i:])
    plt.figure()
    plt.scatter(u1[k],u2[k],marker='x',s=30)
    plt.plot(u1[k][i:],func2(u1[k][i:],*popt))
    plt.loglog([])
    plt.title('skipping first '+str(i)+' points, slope='+str(round(popt[1],2))+'\nnumber of points='+str(len(u1[k][i:])))

plt.show()


slopes = []
for i in range(1,101):
    slopes.append(curve_fit(func2,u1[k][i:],u2[k][i:])[0][1])

plt.plot(range(len(slopes)), slopes, marker='o')
plt.show()



popt, pcov = curve_fit(func,np.log10(u1[k]),np.log10(u2[k]))

plt.scatter(np.log10(u1[k]),np.log10(u2[k]),marker='x',s=30,c='k')
plt.plot(np.log10(u1[k]),func(np.log10(u1[k]),*popt),c='r',lw=.7)
plt.show()





height, width_all = plt.hist(Siz[k], bins=np.logspace(np.log10(1),np.log10(max(Siz[k])), 50))[:-1]
width = width_all[1:]+(width_all[1:]-width_all[:-1])/2




height, width_all = np.histogram(Siz[k])
width = width_all[1:]+(width_all[1:]-width_all[:-1])/2

plt.scatter(width, height)
plt.show()

plt.scatter(u1, u2)
plt.show()


plt.plot(u1[k],func(u1[k],*popt),c='r',lw=.7)
#plt.loglog([])
plt.show()


###


Siz, Dur = [sum(Sizes[p],[]) for p in range(2)], [sum(Durations[p],[]) for p in range(2)]
avgSizes = np.zeros(2, np.ndarray)
for p in range(2):
    u1, u2 = np.unique(Dur[p], return_counts=True)
    avgSizes[p] = np.zeros(len(u1), float)
    for d in range(len(u1)):
        avgSizes[p][d] = np.mean([Siz[p][k] for k in range(len(Siz[p])) if Dur[p][k]==u1[d]])




### Compute Weibull scale and shape


### 3 SCATTERPLOT: SCFC=f(exits), SCFC=f(scale), SCFC=f(shape normalized)





