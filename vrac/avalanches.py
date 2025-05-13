from avalanche_packages import *
#from step import *
from inits import *
#from SC import *
#from FC import *
import sympy as sp
from matplotlib import gridspec
import pandas as pd

def init2(g):
    state, coupling, inter_dist, current_av, branches, new_active, size, dur, ordlist = {}, {}, {}, [], {}, [], [], [], sorted(list(g), key=lambda a:a[1], reverse=True)
    for i in g:
        state[i], coupling[i], inter_dist[i] = 0, {}, {}
    for i,j in itertools.combinations(ordlist, 2):
        if nx.has_path(g,i,j):
            coupling[i][j], inter_dist[i][j] = [], nx.shortest_path_length(g,i,j)
        if i[1] == j[1]:		# for tgrid
            if nx.has_path(g,j,i):
                coupling[j][i], inter_dist[j][i] = [], nx.shortest_path_length(g,j,i)
    return state, coupling, current_av, branches, new_active, size, dur, inter_dist

### PLOTTING ADDITIONAL MAHLERAN FCS (FOR G, G/S AND S/G)
DSCHG, SEDTR = np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=np.ndarray)
for k in range(4):
    DSCHG[k], SEDTR[k] = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'_rainA_highsm_dschg.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'_rainA_highsm_sedtr.asc', skiprows=6)[1:-1,1:-1]

for k in range(4):
    plt.imshow(DSCHG[k], cmap='rainbow')
    plt.title('discharge P'+str(k+1))
    plt.figure()
    plt.imshow(SEDTR[k], cmap='rainbow')
    plt.title('sediment P'+str(k+1))
    plt.figure()

# DSCH to dict
DSCHG2, SEDTR2 = np.zeros(4, dtype=dict), np.zeros(4, dtype=dict)
for k in range(4):
    DSCHG2[k], SEDTR2[k] = {}, {}
    for i,j in G:
        DSCHG2[k][i,59-j], SEDTR2[k][i,59-j] = DSCHG[k][j,i], SEDTR[k][j,i]

VV, SCI, SCO, SC = np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float)
FCI, FCO, FC = np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float)
for i in range(60):
    for j in range(20):
        VV[i,j], SCI[i,j], SCO[i,j] = VP[3][(j,59-i)], SCin[(j,59-i)], SCout[j,59-i]
        SC[i,j] = SCI[i,j]-SCO[i,j]
        FCI[i,j], FCO[i,j] = FCin[(j,59-i)], FCout[j,59-i]
        FC[i,j] = FCI[i,j]-FCO[i,j]

vP, eP, VP, EP = np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=dict), np.zeros(4, dtype=dict)
for k in range(4):
    vP[k], eP[k], VP[k], EP[k] = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'vegcover.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'dem.asc', skiprows=6)[1:-1,1:-1], {}, {}
    for i in range(60):
        for j in range(20):
            VP[k][(j,59-i)], EP[k][(j,59-i)] = vP[k][i,j]/100, eP[k][i,j]

SCFCt, FCFCt = np.zeros((nb_runs, 100), dtype=float), np.zeros((nb_runs, 100), dtype=float)
for run in range(nb_runs):
    print('run n°',run)
    S, Cst, current, branches, active = zeros(Cst)
    for step in range(nb_steps+1):
        S, Cst, current, branches, active, _, _ = sandpilest(Gst, VP[3], S, Cst, current, branches, active, [], [], not_sinks)
        if step % 100 == 0 and step>0:
            #step
            FC = FC_path(Cst)
            SCFCt[run, int(step/100-1)], FCFCt[run, int(step/100-1)] = SCFC(SCst,FC,not_sinks), SCFC(FC,DSCHG2[3],not_sinks)

for run in range(nb_runs):
    print('run n°',run)
    S, Cst, current, branches, active = zeros(Cst)
    for step in range(nb_steps+1):
        S, Cst, current, branches, active, _, _ = sandpilest(Gst, VP[3], S, Cst, current, branches, active, [], [], not_sinks)
    FC = FC_path(Cst)
    corr = SCFC(SCst,FC,not_sinks)
    print(corr)
    if corr>.9:
        break

SCmat, FCmat = np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float)
for i in range(60):
    for j in range(20):
        SCmat[i,j], FCmat[i,j] = SCst[(j,59-i)], FC[(j,59-i)]

fig, ax = plt.subplots(1,4,figsize=(9,6), tight_layout=True)
ax[0].imshow(vP[3], cmap='Greens')
div = make_axes_locatable(ax)
ax[1].imshow(SCmat, cmap='gist_earth')
ax[2].imshow(FCmat, cmap='rainbow')
ax[3].imshow(DSCHG[3], cmap='rainbow')
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[3].set_xticks([])
ax[3].set_yticks([])
ax[0].set_title(' (a)', loc='left', y=1.01)
ax[1].set_title(' (b)', loc='left', y=1.01)
ax[2].set_title(' (c)', loc='left', y=1.01)
ax[3].set_title(' (d)', loc='left', y=1.01)
ticks1, ticks2, ticks3, ticks4 = (vP[3].min(), vP[3].max()), (SCmat.min(), SCmat.max()), (FCmat.min(), FCmat.max()), (DSCHG[3].min(), DSCHG[3].max())
ticks1, ticks2, ticks3, ticks4 = [ticks1[0], (ticks1[1]+ticks1[0])/2, ticks1[1]], [ticks2[0], (ticks2[1]+ticks2[0])/2, ticks2[1]], [ticks3[0], (ticks3[1]+ticks3[0])/2, ticks3[1]], [0, (ticks4[1]+ticks4[0])/2, ticks4[1]]
fig.colorbar(plt.cm.ScalarMappable(cmap='Greens', norm=colors.Normalize(vmax=ticks1[-1], vmin=ticks1[0])), orientation='horizontal', ax=ax[0], fraction=.02, pad=.02, shrink=.82)
fig.colorbar(plt.cm.ScalarMappable(cmap='gist_earth', norm=colors.Normalize(vmax=ticks2[-1], vmin=ticks2[0])), orientation='horizontal', ax=ax[1], fraction=.02, pad=.02, shrink=.82)
fig.colorbar(plt.cm.ScalarMappable(cmap='rainbow', norm=colors.Normalize(vmax=ticks3[-1], vmin=ticks3[0])), orientation='horizontal', ax=ax[2], fraction=.02, pad=.02, shrink=.82)
fig.colorbar(plt.cm.ScalarMappable(cmap='rainbow', norm=colors.Normalize(vmax=ticks4[-1], vmin=ticks4[0])), orientation='horizontal', ax=ax[3], fraction=.02, pad=.02, shrink=.82)
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[0], fraction=.02, pad=.02, shrink=.82)
cbar.ax.set_title('')
plt.subplot_tool()
plt.show()


import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


S, Cst, current, branches, active, Size = init(Gst)
SC = SC_path(Gst, Cst)

for step in range(nb_steps+1):
    S, Cst, current, branches, active, _, _ = sandpilest(Gst, VP[3], S, Cst, current, branches, active, [], [], not_sinks)



SCFCavg, FCFCavg = np.array([np.mean(SCFCt[:,t]) for t in range(100)]), np.array([np.mean(FCFCt[:,t]) for t in range(100)])
SCFCstd, FCFCstd = np.array([np.std(SCFCt[:,t]) for t in range(100)]), np.array([np.std(FCFCt[:,t]) for t in range(100)])


plt.fill_between(range(1,101), FCFCavg-SCFCstd, FCFCavg+SCFCstd, label='FC-FC*: converging towards MAHLERAN', alpha=.7)
plt.fill_between(range(1,101), SCFCavg-SCFCstd, SCFCavg+SCFCstd, label='SC-FC: predicting of dynamics from topology', alpha=.7)
plt.xticks([1,20,40,60,80,100], labels=['100','2000','4000','6000','8000','10000'])
plt.plot(range(1,101), FCFCavg, c='cyan', linewidth=2)
plt.plot(range(1,101), SCFCavg, c='bisque', linewidth=2)
plt.xlabel('step')
plt.ylabel('correlation')
plt.legend(loc='lower right')
plt.show()







# MS1: unscaled // MS2: scaled to p2
MS1 = np.array([[398,344,210,202],[478,380,230,227],[440,369,180,177],[495,408,256,147]])
MS2 = np.array([[399,351,300,454],[427,379,342,472],[373,377,263,364],[444,403,412,282]])

MS1 = np.array([[-.40,-.34,-.20],[-.48,-.38,-.23],[-.49,-.41,-.15]])
MS2 = np.array([[-.40,-.35,-.45],[-.43,-.38,-.47],[-.44,-.40,-.28]])

plt.imshow(MS1, cmap='rainbow')
for i in range(3):
    for j in range(3):
        _ = plt.text(j, i, MS1[i,j], ha='center', va='center', color='w', fontsize=24, weight='bold')

#plt.title('size slope, unscaled vegetation cover')
plt.xlabel('$V_1 \\qquad\\quad V_2 \\qquad\\quad V_3$', fontsize=22)
plt.ylabel('$E_3 \\qquad\\quad E_2 \\qquad\\quad E_1$', fontsize=22)
plt.xticks([])
plt.yticks([])
plt.show()

figure()

plt.imshow(MS2)
for i in range(4):
    for j in range(4):
        _ = plt.text(j, i, MS2[i,j].astype(int), ha='center', va='center', color='w', fontsize=24, weight='bold')

plt.title('size slope, scaled vegetation cover')
plt.xlabel('$V_1 \\qquad\\qquad\\qquad V_2 \\qquad\\qquad\\qquad V_3 \\qquad\\qquad\\qquad V_4$')
plt.ylabel('$E_4 \\qquad\\qquad\\qquad E_3 \\qquad\\qquad\\qquad E_2 \\qquad\\qquad\\qquad E_1$')
plt.xticks([])
plt.yticks([])
plt.show()


plt.title('Durations')
plt.figure()
plt.imshow(MS)
plt.title('Sizes')
plt.show()

for i in G:
    VP[3][i] *= sum(VP[1].values())/sum(VP[3].values())

sum(VP[0].values())

sum(i>0 for i in VP[0].values())
sum(i>0 for i in VP[1].values())
sum(i>0 for i in VP[2].values())
sum(i>0 for i in VP[3].values())



nx.draw(G, dict(zip(G,G)), node_size=80, node_color=list(VP[2].values()), alpha=.8, edgelist=[])
plt.show()


### SEDIGRAPHS ANALYSIS

A = pd.ExcelFile('experimental_data_hydrograph/plot_event_flow_data.xls')
sed_ind = [i for i in range(39) if A.parse(sheet_name=i).shape[1]>16]
for i in sed_ind:
    X = A.parse(sheet_name=i).values
    i, A.parse(sheet_name=i).values[0]

for k in range(len(X)):
    if type(X[k,0])!=float:
        k,X[k,0]
        break

# detecting oddly shaped sheets
events = []
for i in sed_ind:
    X = A.parse(sheet_name=i).values
    for k in range(len(X)):
        if type(X[k,0])!=float:
            i,k,X[k,0]
            events.append(X[k,0])
            break

datapoints = [i for i in range(len(X)) if type(X[i,1])==datatime.time]
np.where(type(X[:,1])==time)

sed_ind
X = A.parse(sheet_name=i).values
header, time = X[0], X[1:,1]
# extract-out last_t
for t in range(len(time)):
    if type(time[t])==float:
        if np.isnan(time[t]):
            last_t=t
            break

t0, tf = 60*time[0].hour + time[0].minute, 60*time[last_t-1].hour + time[last_t-1].minute
T = tf-t0+1
time_adjusted = list(range(T))
rain = X[1:T+1,np.where(header=='Rain (mm/hr)')[0][0]]
rain_avg = X[1:T+1,np.where(header=='Rain (mm/hr)')[0][0]+1]
Q = X[1:T+1,np.where(header=='original Q')[0][0]]
sedc = X[1:T+1,np.where(header=='Susp.  sed. (g/L)')[0][0]]
for i in np.where(np.isnan(sedc.astype(np.float64)))[0]:
    sedc[i] = 0

plt.plot(time_adjusted, sedc*Q/2, linewidth=2, c='tab:orange', markersize=4)
plt.plot(time_adjusted, rain, linewidth=1)
plt.plot(time_adjusted, rain_avg, linestyle='--', linewidth=2, c='grey')
plt.show()

A = pd.ExcelFile('experimental_data_hydrograph/only_sed.xls')
nb_events = 21
data_tables = np.zeros(nb_events, dtype=np.ndarray)
for i in range(nb_events):
    data_tables[i] = A.parse(sheet_name=i).values

timeseries = pickle.load(open('experimental_data_hydrograph/timeseries.txt','rb'))
all_rain, all_sed = timeseries[0], timeseries[1]

A0 = pd.ExcelFile('experimental_data_hydrograph/plot_event_flow_data.xls')

selected_events = [['29.08.06','07.09.06'],['07.09.05','15.08.06','29.08.06'],['31.07.06','01.08.06','11.08.06','29.08.06','07.09.06'],['05.07.06',

names0 = np.array(A0.sheet_names)
names1 = np.array([data_tables[i][0,0] for i in range(nb_events)])
P4_loc0 = [np.where(names0==i)[0][0] for i in names1 if i[:2]=='P4']
P4_loc1 = [i for i in range(nb_events) if data_tables[i][0,0][:2]=='P4']
names0[P4_loc0]
folder_names = ['P4/plot 4 '+i[3:5]+i[6:8]+i[9:] for i in names1[P4_loc1]]
folder_names[8] += ' calibrated'
folder_inds = P4_loc1

# triple comparison: sandpile vs experimental vs mahleran time-rescaled
for i in range(len(P4_loc1)):
    k = P4_loc1[i]
    fig, ax = plt.subplots(1,3,figsize=(14,5), tight_layout=True)
    ax[0].plot(range(len(all_rain[k])), total_rain[k], label='dropped particles')
    ax[0].plot(range(len(all_rain[k])), total_exits[k], label='exiting particles')
    ax[0].legend()
    ax[0].set_title('simulated')
    ax[0].set_ylabel('rainfall')
    ax[0].set_xticks([])
    ax[0].legend(title='rainfall (mm.h$^{-1}$)', loc='upper left')
    ax[1].plot(range(len(all_rain[k])), all_rain[k], label='rainfall amount')
    ax[1].plot(range(len(all_rain[k])), all_sed[k], label='sediment flow')
    
    ax[1].set_xlabel('minute')
    ax[1].legend(loc='upper left')
    
    event = pd.read_excel('experimental_data_hydrograph/only_sed.xls', sheet_name='EVENT_'+str(k))
    sedtr, hydro = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/'+folder_names[i]+'/sedtr000.dat', float), np.loadtxt('experimental_data_hydrograph/plot_4_outputs/'+folder_names[i]+'/hydro000.dat', float)
    sed = sedtr[1:,1][:np.where(sedtr[1:,1]==0)[0][0]]
    hyd = hydro[1:len(sed)+1,1]
    sedM, hydM = [np.mean(sed[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)], [np.mean(hyd[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)]
    
    ax[2].plot(range(len(hydM)), hydM, label='hydro001 input')
    ax[2].plot(range(len(sedM)), sedM, label='sediment discharge')
    ax[2].legend()
    ax[2].set_title('time rescaled MAHLERAN')
    fig.suptitle('EVENT_'+str(k1)+': '+str(names1[k]))

plt.show()



test1 = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot 4 290706/sedtr001.dat', dtype=float)[:2832,1]
test2 = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot 4 290706/hz = 0.004/sedtr001.dat', dtype=float)[:2782,1]

plt.plot(range(len(test1)), test1)
plt.plot(range(len(test2)), test2)
plt.show()




all_rain, all_sed = np.zeros(nb_events, np.ndarray), np.zeros(nb_events, np.ndarray)
#for i in range(len(data_tables)):
for i in names1+names4:
    X = deepcopy(data_tables[i])
    header, time = X[0], X[1:,1]
    for t in range(len(time)):
        if type(time[t])==float:
            if np.isnan(time[t]):
                last_t=t
                break
    rain = X[1:last_t+1,np.where(header=='Rain (mm/hr)')[0][0]]
    rain_avg = X[1:last_t+1,np.where(header=='Rain (mm/hr)')[0][0]+1]	
    Q = X[1:last_t+1,np.where(header=='scaled original Q')[0][0]]
    sedc = X[1:last_t+1,np.where(header=='Susp.  sed. (g/L)')[0][0]]
    for k in range(len(sedc)):
        if type(sedc[k])==float:
            if np.isnan(sedc[k]):
                sedc[k]=0
        else:
            sedc[k]=0
    #total_rain[i], total_sed[i] = sum(rain), sum(sedc*Q)
    all_rain[i], all_sed[i] = rain, sedc*Q
    
    fig, ax = plt.subplots(2,1,figsize=(10,7), tight_layout=True)
    ax[0].plot(list(range(last_t)), rain, linewidth=1, label='instant')
    ax[0].plot(list(range(last_t)), rain_avg, c='grey', linestyle='--', linewidth=2, label='5-min intensity')
    ax[0].set_xticks([])
    ax[0].legend(title='rainfall (mm.h$^{-1}$)', loc='upper left')
    ax[1].plot(list(range(last_t)), sedc*Q, linewidth=2, c='tab:orange', markersize=4, label='sediment flow (l$^{-1}$)')
    ax[1].set_xlabel('minute')
    ax[1].legend(loc='upper left')
    plt.suptitle('Event n°'+str(i)+': '+str(header[0]))
    plt.savefig('experimental_data_hydrograph/sedigraphs/EVENT_'+str(i)+'.jpg')
    plt.close('all')


plt.plot(range(nb_events), total_rain, marker='o', label='total rain (mm.h$^{-1}$)', markersize=6)
plt.plot(range(nb_events), total_sed, marker='o', label='total sediment flow (l$^{-1}$)', markersize=6)
plt.xlabel('event index')
plt.xticks(list(range(nb_events))[::2])
plt.legend()
plt.tight_layout()
plt.show()


# computing scaled original Q
missingQscaled = [0,1,3,10,14,15,16,20]
Qtot = [5275.218, 2820.3, 2294.442, 6549.66, 1503.864, 452.292, 803.676, 6366.168]
ratio = 2181.18/np.array(Qtot)

# checking for wells
sum([G.degree(i)==0 for i in not_sinks])
nopbc = [(i,j) for i,j in Gst.edges if (i[0],j[0]) not in [(0,19),(19,0)]]

nx.draw(Gst, dict(zip(G,G)), node_size=0, with_labels=True, edgelist=nopbc)

tracker, ongoing = [], False
S, C, current, branches, active, Dst, Sst = zeros(Cst)
for step in range(nb_steps+1):
    saved_branches = branches
    if len(branches):
        ongoing = True
    S, Cst, current, branches, active, Dst, Sst = sandpilest(Gst, VP[3], S, Cst, current, branches, active, Dst, Sst, not_sinks)
    if len(branches) == 0 and ongoing:
        tracker.append(step-1)# avalanche stopped in previous step
        ongoing = False
        step, Sst[-1], list(saved_branches.values())
    if len(branches)==1:
        branches
        break
        tracker.append(step)

# with stochastic descent
Ep = eprob(Geff, EP[3])
tracker, ongoing = [], False
S, C, current, branches, active, Seff, Deff = zeros(Ceff)
for step in range(nb_steps+1):
    saved_ava = current
    if len(branches):
        ongoing = True
    S, current, branches, active, Seff, Deff = simple_sandpile(Geff, VP[3], Ep, S, current, branches, active, Seff, Deff, not_sinks)
    if len(current) == 0 and ongoing:
        tracker.append(step-1)# avalanche stopped in previous step
        ongoing = False
        #saved_ava, Deff[-1], Seff[-1]

### 


for k in range(1):#len(all_rain)):
    T_ratio, total_rain, total_exits = int(max(all_rain[k])), np.zeros(len(all_rain[k]), dtype=int), np.zeros(len(all_rain[k]), dtype=int)
    S, C, current, branches, active, Seff, Deff = zeros(Cst)
    for i in S:
        S[i] = s0
    for T in range(len(all_rain[k])):
        nb_rain, exiting, p_rain = 0, 0, all_rain[k][T]/T_ratio
        while nb_rain < all_rain[k][T]:
            S, current, branches, active, Seff, Deff, nb_rain, exiting = simple_sandpilest_sedigraph_with_rain_intensity(Gst, VP[3], S, current, branches, active, Seff, Deff, nb_rain, exiting, p_rain, not_sinks)
            #S, current, branches, active, Seff, Deff, nb_rain, exiting = simple_sandpile_sedigraph_with_rain_intensity(Geff, VP[3], ep, S, current, branches, active, Seff, Deff, nb_rain, exiting, p_rain, not_sinks)
        total_rain[T] = nb_rain
        total_exits[T] = exiting
    fig, ax = plt.subplots(1, 2, figsize=(12,5), tight_layout=True)
    ax[0].plot(range(len(all_rain[k])), total_rain, label='dropped particles')
    ax[0].plot(range(len(all_rain[k])), total_exits, label='exiting particles')
    ax[0].legend()
    ax[0].set_title('simulated')
    ax[1].plot(range(len(all_rain[k])), all_rain[k], label='rainfall amount')
    ax[1].plot(range(len(all_rain[k])), all_sed[k], label='sediment flow')
    ax[1].legend()
    ax[1].set_title('experimental')
    fig.suptitle('EVENT_'+str(k+1)+', initialized at '+str(s0))
    plt.savefig('simulated_EVENT_'+str(k+1)+'_init_'+str(s0)+'.png')
    plt.close('all')

plt.plot(range(len(all_rain[k])), total_rain)
plt.plot(range(len(all_rain[k])), total_exits)
plt.figure()
plt.plot(range(len(all_rain[k])), all_rain[k])
plt.plot(range(len(all_rain[k])), all_sed[k])
plt.show()


## without stochastic rain

total_rain, total_exits = np.zeros(21,np.ndarray), np.zeros(21,np.ndarray)
for i in range(len(P4_loc1)):
    k = P4_loc1[i]
    T_ratio, total_rain[k], total_exits[k] = int(max(all_rain[k])), np.zeros(len(all_rain[k]), dtype=int), np.zeros(len(all_rain[k]), dtype=int)
    S, C, current, branches, active, Seff, Deff = zeros(Cst)
    for T in range(len(all_rain[k])):
        nb_rain, exiting = 0, 0
        while nb_rain < all_rain[k][T]:
            #S, current, branches, active, Seff, Deff, nb_rain, exiting = simple_sandpilest_sedigraph(Gst, VP[3], S, current, branches, active, Seff, Deff, nb_rain, exiting, not_sinks)
            S, current, branches, active, Seff, Deff, nb_rain, exiting = simple_sandpilest_sedigraph_probabilistic_toppling(Gst, VP[3], S, current, branches, active, Seff, Deff, nb_rain, exiting, not_sinks)
        total_rain[k][T] = nb_rain
        total_exits[k][T] = exiting
    
    sedtr, hydro = np.loadtxt('experimental_data_hydrograph/mahleran_outputs/'+folder_names[i]+'/sedtr000.dat', float), np.loadtxt('experimental_data_hydrograph/mahleran_outputs/'+folder_names[i]+'/hydro000.dat', float)
    sed = sedtr[1:,1][:np.where(sedtr[1:,1]==0)[0][0]]
    hyd = hydro[1:len(sed)+1,1]
    sedM, hydM = np.array([np.mean(sed[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)]), np.array([np.mean(hyd[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)])
    slen = len(all_rain[k])
    
    fig = plt.figure(figsize=(8,6), tight_layout=True)
    plt.plot(range(slen), all_rain[k], linestyle='--', linewidth=1, label='rainfall amount')
    plt.plot(range(slen), all_sed[k], linestyle='--', linewidth=2, c='tab:green', label='experiment')
    plt.plot(range(slen), total_exits[k], c='tab:orange', linewidth=2, alpha=.5, label='sandpile')
    plt.plot(range(slen), 1e3*sedM[:slen], c='tab:red', linewidth=2, alpha=.5, label='Mahleran (*1e3)')
    plt.legend()
    plt.title(str(names0[P4_loc0][i]))


    plt.savefig('Tricomparison on '+str(names0[P4_loc0][i])+', with_storage.png')
    plt.close('all')

plt.show()

slen = len(all_rain[k])

fig = plt.figure(figsize=(8,6), tight_layout=True)
plt.plot(range(slen), all_rain[k], linestyle='--', label='rainfall amount')
plt.plot(range(slen), all_sed[k], c='tab:green', label='observed sediment')
plt.plot(range(slen), total_exits[k], c='tab:orange', label='sandpile sediment')
plt.legend()
plt.show()

plt.plot(range(slen), hydM[:slen], c='tab:blue', label='Mahleran rain')
plt.plot(range(slen), 1e4*sedM[:slen], c='tab:red', label='Mahleran sediment')
plt.legend()
plt.show()


event = pd.read_excel('experimental_data_hydrograph/only_sed.xls', sheet_name='EVENT_'+k1)
tr = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/'+folder_names[k]+'/sedtr001.dat', dtype=float)
hydro = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/'+folder_names[k]+'/hz = 0.004/hydro001.dat', dtype=float)



tr = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot_4_010806/hz = 0.004/sedtr001.dat', dtype=float)
disch = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot_4_010806/hz = 0.004/seddisch001.dat', dtype=float)
event = pd.read_excel('experimental_data_hydrograph/plot_event_flow_data.xls', sheet_name='P4 01.08.06')
event = pd.read_excel('experimental_data_hydrograph/only_sed.xls', sheet_name='EVENT_17')
timeseries = pickle.load(open('experimental_data_hydrograph/timeseries.txt','rb'))
all_rain, all_sed = timeseries[0], timeseries[1]

rain, sed = all_rain[17], all_sed[17]
plt.plot(range(107), rain)
plt.plot(range(107), sed)
plt.show()

sedM = tr[1:7113,1]
plt.plot(range(len(sedM)), sedM)
plt.show()

hydro = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot_4_010806/hz = 0.004/hydro001.dat', dtype=float)
#dschg = np.loadtxt('experimental_data_hydrograph/plot_4_outputs/plot_4_010806/hz = 0.004/dschg001.dat', dtype=float)
hydroM = hydro[1:7113,1]
plt.plot(range(len(hydroM)), hydroM, label='hydro001 input')
plt.plot(range(len(sedM)), sedM, label='sediment discharge')
plt.xlabel('time')
plt.title('EVENT_18, P4: 01.08.06, MAHLERAN')
plt.legend()
plt.show()

sedM2 = [np.mean(sedM[t*60:(t+1)*60]) for t in range(120)]
hydroM2 = [np.mean(hydroM[t*60:(t+1)*60]) for t in range(120)]
plt.plot(range(len(hydroM2)), hydroM2, label='hydro001 input')
plt.plot(range(len(sedM2)), sedM2, label='sediment discharge')
plt.xlabel('time')
plt.title('EVENT_18, P4: 01.08.06, MAHLERAN')
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 2, figsize=(11,5), tight_layout=True)
ax[0].plot(range(len(hydroM2)), hydroM2, label='hydro001 input')
ax[0].plot(range(len(sedM2)), sedM2, label='sediment discharge')
ax[0].legend()
ax[0].set_title('time rescaled MAHLERAN')
ax[1].plot(range(len(all_rain[k])), all_rain[k], label='rainfall amount')
ax[1].plot(range(len(all_rain[k])), all_sed[k], label='sediment flow')
ax[1].legend()
ax[1].set_title('experimental')
fig.suptitle('EVENT_18, P4: 01.08.06')
plt.show()

ax[1].plot(range(len(all_rain[k])), all_rain[k], label='rainfall amount')
ax[1].plot(range(len(all_rain[k])), all_sed[k], label='sediment flow')
ax[1].legend()
ax[1].set_title('experimental')
plt.show()

### COMPILING VCOVER=10 AND VCOVER=20 SIZE BEHAVIOR

grid, vc = 'S', '10'

S10_Sizes = np.vstack(tuple(pickle.load(open('S/vc_10/pc_'+pstr+'/Size_eff.txt','rb')) for pstr in ['2','5','8']))
S40_Sizes = np.vstack(tuple(pickle.load(open('S/vc_40/pc_'+pstr+'/Size_eff.txt','rb')) for pstr in ['2','5','8']))
GS10_Sizes = np.vstack(tuple(pickle.load(open('GS/vc_10/pc_'+pstr+'/Size_eff.txt','rb')) for pstr in ['2','5','8']))
GS40_Sizes = np.vstack(tuple(pickle.load(open('GS/vc_40/pc_'+pstr+'/Size_eff.txt','rb')) for pstr in ['2','5','8']))

grid, vc = 'GS', '40'
GS40_Corrs, GS40_exits = np.vstack(tuple(pickle.load(open(grid+'/vc_'+vc+'/pc_'+pstr+'/Cr_eff.txt','rb')) for pstr in ['2','5','8'])), np.vstack(tuple(pickle.load(open(grid+'/vc_'+vc+'/pc_'+pstr+'/exits_eff.txt','rb')) for pstr in ['2','5','8']))
with open(grid+'/vc_'+vc+'/Cr_eff.txt','wb') as ll:
    pickle.dump(GS40_Corrs, ll)

with open(grid+'/vc_'+vc+'/exits_eff.txt','wb') as ll:
    pickle.dump(GS40_exits, ll)


S10_Sizes, GS10_Sizes, S40_Sizes, GS40_Sizes = pickle.load(open('S/vc_10/Size_eff.txt','rb')), pickle.load(open('GS/vc_10/Size_eff.txt','rb')), pickle.load(open('S/vc_40/Size_eff.txt','rb')), pickle.load(open('GS/vc_40/Size_eff.txt','rb'))
S10_Corrs1, GS10_Corrs1, S40_Corrs1, GS40_Corrs1 = pickle.load(open('S/vc_10/Cr_len.txt','rb')), pickle.load(open('GS/vc_10/Cr_len.txt','rb')), pickle.load(open('S/vc_40/Cr_len.txt','rb')), pickle.load(open('GS/vc_40/Cr_len.txt','rb'))
S10_Corrs3, GS10_Corrs3, S40_Corrs3, GS40_Corrs3 = pickle.load(open('S/vc_10/Cr_var.txt','rb')), pickle.load(open('GS/vc_10/Cr_var.txt','rb')), pickle.load(open('S/vc_40/Cr_var.txt','rb')), pickle.load(open('GS/vc_40/Cr_var.txt','rb'))
S10_exits, GS10_exits, S40_exits, GS40_exits = pickle.load(open('S/vc_10/exits_eff.txt','rb')), pickle.load(open('GS/vc_10/exits_eff.txt','rb')), pickle.load(open('S/vc_40/exits_eff.txt','rb')), pickle.load(open('GS/vc_40/exits_eff.txt','rb'))
S10_alpha, GS10_alpha, S40_alpha, GS40_alpha = pickle.load(open('S/vc_10/shape.txt','rb')), pickle.load(open('GS/vc_10/shape.txt','rb')), pickle.load(open('S/vc_40/shape.txt','rb')), pickle.load(open('GS/vc_40/shape.txt','rb'))
S10_beta, GS10_beta, S40_beta, GS40_beta = pickle.load(open('S/vc_10/scale.txt','rb')), pickle.load(open('GS/vc_10/scale.txt','rb')), pickle.load(open('S/vc_40/scale.txt','rb')), pickle.load(open('GS/vc_40/scale.txt','rb'))

grid, vc = 'S', '10'
Corrs, exits, betas = S10_Corrs1, S10_exits, S10_alpha
fig, ax = plt.subplots(1,3,figsize=(14,4),tight_layout=True)
xmin, xmax, ymin, ymax = min([min([np.mean(exits[p][t]) for t in range(100)]) for p in range(3)])-200, max([max([np.mean(exits[p][t]) for t in range(100)]) for p in range(3)])+200, min([min(Corrs[p]) for p in range(3)])-.01, max([max(Corrs[p]) for p in range(3)])+.01
_ = [ax[p].scatter([np.mean(exits[p][t]) for t in range(100)], Corrs[p], s=20) for p in range(3)]
_ = [ax[p].set_title('$\\rho$='+str(round(pearsonr([np.mean(exits[p][t]) for t in range(100)], Corrs[p])[0],2))) for p in range(3)]
_ = [ax[p].set_xlim(xmin,xmax) for p in range(3)]
_ = [ax[p].set_ylim(ymin,ymax) for p in range(3)]
_ = [ax[p].set_yticklabels([]) for p in range(1,3)]
_ = [ax[p].set_xlabel('average number of exiting particles', size=14) for p in range(3)]
_ = ax[0].set_ylabel('SC-FC correlation', size=14)
plt.suptitle(grid+', $v_c$='+vc)
#plt.show()

fig, ax = plt.subplots(1,3,figsize=(14,4),tight_layout=True)
xmin, xmax, ymin, ymax = min([min([np.mean(exits[p][t]) for t in range(100)]) for p in range(3)])-200, max([max([np.mean(exits[p][t]) for t in range(100)]) for p in range(3)])+200, min([min(betas[p]) for p in range(3)])-.01, max([max(betas[p]) for p in range(3)])+.01
_ = [ax[p].scatter([np.mean(exits[p][t]) for t in range(100)], betas[p], s=20) for p in range(3)]
_ = [ax[p].set_title('$\\rho$='+str(round(pearsonr([np.mean(exits[p][t]) for t in range(100)], betas[p])[0],2))) for p in range(3)]
_ = [ax[p].set_xlim(xmin,xmax) for p in range(3)]
_ = [ax[p].set_ylim(ymin,ymax) for p in range(3)]
_ = [ax[p].set_yticklabels([]) for p in range(1,3)]
_ = [ax[p].set_xlabel('average number of exiting particles', size=14) for p in range(3)]
_ = ax[0].set_ylabel('scale parameter', size=14)
plt.suptitle(grid+', $v_c$='+vc)
plt.show()

fig, ax = plt.subplots(1,3,figsize=(14,4),tight_layout=True)
xmin, xmax, ymin, ymax = min([min(Corrs[p]) for p in range(3)]), max([max(Corrs[p]) for p in range(3)]), min([min(betas[p]) for p in range(3)]), max([max(betas[p]) for p in range(3)])
_ = [ax[p].scatter(Corrs[p], betas[p], s=20) for p in range(3)]
_ = [ax[p].set_title('$\\rho$='+str(round(pearsonr(Corrs[p], betas[p])[0],2))) for p in range(3)]
_ = [ax[p].set_xlim(xmin,xmax) for p in range(3)]
_ = [ax[p].set_ylim(ymin,ymax) for p in range(3)]
_ = [ax[p].set_yticklabels([]) for p in range(1,3)]
plt.show()

## SLOPE CALCULATION
shape, scale, grid, vc = np.zeros((3,100), float), np.zeros((3,100), float), 'S', '40'
Sizes = pickle.load(open(grid+'/vc_'+vc+'/Size_eff.txt','rb'))
for p in range(3):
    for t in range(100):
        print('TRIAL N°',t)
        wb = Fit_Weibull_3P(sum(Sizes[p][t],[]))
        shape[p,t], scale[p,t] = wb.alpha, wb.beta

with open(grid+'/vc_'+vc+'/shape.txt','wb') as ll:
    pickle.dump(shape,ll)

with open(grid+'/vc_'+vc+'/scale.txt','wb') as ll:
    pickle.dump(scale,ll)


### GENERATE CR_EFF (SCw, FClen)

Corrs1, Corrs3, grid, vc = np.zeros((3,100), float), np.zeros((3,100), float), 'S', '40'
for p in range(3):
    pstr = str(p*3+2)
    SCw, FC1eff, FC3eff = pickle.load(open(grid+'/vc_'+vc+'/pc_'+pstr+'/SCw.txt','rb')), pickle.load(open(grid+'/vc_'+vc+'/pc_'+pstr+'/FC1eff.txt','rb')), pickle.load(open(grid+'/vc_'+vc+'/pc_'+pstr+'/FC3eff.txt','rb'))
    for t in range(100):
        Corrs1[p,t], Corrs3[p,t] = SCFC(SCw[t], FC1eff[t], not_sinks), SCFC(SCw[t], FC3eff[t], not_sinks)

with open(grid+'/vc_'+vc+'/Cr_len.txt','wb') as ll:
    pickle.dump(Corrs1, ll)

with open(grid+'/vc_'+vc+'/Cr_var.txt','wb') as ll:
    pickle.dump(Corrs3, ll)

## AVERAGES


Size_GS10, Size_S10 = pickle.load(open('old/GS/vc_10/Size_eff.txt','rb')), pickle.load(open('old/S/vc_10/Size_eff.txt','rb'))
Size_GS40, Size_S40 = pickle.load(open('old/GS/vc_40/Size_eff.txt','rb')), pickle.load(open('old/S/vc_40/Size_eff.txt','rb'))

Size_GS.shape

avg_GS10, avg_S10, avg_GS40, avg_S40 = np.zeros((3,100), float), np.zeros((3,100), float), np.zeros((3,100), float), np.zeros((3,100), float)
for p in range(3):
    p
    for t in range(100):
        avg_GS10[p,t], avg_S10[p,t], avg_GS40[p,t], avg_S40[p,t] = np.mean([np.mean(Size_GS10[p,t][run]) for run in range(100)]), np.mean([np.mean(Size_S10[p,t][run]) for run in range(100)]), np.mean([np.mean(Size_GS40[p,t][run]) for run in range(100)]), np.mean([np.mean(Size_S40[p,t][run]) for run in range(100)])

fig, ax = plt.subplots(1,2,figsize=(10,4),tight_layout=True)
for p in range(3):
    ax[0].violinplot([avg_GS10[p] for p in range(3)], [.2,.5,.8])
    ax[1].violinplot([avg_S10[p] for p in range(3)], [.2,.5,.8])

ax[0].plot([.2,.5,.8], [np.mean(avg_GS10[p]) for p in range(3)], marker='o')
ax[1].plot([.2,.5,.8], [np.mean(avg_S10[p]) for p in range(3)], marker='o')
_ = [ax[i].set_ylim(20,35) for i in range(3)]
ax[0].set_title('grassland vs=10%')
ax[1].set_title('shrubland vs=10%')
#plt.show()

fig, ax = plt.subplots(1,2,figsize=(10,4),tight_layout=True)
for p in range(3):
    ax[0].violinplot([avg_GS40[p] for p in range(3)], [.2,.5,.8])
    ax[1].violinplot([avg_S40[p] for p in range(3)], [.2,.5,.8])

ax[0].plot([.2,.5,.8], [np.mean(avg_GS40[p]) for p in range(3)], marker='o')
ax[1].plot([.2,.5,.8], [np.mean(avg_S40[p]) for p in range(3)], marker='o')
_ = [ax[i].set_ylim(8,16) for i in range(3)]
ax[0].set_title('grassland vs=40%')
ax[1].set_title('shrubland vs=40%')
plt.show()



### REPRESENTING VEG DISTRIBUTIONS AS SCATTERPLOT (PROBABILITY DENSITY FUNCTION)
veg1 = vP[0][np.where(vP[0])]/100
veg4 = vP[3][np.where(vP[3])]/100

interval = np.array([round(i*.1,1) for i in range(12)])
nb1, _ = np.histogram(veg1, bins=interval, density=True)
nb4, _ = np.histogram(veg4, bins=interval, density=True)
nb1 *= 10
nb4 *= 10

fig = plt.figure(figsize=(5,3.5), constrained_layout=True)
plt.plot(interval[:-1], nb1, marker='o', markersize=8, color='green', alpha=.8, label='grassland')
plt.plot(interval[:-1], nb4, marker='o', markersize=8, color='darkorange', alpha=.8, label='shrubland')
plt.xlabel('vegetation density')
plt.ylabel('count (%): probability density function')
plt.legend()
plt.xticks(interval[:-1])
plt.yticks([6,10,14,18,22])
plt.show()

### DEGREE REPRESENTATION
d0 = [Geff[0].out_degree(i) for i in G]
d3 = [Geff[3].out_degree(i) for i in G]

b0, c0 = np.unique(d0, return_counts=True)
b3, c3 = np.unique(d3, return_counts=True)

fig = plt.figure(constrained_layout=True)
plt.scatter(b0, c0, color='green', marker='x', label='grassland', alpha=.8, s=40, linewidth=3)
plt.scatter(b3, c3, color='orange', marker='x', label='shrubland', alpha=.8, s=40, linewidth=3)
plt.yscale('log')
fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('out-degree')
plt.ylabel('count')
plt.legend(loc='lower right')
plt.yticks([10,100,1000])
plt.show()

from matplotlib.ticker import MaxNLocator


plt.scatter(b0, c0, color='green')
plt.scatter(b3, c3, color='orange')
plt.yscale('log')
plt.show()

deg_distrib = np.zeros((2,3,100), list)
### IMPORT THIS DATA FROM OLD LAPTOP ONCE AT HOME
Vs_grass, Es_grass = pickle.load(open('grids/grassland/Vs.txt','rb')), pickle.load(open('grids/grassland/Es.txt','rb'))
Vs_shrub, Es_shrub = pickle.load(open('grids/shrubland/Vs.txt','rb')), pickle.load(open('grids/shrubland/Es.txt','rb'))

Es_grass, Es_shrub = pickle.load(open('grids/grassland/Egrass.txt','rb')), pickle.load(open('grids/shrubland/Eshrub.txt','rb'))
grass_len, shrub_len = np.zeros((3,100), float), np.zeros((3,100), float)
counter0, counter3 = np.zeros((3,100), np.ndarray), np.zeros((3,100), np.ndarray)
v_sample = vP[3][np.where(vP[3])]/100
Vs, Es_grass, Es_shrub = np.zeros((3,100), dict), np.zeros((3,100), dict), np.zeros((3,100), dict)
for p in range(3):
    for t in range(100):
        #V = tclustering_square(G, .3, p*.3+.2, random.choice(v_sample))
        #Vs[p,t] = V
        Es_grass[p,t] = generate_elevation_new(G, not_sinks, V, EP[3], alpha_grass, alpha_std_grass, beta_grass)[0]
        Es_shrub[p,t] = generate_elevation_new(G, not_sinks, V, EP[3], alpha_shrub, alpha_std_shrub, beta_shrub)[0]
        G_grass, G_shrub = eff_lattice(G, Es_grass[p,t]), eff_lattice(G, Es_shrub[p,t])
        deg_distrib[0,p,t], deg_distrib[1,p,t] = [G_grass.out_degree(i) for i in G], [G_shrub.out_degree(i) for i in G]
        pp0, pp3 = list(nx.all_pairs_shortest_path_length(G_grass)), list(nx.all_pairs_shortest_path_length(G_shrub))
        grass_len, shrub_len = np.mean([max(pp0[k][1].values()) for k in range(N)]), np.mean([max(pp3[k][1].values()) for k in range(N)])
        counter0[p,t], counter3[p,t] = np.unique(deg_distrib[0,p,t], return_counts=True)[1]/N, np.unique(deg_distrib[1,p,t], return_counts=True)[1]/N

with open('data2025/Vs.txt','wb') as f:
    pickle.dump(Vs,f)

with open('data2025/Es_grass.txt','wb') as f:
    pickle.dump(Es_grass,f)

with open('data2025/Es_shrub.txt','wb') as f:
    pickle.dump(Es_shrub,f)

avg0_length, avg3_length = np.zeros((3,100), float), np.zeros((3,100), float)
for p in range(3):
    print('p=',p)
    for t in range(100):
        t
        G_grass, G_shrub = eff_lattice(G, Es_grass[p,t]), eff_lattice(G, Es_shrub[p,t])
        pp0 = list(nx.all_pairs_shortest_path_length(G_grass))
        pp3 = list(nx.all_pairs_shortest_path_length(G_shrub))
        avg0_length[p,t] = np.mean([max(pp0[k][1].values()) for k in range(len(pp0))])
        avg3_length[p,t] = np.mean([max(pp3[k][1].values()) for k in range(len(pp3))])



fig = plt.figure(constrained_layout=True)
plt.scatter([0,1,2,3], [100*np.mean([counter0[p,t][i] for t in range(100) for p in range(3)]) for i in range(4)], color='green', marker='x', label='grassland', alpha=.8, s=40, linewidth=3)
plt.scatter([0,1,2,3], [100*np.mean([counter3[p,t][i] for t in range(100) for p in range(3)]) for i in range(4)], color='orange', marker='x', label='shrubland', alpha=.8, s=40, linewidth=3)
#plt.yscale('log')
fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('out-degree')
plt.ylabel('percentage of nodes (%)')
plt.legend(loc='lower right')
plt.show()

fig = plt.figure(constrained_layout=True)
vp = plt.violinplot([[counter0[p,t][i] for t in range(100) for p in range(3)] for i in range(4)], [0,1,2,3])
plt.scatter([0,1,2,3], [np.mean([counter0[p,t][i] for t in range(100) for p in range(3)]) for i in range(4)], color='green', marker='x', label='grassland', alpha=.8, s=40, linewidth=3)
for v in vp['bodies']:
    v.set_facecolor('green')
    v.set_edgecolor('black')
    v.set_linewidth(.5)
    v.set_alpha(.5)

vp = plt.violinplot([[counter3[p,t][i] for t in range(100) for p in range(3)] for i in range(4)], [0,1,2,3])
plt.scatter([0,1,2,3], [np.mean([counter3[p,t][i] for t in range(100) for p in range(3)]) for i in range(4)], color='orange', marker='x', label='shrubland', alpha=.8, s=40, linewidth=3)
for v in vp['bodies']:
    v.set_facecolor('orange')
    v.set_edgecolor('black')
    v.set_linewidth(.5)
    v.set_alpha(.5)

plt.show()

### Shortest paths to evaluate effective lengths

pp0 = list(nx.all_pairs_shortest_path_length(Geff[0]))
pp3 = list(nx.all_pairs_shortest_path_length(Geff[3]))

pp0_maxes=[max(pp0[k][1].values()) for k in range(len(pp0))]
pp3_maxes=[max(pp3[k][1].values()) for k in range(len(pp3))]

np.mean(pp0_maxes) == 38.5
np.mean(pp3_maxes) == 39

### RECREATING SCFCT

# COMPUTING SC_EDGE

epmap0 = eprob_map(Geff[0], EP[0])
sce0 = SC_edge(Geff[0], epmap0)
plt.imshow(nx_to_mat(sc0), cmap='Spectral_r')
plt.show()

epmap3 = eprob_map(Geff[3], EP[3])
sce3 = SC_edge(Geff[3], epmap3)
plt.imshow(nx_to_mat(sc3), cmap='Spectral_r')
plt.show()

_, ceff0, _, _, _, _ = init(Geff[0])
_, ceff3, _, _, _, _ = init(Geff[3])

scp0 = SCw_path(Geff[0], ceff0, epmap0)
scp3 = SCw_path(Geff[3], ceff3, epmap3)

plt.imshow(nx_to_mat(scp3), cmap='Spectral_r')
plt.show()


nb_runs, nb_steps = 100, int(4*1e5)
eprob0, eprob3 = eprob(Geff[0], EP[0]), eprob(Geff[3], EP[3])
FCt, Sizes = np.zeros((2,100,101), dict), np.zeros((2,100), list)
for run in range(nb_runs):
    run
    S0, C0, current0, branch0, active0, Sizes[0,run] = zeros(ceff0)
    S3, C3, current3, branch3, active3, Sizes[1,run] = zeros(ceff3)
    for t in range(1,nb_steps+1):
        S0, C0, current0, branch0, active0, Sizes[0,run] = sandpile(Geff[0], VP[0], eprob0, S0, C0, current0, branch0, active0, Sizes[0,run], not_sinks)
        S3, C3, current3, branch3, active3, Sizes[1,run] = sandpile(Geff[3], VP[3], eprob3, S3, C3, current3, branch3, active3, Sizes[1,run], not_sinks)
        if t==1000:
            FCt[0,run,0], FCt[1,run,0] = FC_path(C0,'var'), FC_path(C3,'var')
            k = 1
        if t%(4*1e3)==0:
            FCt[0,run,k], FCt[1,run,k] = FC_path(C0,'var'), FC_path(C3,'var')
            k += 1

SCFCt = np.zeros((2,100,101), float)
for run in range(nb_runs):
    for k in range(101):
        SCFCt[0,run,k], SCFCt[1,run,k] = SCFC(sce0,FCt[0,run,k],not_sinks), SCFC(sce3,FCt[1,run,k],not_sinks)

avgSCFC0, avgSCFC1 = np.array([np.mean([SCFCt[0,run,k] for run in range(nb_runs)]) for k in range(101)]), np.array([np.mean([SCFCt[1,run,k] for run in range(nb_runs)]) for k in range(101)])
stdSCFC0, stdSCFC1 = np.array([np.std([SCFCt[0,run,k] for run in range(nb_runs)]) for k in range(101)]), np.array([np.std([SCFCt[1,run,k] for run in range(nb_runs)]) for k in range(101)])

with open('data2025/SCFCt_stochastic_natural_landscapes.txt','wb') as f:
    pickle.dump(SCFCt, f)

sizes0, sizes1 = sum(Sizes[0],[]), sum(Sizes[1],[])
r0, r1 = powerlaw.Fit(sizes0, xmin=2), powerlaw.Fit(sizes1, xmin=2)
max_size = max(max(sizes0), max(sizes1))
binning = np.logspace(0, log10(max_size), 30)
bx0, by0 = log_binning_with_input(Counter(sizes0),binning)
bx0, by0 = bx0[~np.isnan(bx0)], by0[~np.isnan(by0)]
bx1, by1 = log_binning_with_input(Counter(sizes1),binning)
bx1, by1 = bx1[~np.isnan(bx1)], by1[~np.isnan(by1)]

fig, ax = plt.subplots(1, 2, figsize=(8.8,3), constrained_layout=True)
ax[0].fill_between(range(101), avgSCFC0-stdSCFC0, avgSCFC0+stdSCFC0, color='forestgreen', alpha=.8, label='grassland')
ax[0].fill_between(range(101), avgSCFC1-stdSCFC1, avgSCFC1+stdSCFC1, color='darkorange', alpha=.8, label='shrubland')
ax[0].plot(range(101), avgSCFC0, lw=1, c='lime')#or honeydew
ax[0].plot(range(101), avgSCFC1, lw=1, c='bisque')#or yellow
#ax[0].set_ylim(bottom=.1)
ax[0].legend()
ax[0].set_xticks([0]+[25*k for k in range(1,5)], labels=[1]+[100*k for k in range(1,5)])
ax[0].set_xlabel('time step (10$^3$)')
ax[0].set_ylabel('SC-FC correlation')
ax[0].set_title('(a)', loc='left')
r0.truncated_power_law.plot_pdf(color='lightgreen', linestyle='--', ax=ax[1], alpha=.9, label='grassland')
ax[1].scatter(bx0, by0/sum(by0), marker='.', color='forestgreen', alpha=.8, s=60)
r1.truncated_power_law.plot_pdf(color='sandybrown', linestyle='--', ax=ax[1], alpha=.9, label='shrubland')
ax[1].scatter(bx1, by1/sum(by1), marker='.', color='darkorange', alpha=.8, s=60)
ax[1].legend(title='truncated power law fit', alignment='left', loc='lower left')
ax[1].set_xlabel('avalanche size')
ax[1].set_ylabel('frequency')
ax[1].set_title('(b)', loc='left')
plt.show()

fc0_total, fc3_total = {i:0 for i in G}, {i:0 for i in G}
for run in range(100):
    for i in G:
        fc0_total[i] += FCt[0,run,-1][i]
        fc3_total[i] += FCt[1,run,-1][i]

# which fc0 correlates best?
best_fc0_idx, best_fc3_idx = np.argmax(SCFCt[0,:,-1]), np.argmax(SCFCt[1,:,-1])
SCFCt[0,best_fc0_idx,-1], SCFCt[1,best_fc3_idx,-1]
#best_fc0_idx = -9
best_fc0, best_fc3 = FCt[0,best_fc0_idx,-1], FCt[1,best_fc3_idx,-1]

## SC_FC_SC_FC

MT0, MT3 = deepcopy(eP[0]), deepcopy(eP[3])
for i in range(60):
    MT0[i] /= np.mean(MT0[i])
    MT3[i] /= np.mean(MT3[i])


selected_fc0, selected_fc3 = best_fc0, best_fc3
#selected_fc0, selected_fc3 = fc0_total, fc3_total
fig, ax = plt.subplots(1,6, figsize=(8.8,4.5), constrained_layout=True)
#p = ax[0].imshow(nx_to_mat(sce0), cmap='terrain'), ax[1].imshow(nx_to_mat(fc0_total), cmap='terrain'), ax[2].imshow(nx_to_mat(sce3), cmap='terrain'), ax[3].imshow(nx_to_mat(fc3_total), cmap='terrain')
p = ax[0].imshow(vP[0]/100, cmap='Greens'), ax[1].imshow(nx_to_mat(sce0), cmap='terrain_r'), ax[2].imshow(nx_to_mat(selected_fc0), cmap='terrain_r'), ax[3].imshow(vP[3]/100, cmap='Greens'), ax[4].imshow(nx_to_mat(sce3), cmap='terrain_r'), ax[5].imshow(nx_to_mat(selected_fc3), cmap='terrain_r')
#p[0], p[3] = ax[0].imshow(MT0, cmap='terrain_r'), ax[3].imshow(MT3, cmap='terrain_r')
_ = [ax[i].set_title(['(c)','(d)','(e)','(f)','(g)','(h)'][i], loc='left', y=1.0) for i in range(6)]
_ = [ax[i].set_xticks([]) for i in range(6)]
_ = [ax[i].set_yticks([]) for i in range(6)]
ticks1, ticks2, ticks3, ticks4, ticks5, ticks6 = (0,1), (min(sce0.values()), max(sce0.values())), (min(selected_fc0.values()), max(selected_fc0.values())), (0,1), (min(sce3.values()), max(sce3.values())), (min(selected_fc3.values()), max(selected_fc3.values()))
ticks1, ticks2, ticks3, ticks4, ticks5, ticks6 = [ticks1[0], (ticks1[1]+ticks1[0])/2, ticks1[1]], [ticks2[0], (ticks2[1]+ticks2[0])/2, ticks2[1]], [ticks3[0], (ticks3[1]+ticks3[0])/2, ticks3[1]], [ticks4[0], (ticks4[1]+ticks4[0])/2, ticks4[1]], [ticks5[0], (ticks5[1]+ticks5[0])/2, ticks5[1]], [ticks6[0], (ticks6[1]+ticks6[0])/2, ticks6[1]]
ticks = ticks1, ticks2, ticks3, ticks4, ticks5, ticks6
_ = [fig.colorbar(plt.cm.ScalarMappable(cmap='Greens', norm=colors.Normalize(vmax=ticks[i][-1], vmin=ticks[i][0])), orientation='horizontal', ax=ax[i], fraction=.03, pad=.02, shrink=.9) for i in [0,3]]
_ = [fig.colorbar(plt.cm.ScalarMappable(cmap='terrain_r', norm=colors.Normalize(vmax=ticks[i][-1], vmin=ticks[i][0])), orientation='horizontal', ax=ax[i], fraction=.03, pad=.02, shrink=.9) for i in [1,2,4,5]]
#ax[1].set_title(str(round(SCFC(sce0, selected_fc0, not_sinks),2)))
#ax[4].set_title(str(round(SCFC(sce3, selected_fc3, not_sinks),2)))
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[0], fraction=.03, pad=.02, shrink=.9)
cbar.ax.set_title('vegetation density', size=10)
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[1], fraction=.03, pad=.02, shrink=.9)
cbar.ax.set_title(' ', size=10)
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[2], fraction=.03, pad=.02, shrink=.9)
cbar.ax.set_title('structural connectivity', size=10)
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[3], fraction=.03, pad=.02, shrink=.9)
cbar.ax.set_title(' ', size=9.5)
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens'), orientation='horizontal', ax=ax[4], fraction=.03, pad=.02, shrink=.9)
cbar.ax.set_title('microtopography (µm)', size=10)
plt.show()

0.61, 0.78

ax[0].plot(range(201), avgSCFC0, lw=1, c='lime')#or honeydew
ax[0].plot(range(201), avgSCFC1, lw=1, c='bisque')#or yellow
ax[0].set_ylim(bottom=.15)
ax[0].legend()
ax[0].set_xticks([0]+[40*k for k in range(1,6)], labels=[1]+[200*k for k in range(1,6)])

ax[0].set_title('(a)', loc='left')

r0.truncated_power_law.plot_pdf(color='lightgreen', linestyle='--', ax=ax[1], alpha=.9)
ax[1].scatter(bx0, by0/sum(by0), marker='.', color='forestgreen', alpha=.8, s=60, label='grassland')
r1.truncated_power_law.plot_pdf(color='sandybrown', linestyle='--', ax=ax[1], alpha=.9)
ax[1].scatter(bx1, by1/sum(by1), marker='.', color='darkorange', alpha=.8, s=60, label='shrubland')
ax[1].legend(title='observed distribution', alignment='left', loc='lower left')
ax[1].set_xlabel('avalanche size')
ax[1].set_ylabel('frequency')
ax[1].set_title('(b)', loc='left')
plt.show()




from math import log10
t1 = sum(Sizes[1],[])
binning = np.logspace(0, log10(max(t1)), 30)
bx_1, by_1 = log_binning_with_input(Counter(t1), binning)
bx_1, by_1 = drop_zeros(bx_1), drop_zeros(by_1)

#binning = log_binning(Counter(t2))[0]
#binning = binning[np.where(~np.isnan(binning))]
t0 = sum(Sizes[0],[])
bx_0, by_0 = log_binning_with_input(Counter(t0), binning)
bx_0, by_0 = bx_0[~np.isnan(bx_0)], by_0[~np.isnan(by_0)]

plt.plot(bx_0, by_0, marker='.', color='darkgreen')
plt.plot(bx_1, by_1, marker='.', color='darkorange')
plt.yscale('log')
plt.xscale('log')
plt.show()

FCl_eff,FCv_eff, SCpFCl_eff,SCeFCv_eff,SCeFCl_eff,SCpFCv_eff, Sizes_eff_prev = pickle.load(open('time_convergence_stochastic_dyn.txt','rb'))



def SC_edge(g,epmap):
    gconn = deepcopy(g)
    gconn.remove_edges_from(list(g.edges()))
    for j in range(60):
        next_layer_successors, same_layer_successors = {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]-1] for v1 in [(i,j) for i in range(20)]}, {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]] for v1 in [(i,j) for i in range(20)]}
        for i in range(20):
            v1 = (i,j)
            gconn.add_weighted_edges_from([(v1,v2,epmap[(v1,v2)]) for v2 in next_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,gconn[v2][v3]['weight']) for v3 in gconn.successors(v2)] for v2 in next_layer_successors[v1]],[]))
        for v1 in sort_layer(same_layer_successors):
            gconn.add_weighted_edges_from([(v1,v2,epmap[(v1,v2)]) for v2 in same_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,gconn[v2][v3]['weight']) for v3 in gconn.successors(v2)] for v2 in same_layer_successors[v1]],[]))
    return {node:gconn.in_degree(node, weight='weight') for node in g}

def eprob_map(g, e):
    ep = {}
    for i in g:
        epi = {j: max(0,e[i]-e[j]) for j in g.successors(i)}
        sum_probs = sum(epi.values())
        for j in epi:
            if epi[j]:
                ep[(i,j)] = epi[j]/sum_probs
        if epi==[0] and len(list(g.successors(i)))==1:
            ep[(i,list(g.successors(i))[0])] = 1
    return ep

def SCw_path(g, c, epmap):
    sc, i, gconn = {}, 0, deepcopy(g)
    gconn.remove_edges_from(list(g.edges))
    for v1 in c:
        print('NODE N°',i,':',v1,'HAS',len(c[v1]),' possible targets...')
        i+=1
        for v2 in c[v1]:
            nb_steps = nx.shortest_path_length(g,v1,v2)
            if nb_steps<=30:
                paths = list(nx.all_shortest_paths(g,v1,v2))
                prod = np.sum([np.prod([epmap[(p[i],p[i+1])] for i in range(len(p)-1)]) for p in paths])
                gconn.add_edge(v1,v2,weight=prod)
    for v1 in c:
        sc[v1] = gconn.in_degree(v1,weight='weight') 
    return sc

def drop_zeros(a_list):
    #   return [i for i in a_list if i>0]
    return a_list[~np.isnan(a_list)]

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

def log_binning_with_input(counter_dict, overall_binning):
    bin_means_y = np.histogram(list(counter_dict.keys()), bins=overall_binning, weights=list(counter_dict.values()))[0] / np.histogram(list(counter_dict.keys()), bins=overall_binning)[0]
    bin_means_x = np.histogram(list(counter_dict.keys()), bins=overall_binning, weights=list(counter_dict.keys()))[0] / np.histogram(list(counter_dict.keys()), bins=overall_binning)[0]
    return bin_means_x, bin_means_y

### RERUN STOCHASTIC ON 600 GRIDS
nb_steps = 100000
Vs, Es_grass, Es_shrub = pickle.load(open('data2025/Vs.txt','rb')), pickle.load(open('data2025/Es_grass.txt','rb')), pickle.load(open('data2025/Es_shrub.txt','rb'))
for p in range(3):
    print('P='+str(p))
    for t in range(100):
        print('.. trial n°'+str(t))
        V, E_grass, E_shrub = Vs[p,t], Es_grass[p,t], Es_shrub[p,t]
        Geff_grass, Geff_shrub = eff_lattice(G,E_grass), eff_lattice(G,E_shrub)
        Ep_grass, Ep_shrub = eprob(Geff_grass,E_grass), eprob(Geff_shrub,E_shrub)
        Sizes_grass, Durations_grass = np.zeros(100, list), np.zeros(100, list)
        Sizes_shrub, Durations_shrub = np.zeros(100, list), np.zeros(100, list)
        for run in range(100):
            print('... run n°'+str(run)+'/100')
            S_grass, current_grass, branches_grass, active_grass, Sizes_grass[run], Durations_grass[run] = simple_init(Geff_grass)
            S_shrub, current_shrub, branches_shrub, active_shrub, Sizes_shrub[run], Durations_shrub[run] = simple_init(Geff_shrub)
            for step in range(nb_steps):
                S_grass, current_grass, branches_grass, active_grass, Sizes_grass[run], Durations_grass[run] = simple_sandpile(Geff_grass, V, Ep_grass, S_grass, current_grass, branches_grass, active_grass, Sizes_grass[run], Durations_grass[run], not_sinks)
                S_shrub, current_shrub, branches_shrub, active_shrub, Sizes_shrub[run], Durations_shrub[run] = simple_sandpile(Geff_shrub, V, Ep_shrub, S_shrub, current_shrub, branches_shrub, active_shrub, Sizes_shrub[run], Durations_shrub[run], not_sinks)
        with open('data2025/Sizes_grass_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Sizes_grass, f)
        with open('data2025/Durations_grass_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Durations_grass, f)
        with open('data2025/Sizes_shrub_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Sizes_shrub, f)
        with open('data2025/Durations_shrub_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Durations_shrub, f)
        del Sizes_grass, Sizes_shrub, Durations_grass, Durations_shrub


def simple_sandpile(g, v, ep, state, current_av, branches, old_active, size, duration, not_sinks): # produces properly long avalanches
    new_active, unstable = [], [node for node in old_active if state[node] >= g.out_degree(node)]
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                spreading_scheme = Counter(random.choices(list(g.successors(node1)), weights=ep[(node1)], k=nb_partcl))
                for node2 in spreading_scheme:
                    if node2 in not_sinks:
                        state[node2] += np.sum(np.random.rand(spreading_scheme[node2]) > v[node2]/2)
                    current_av[node2] = current_av[node1]+1
                    if node2 not in branches:
                        branches[node2] = [node1]
                    elif node1 not in branches[node2]:
                        branches[node2].append(node1)
                    if node2 not in new_active:
                        new_active.append(node2)
    else:
        if len(current_av):
            size.append(len(current_av))
            duration.append(max(list(current_av.values())))
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, current_av, branches, new_active, size, duration


def simple_init(g):
    state, current_av, branches, new_active, size, dur = {}, [], {}, [], [], []
    for i in g:
        state[i] = 0
    return state, current_av, branches, new_active, size, dur

### ANALYZING PRELIMINARY AVALANCHE SLOPES

import pickle
import numpy as np
from math import log10
max_t, nb_runs = 100, 100

Sizes = np.zeros((2, max_t, nb_runs), list)
for t in range(100):
    t
    Sizes[0,t,:], Sizes[1,t,:] = pickle.load(open('data2025/Sizes_grass_p0_t'+str(t)+'.txt','rb')), pickle.load(open('data2025/Sizes_shrub_p0_t'+str(t)+'.txt','rb'))


t0, t1 = sum(Sizes[1,t],[]), sum(Sizes[0,t],[])
t0, t1 = np.array(big_sum0)-1, np.array(big_sum1)-1
binning = np.logspace(0, log10(max(t1)), 30)
bx_1, by_1 = log_binning_with_input(Counter(t1), binning)
bx_1, by_1 = bx_1[~np.isnan(bx_1)], by_1[~np.isnan(by_1)]

bx_0, by_0 = log_binning_with_input(Counter(t0), binning)
bx_0, by_0 = bx_0[~np.isnan(bx_0)], by_0[~np.isnan(by_0)]

plt.plot(bx_0, by_0, marker='.', color='darkgreen')
plt.plot(bx_1, by_1, marker='.', color='darkorange')
plt.yscale('log')
plt.xscale('log')
plt.show()
big_sum0 = sum([sum(list(Sizes[0,t,:]),[]) for t in range(50)],[])
big_sum1 = sum([sum(list(Sizes[1,t,:]),[]) for t in range(50)],[])

BX, BY = np.zeros((2,3,100), np.ndarray), np.zeros((2,3,100), np.ndarray)
for p in range(3):
    for t in range(100):
        Sizes_grass, Sizes_shrub = np.array(sum(pickle.load(open('data2025/Sizes_grass_p'+str(p)+'_t'+str(t)+'.txt','rb')),[]))-1, np.array(sum(pickle.load(open('data2025/Sizes_shrub_p'+str(p)+'_t'+str(t)+'.txt','rb')),[]))-1
        bx, by = log_binning_with_input(Counter(Sizes_grass), binning)
        BX[0,p,t], BY[0,p,t] = bx[~np.isnan(bx)], by[~np.isnan(by)]
        bx, by = log_binning_with_input(Counter(Sizes_shrub), binning)
        BX[1,p,t], BY[1,p,t] = bx[~np.isnan(bx)], by[~np.isnan(by)]

p = 0
for t in range(3):
    _ = plt.plot(BX[0,p,t], BY[0,p,t], marker='.', color='darkgreen')
    _ = plt.plot(BX[1,p,t], BY[1,p,t], marker='.', color='darkorange')

plt.loglog([])
plt.show()

by_adjusted0 = np.array([np.hstack((BY[0,0,t], np.zeros(26-len(BX[0,0,t])))) for t in range(100)])
by_adjusted1 = np.array([np.hstack((BY[1,0,t], np.zeros(26-len(BX[1,0,t])))) for t in range(100)])
by_mean0, by_std0 = np.array([np.mean([by_adjusted0[t][k] for t in range(100)]) for k in range(26)]), np.array([np.std([by_adjusted0[t][k] for t in range(100)]) for k in range(26)])
by_mean1, by_std1 = np.array([np.mean([by_adjusted1[t][k] for t in range(100)]) for k in range(26)]), np.array([np.std([by_adjusted1[t][k] for t in range(100)]) for k in range(26)])

plt.fill_between(BX[1,0,1], by_mean0-by_std0, by_mean0+by_std0, color='darkgreen', alpha=.7)
plt.fill_between(BX[1,0,1], by_mean1-by_std1, by_mean1+by_std1, color='darkorange', alpha=.7)
plt.loglog([])
plt.show()



np.argmax([sum(list(Sizes[1,t,:]),[]) for t in range(100)])
_ = [max(sum(Sizes[1,t],[])) for t in range(100)]

nb_steps = 100000
Vs, Es_grass, Es_shrub = pickle.load(open('data2025/Vs.txt','rb')), pickle.load(open('data2025/Es_grass.txt','rb')), pickle.load(open('data2025/Es_shrub.txt','rb'))
for p in range(3):
    print('P='+str(p))
    for t in range(100):
        print('.. trial n°'+str(t))
        V, E_grass, E_shrub = Vs[p,t], Es_grass[p,t], Es_shrub[p,t]
        G_grass, G_shrub = eff_lattice(G,E_grass), eff_lattice(G,E_shrub)
        Ep_grass, Ep_shrub = eprob(G_grass,E_grass), eprob(G_shrub,E_shrub)
        Epm_grass, Epm_shrub = eprob_map(G_grass, E_grass), eprob_map(G_shrub, E_shrub)
        _, Coupling_grass, _, _, _, _ = init(G_grass)
        _, Coupling_shrub, _, _, _, _ = init(G_shrub)
        #sce_grass, sce_shrub = SC_edge(G_grass, Epm_grass), SC_edge(G_shrub, Epm_shrub)
        with open('data2025/Coupling_grass_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Coupling_grass, f)
        with open('data2025/Coupling_shrub_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(Coupling_shrub, f)

Vs, Es_grass, Es_shrub = pickle.load(open('data2025/Vs.txt','rb')), pickle.load(open('data2025/Es_grass.txt','rb')), pickle.load(open('data2025/Es_shrub.txt','rb'))
SCFCs = np.zeros((2,3,100), float)
for p in range(3):
    print('P='+str(p))
    for t in range(100):
        print('.. trial n°'+str(t))
        Coupling_grass, Coupling_shrub = pickle.load(open('data2025/Coupling_grass_p'+str(p)+'_t'+str(t)+'.txt','rb')), pickle.load(open('data2025/Coupling_shrub_p'+str(p)+'_t'+str(t)+'.txt','rb'))
        V, E_grass, E_shrub = Vs[p,t], Es_grass[p,t], Es_shrub[p,t]
        G_grass, G_shrub = eff_lattice(G,E_grass), eff_lattice(G,E_shrub)
        Ep_grass, Ep_shrub = eprob(G_grass,E_grass), eprob(G_shrub,E_shrub)
        Epm_grass, Epm_shrub = eprob_map(G_grass, E_grass), eprob_map(G_shrub, E_shrub)
        sce_grass, sce_shrub = SC_edge(G_grass, Epm_grass), SC_edge(G_shrub, Epm_shrub)
        FCv_grass, FCv_shrub = {i:0 for i in G_grass}, {i:0 for i in G_shrub}
        for run in range(100):
            print('... run°'+str(run))
            S_grass, C_grass, current_grass, branches_grass, active_grass = zeros_no_avalanches(Coupling_grass)
            S_shrub, C_shrub, current_shrub, branches_shrub, active_shrub = zeros_no_avalanches(Coupling_shrub)
            for step in range(nb_steps):
                S_grass, C_grass, current_grass, branches_grass, active_grass = sandpile_only_scfc(G_grass, V, Ep_grass, S_grass, C_grass, current_grass, branches_grass, active_grass, not_sinks)
                S_shrub, C_shrub, current_shrub, branches_shrub, active_shrub = sandpile_only_scfc(G_shrub, V, Ep_shrub, S_shrub, C_shrub, current_shrub, branches_shrub, active_shrub, not_sinks)
            fc_grass, fc_shrub = FC_transport(C_grass), FC_transport(C_shrub)
            for i in G:
                FCv_grass[i] += fc_grass[i]
                FCv_shrub[i] += fc_shrub[i]
        SCFCs[0,p,t] = SCFC(sce_grass, FCv_grass, not_sinks)
        SCFCs[1,p,t] = SCFC(sce_shrub, FCv_shrub, not_sinks)
        with open('data2025/SCFC_grass_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(SCFCs[0,p,t], f)
        with open('data2025/SCFC_shrub_p'+str(p)+'_t'+str(t)+'.txt','wb') as f:
            pickle.dump(SCFCs[1,p,t], f)



### COMPUTING SLOPES (POWERLAW WITH CUTOFF)

import powerlaw

fit_exp_grass, fit_exp_shrub = np.zeros((2,100),float), np.zeros((2,100),float)
fit_pl_grass, fit_pl_shrub = np.zeros((2,100),float), np.zeros((2,100),float)
for k in range(100):
    k
    r0, r1 = powerlaw.Fit(Sizes_grass[k], xmin=xmin), powerlaw.Fit(Sizes_shrub[k], xmin=xmin)
    fit_exp_grass[:,k] = r0.stretched_exponential.parameter1, r0.stretched_exponential.parameter2
    fit_exp_shrub[:,k] = r1.stretched_exponential.parameter1, r1.stretched_exponential.parameter2
    fit_pl_grass[:,k] = r0.truncated_power_law.parameter1, r0.truncated_power_law.parameter2
    fit_pl_shrub[:,k] = r1.truncated_power_law.parameter1, r1.truncated_power_law.parameter2

# stretch exp
print('GRASS exp slope ='+str(round(np.mean(fit_exp_grass[0]),2))+' (+-'+str(round(np.std(fit_exp_grass[0]),2))+')')
print('SHRUB exp slope ='+str(round(np.mean(fit_exp_shrub[0]),2))+' (+-'+str(round(np.std(fit_exp_shrub[0]),2))+')')
print('GRASS scaling x0 ='+str(round(np.mean(fit_exp_grass[1]),2))+' (+-'+str(round(np.std(fit_exp_grass[1]),2))+')')
print('SHRUB scaling x0 ='+str(round(np.mean(fit_exp_shrub[1]),2))+' (+-'+str(round(np.std(fit_exp_shrub[1]),2))+')')

# trunc powerlaw
print('GRASS pl slope ='+str(round(np.mean(fit_pl_grass[0]),2))+' (+-'+str(round(np.std(fit_pl_grass[0]),2))+')')
print('SHRUB pl slope ='+str(round(np.mean(fit_pl_grass[0]),2))+' (+-'+str(round(np.std(fit_pl_grass[0]),2))+')')
print('GRASS cutoff ='+str(round(np.mean(fit_pl_shrub[1]),2))+' (+-'+str(round(np.std(fit_pl_shrub[1]),2))+')')
print('SHRUB cutoff ='+str(round(np.mean(fit_pl_shrub[1]),2))+' (+-'+str(round(np.std(fit_pl_shrub[1])),2)+')')





test0, test1 = sum(Sizes_grass,[]), sum(Sizes_shrub,[])
r0, r1 = powerlaw.Fit(Sizes_grass[k], xmin=xmin), powerlaw.Fit(Sizes_shrub[k], xmin=xmin)
print(r0.truncated_power_law.alpha, r1.truncated_power_law.alpha)
print(r0.truncated_power_law.D, r1.truncated_power_law.D)

r0, r1 = powerlaw.Fit(Sizes_grass[k], xmin=xmin), powerlaw.Fit(Sizes_shrub[k], xmin=xmin)
r0.distribution_compare('truncated_power_law', 'stretched_exponential')[0], r1.distribution_compare('truncated_power_law', 'stretched_exponential')[0]

p, t = 0, 0
Sizes_grass_syn, Sizes_shrub_syn = pickle.load(open('data2025/Sizes_grass_p'+str(p)+'_t'+str(t)+'.txt','rb')), pickle.load(open('data2025/Sizes_shrub_p'+str(p)+'_t'+str(t)+'.txt','rb'))
r0, r1 = powerlaw.Fit(Sizes_grass_syn[k], xmin=xmin), powerlaw.Fit(Sizes_shrub_syn[k], xmin=xmin)
r0.distribution_compare('truncated_power_law', 'stretched_exponential')[0], r1.distribution_compare('truncated_power_law', 'stretched_exponential')[0]

fit_grass, fit_shrub = np.zeros((2,3,100,100),object), np.zeros((2,3,100,100),object)
for p in range(3):
    print('P='+str(p))
    for t in range(100):
        print(' t='+str(t))
        Sizes_grass_syn, Sizes_shrub_syn = pickle.load(open('data2025/Sizes_grass_p'+str(p)+'_t'+str(t)+'.txt','rb')), pickle.load(open('data2025/Sizes_shrub_p'+str(p)+'_t'+str(t)+'.txt','rb'))
        for k in range(100):
            for xmin in [0,1]:
                fit_grass[xmin,p,t,k], fit_shrub[xmin,p,t,k] = powerlaw.Fit(Sizes_grass_syn[k], xmin=xmin+1), powerlaw.Fit(Sizes_shrub_syn[k], xmin=xmin+1)

SCFC_grass, SCFC_shrub = np.array([[pickle.load(open('data2025/SCFC_grass_p'+str(p)+'_t'+str(t)+'.txt','rb')) for t in range(100)] for p in range(3)]), np.array([[pickle.load(open('data2025/SCFC_shrub_p'+str(p)+'_t'+str(t)+'.txt','rb')) for t in range(100)] for p in range(3)])

exp_params, truncpl_params = np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float)

stretchexp_slope_grass, stretchexp_cutoff_grass, stretchexp_slope_shrub, stretchexp_cutoff_shrub = np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float)
truncpl_slope_grass, truncpl_scale_grass, truncpl_slope_shrub, truncpl_scale_shrub = np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float), np.zeros((2,3,100,100),float)
for p in range(1):
    print('P='+str(p))
    for t in range(100):
        print(' t='+str(t))
        for k in range(100):
            #k
            for xmin in [1]:
            #for xmin in [0]:
                #fit_grass[xmin,p,t,k], fit_shrub[xmin,p,t,k] = powerlaw.Fit(Sizes_grass_syn[k], xmin=xmin+1), powerlaw.Fit(Sizes_shrub_syn[k], xmin=xmin+1)
                #stretchexp_slope_grass[xmin,p,t,k], stretchexp_cutoff_grass[xmin,p,t,k] = fit_grass[xmin,p,t,k].stretched_exponential.parameter1, fit_grass[xmin,p,t,k].stretched_exponential.parameter2 # done for xmin=0
                #stretchexp_slope_shrub[xmin,p,t,k], stretchexp_cutoff_shrub[xmin,p,t,k] = fit_shrub[xmin,p,t,k].stretched_exponential.parameter1, fit_shrub[xmin,p,t,k].stretched_exponential.parameter2 # done for xmin=0
                truncpl_slope_grass[xmin,p,t,k], truncpl_scale_grass[xmin,p,t,k] = fit_grass[xmin,p,t,k].truncated_power_law.parameter1, fit_grass[xmin,p,t,k].truncated_power_law.parameter2
                truncpl_slope_shrub[xmin,p,t,k], truncpl_scale_shrub[xmin,p,t,k] = fit_shrub[xmin,p,t,k].truncated_power_law.parameter1, fit_shrub[xmin,p,t,k].truncated_power_law.parameter2

_ = stretchexp_slope_grass, stretchexp_cutoff_grass, stretchexp_slope_shrub, stretchexp_cutoff_shrub
with open('data2025/stretched_exponential_analysis.txt','wb') as f:
    pickle.dump(_,f)

_ = truncpl_slope_grass, truncpl_scale_grass, truncpl_slope_shrub, truncpl_scale_shrub
with open('data2025/truncated_powerlaw_analysis.txt','wb') as f:
    pickle.dump(_,f)

stretchexp_slope_grass, stretchexp_cutoff_grass, stretchexp_slope_shrub, stretchexp_cutoff_shrub = pickle.load(open('data2025/stretched_exponential_analysis.txt','rb'))
truncpl_slope_grass, truncpl_scale_grass, truncpl_slope_shrub, truncpl_scale_shrub = pickle.load(open('data2025/truncated_powerlaw_analysis.txt','rb'))


## STRETCHED EXPONENTIAL VISUALISATION

xlim_left, xlim_right = .3,.82
ylim_bottom, ylim_top = .085,.26
c_avg = 'k'
fig0, ax = plt.subplots(2,3, figsize=(8.5,5), constrained_layout=True)
for p in range(3):
    title_grass = '(a)'*(p==0) + '(b)'*(p==1) + '(c)'*(p==2)
    title_shrub = '(d)'*(p==0) + '(e)'*(p==1) + '(f)'*(p==2)
    _ = ax[0,p].scatter(SCFC_grass[p], [np.mean(stretchexp_slope_grass[xmin,p,t]) for t in range(100)], marker='.', c='green', alpha=.7)
    _ = ax[1,p].scatter(SCFC_shrub[p], [np.mean(stretchexp_slope_shrub[xmin,p,t]) for t in range(100)], marker='.', c='darkorange', alpha=.7)
    SCFC_grass_avg, SCFC_shrub_avg = np.mean(SCFC_grass[p]), np.mean(SCFC_shrub[p])
    slope_grass_avg, slope_shrub_avg = np.mean([np.mean(stretchexp_slope_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(stretchexp_slope_shrub[xmin,p,t]) for t in range(100)])
    corr_grass, corr_shrub = round(pearsonr(SCFC_grass[p], [np.mean(stretchexp_slope_grass[xmin,p,t]) for t in range(100)])[0],2), round(pearsonr(SCFC_shrub[p], [np.mean(stretchexp_slope_shrub[xmin,p,t]) for t in range(100)])[0],2)
    _ = ax[0,p].scatter(SCFC_grass_avg, slope_grass_avg, marker='x', c=c_avg, s=30)
    _ = ax[1,p].scatter(SCFC_shrub_avg, slope_shrub_avg, marker='x', c=c_avg, s=30)
    _ = ax[0,p].set_xlim(xlim_left, xlim_right)
    _ = ax[0,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[1,p].set_xlim(xlim_left, xlim_right)
    _ = ax[1,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[0,p].plot([xlim_left, SCFC_grass_avg], [slope_grass_avg, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([xlim_left, SCFC_shrub_avg], [slope_shrub_avg, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].plot([SCFC_grass_avg, SCFC_grass_avg], [ylim_bottom, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([SCFC_shrub_avg, SCFC_shrub_avg], [ylim_bottom, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].set_title(title_grass+' p$_c$='+str(.2+p*.3)+', r='+str(corr_grass), loc='left')
    _ = ax[1,p].set_title(title_shrub+' p$_c$='+str(.2+p*.3)+', r='+str(corr_shrub), loc='left')
    _ = ax[0,p].set_xlabel('SC-FC correlation')
    _ = ax[1,p].set_xlabel('SC-FC correlation')

ax[0,0].set_ylabel('grassland\nstretching exponent $\\beta$')
ax[1,0].set_ylabel('shrubland\nstretching exponent $\\beta$')
fig0.show()

xlim_left, xlim_right = .3,.815
ylim_bottom, ylim_top = .375,.98
fig1, ax = plt.subplots(2,3, figsize=(8.5,5), constrained_layout=True)
for p in range(3):
    title_grass = '(a)'*(p==0) + '(b)'*(p==1) + '(c)'*(p==2)
    title_shrub = '(d)'*(p==0) + '(e)'*(p==1) + '(f)'*(p==2)
    _ = ax[0,p].scatter(SCFC_grass[p], [np.mean(stretchexp_scale_grass[xmin,p,t]) for t in range(100)], marker='.', c='darkgreen', alpha=.7)
    _ = ax[1,p].scatter(SCFC_shrub[p], [np.mean(stretchexp_scale_shrub[xmin,p,t]) for t in range(100)], marker='.', c='darkorange', alpha=.7)
    SCFC_grass_avg, SCFC_shrub_avg = np.mean(SCFC_grass[p]), np.mean(SCFC_shrub[p])
    scale_grass_avg, scale_shrub_avg = np.mean([np.mean(stretchexp_scale_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(stretchexp_scale_shrub[xmin,p,t]) for t in range(100)])
    corr_grass, corr_shrub = round(pearsonr(SCFC_grass[p], [np.mean(stretchexp_scale_grass[xmin,p,t]) for t in range(100)])[0],2), round(pearsonr(SCFC_shrub[p], [np.mean(stretchexp_scale_shrub[xmin,p,t]) for t in range(100)])[0],2)
    _ = ax[0,p].scatter(SCFC_grass_avg, scale_grass_avg, marker='x', c=c_avg, s=30)
    _ = ax[1,p].scatter(SCFC_shrub_avg, scale_shrub_avg, marker='x', c=c_avg, s=30)
    _ = ax[0,p].set_xlim(xlim_left, xlim_right)
    _ = ax[0,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[1,p].set_xlim(xlim_left, xlim_right)
    _ = ax[1,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[0,p].plot([xlim_left, SCFC_grass_avg], [scale_grass_avg, scale_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([xlim_left, SCFC_shrub_avg], [scale_shrub_avg, scale_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].plot([SCFC_grass_avg, SCFC_grass_avg], [ylim_bottom, scale_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([SCFC_shrub_avg, SCFC_shrub_avg], [ylim_bottom, scale_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].set_title(title_grass+' p$_c$='+str(.2+p*.3)+', r='+str(corr_grass), loc='left')
    _ = ax[1,p].set_title(title_shrub+' p$_c$='+str(.2+p*.3)+', r='+str(corr_shrub), loc='left')
    _ = ax[0,p].set_xlabel('SC-FC correlation')
    _ = ax[1,p].set_xlabel('SC-FC correlation')

ax[0,0].set_ylabel('grassland\ncharacteristic scale s$_0$')
ax[1,0].set_ylabel('shrubland\ncharacteristic scale s$_0$')
fig1.show()

### TRUNCATED POWERLAW FIT VISUALISATION


xlim_left, xlim_right = .22, .82
ylim_bottom, ylim_top = 1.11, 1.83
c_avg = 'k'
fig0, ax = plt.subplots(2,3, figsize=(8.5,5), constrained_layout=True)
for p in range(3):
    title_grass = '(a)'*(p==0) + '(b)'*(p==1) + '(c)'*(p==2)
    title_shrub = '(d)'*(p==0) + '(e)'*(p==1) + '(f)'*(p==2)
    _ = ax[0,p].scatter(SCFC_grass[p], [np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)], marker='.', c='green', alpha=.7)
    _ = ax[1,p].scatter(SCFC_shrub[p], [np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)], marker='.', c='darkorange', alpha=.7)
    SCFC_grass_avg, SCFC_shrub_avg = np.mean(SCFC_grass[p]), np.mean(SCFC_shrub[p])
    slope_grass_avg, slope_shrub_avg = np.mean([np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])
    corr_grass, corr_shrub = round(pearsonr(SCFC_grass[p], [np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)])[0],2), round(pearsonr(SCFC_shrub[p], [np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])[0],2)
    _ = ax[0,p].scatter(SCFC_grass_avg, slope_grass_avg, marker='x', c=c_avg, s=30)
    _ = ax[1,p].scatter(SCFC_shrub_avg, slope_shrub_avg, marker='x', c=c_avg, s=30)
    _ = ax[0,p].set_xlim(xlim_left, xlim_right)
    _ = ax[0,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[1,p].set_xlim(xlim_left, xlim_right)
    _ = ax[1,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[0,p].plot([xlim_left, SCFC_grass_avg], [slope_grass_avg, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([xlim_left, SCFC_shrub_avg], [slope_shrub_avg, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].plot([SCFC_grass_avg, SCFC_grass_avg], [ylim_bottom, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([SCFC_shrub_avg, SCFC_shrub_avg], [ylim_bottom, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].set_title(title_grass+' p$_c$='+str(.2+p*.3)+', r='+str(corr_grass), loc='left')
    _ = ax[1,p].set_title(title_shrub+' p$_c$='+str(.2+p*.3)+', r='+str(corr_shrub), loc='left')
    _ = ax[0,p].set_xlabel('SC-FC correlation')
    _ = ax[1,p].set_xlabel('SC-FC correlation')

ax[0,0].set_ylabel('grassland\npower law slope $\\alpha$')
ax[1,0].set_ylabel('shrubland\npower law slope $\\alpha$')
fig0.show()

xlim_left, xlim_right = .22, .82
ylim_bottom, ylim_top = .001, .055
#ylim_bottom, ylim_top = 1.13, 1.82
clim_grass_low, clim_grass_high = 1.14, 1.62
clim_shrub_low, clim_shrub_high = 1.36, 1.81
clim_both_low, clim_both_high = 1.14, 1.81
#clim_both_low, clim_both_high = 0.005, 0.045
fig1, ax = plt.subplots(2,3, figsize=(8.5,5), constrained_layout=True)
for p in range(3):
    title_grass = '(a)'*(p==0) + '(b)'*(p==1) + '(c)'*(p==2)
    title_shrub = '(d)'*(p==0) + '(e)'*(p==1) + '(f)'*(p==2)
    #_ = ax[0,p].scatter(SCFC_grass[p], [np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)], marker='.', c='darkgreen', alpha=.7)
    #_ = ax[1,p].scatter(SCFC_shrub[p], [np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)], marker='.', c='darkorange', alpha=.7)
    c_grass = ax[0,p].scatter(SCFC_grass[p], [np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)], marker='.', c=[np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)], cmap='terrain_r' , alpha=.7, s=40)
    c_shrub = ax[1,p].scatter(SCFC_shrub[p], [np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)], marker='.', c=[np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)], cmap='terrain_r', alpha=.7, s=40)
    #c_grass.set_clim(clim_grass_low, clim_grass_high)
    #c_shrub.set_clim(clim_shrub_low, clim_shrub_high)
    c_grass.set_clim(clim_both_low, clim_both_high)
    c_shrub.set_clim(clim_both_low, clim_both_high)
    if p==2:
        cbar_grass = plt.colorbar(c_grass, ax=ax[0,p])
        cbar_shrub = plt.colorbar(c_shrub, ax=ax[1,p])
        cbar_grass.ax.set_ylabel('scaling exponent $\\alpha$', rotation=270, labelpad=15)
        cbar_shrub.ax.set_ylabel('scaling exponent $\\alpha$', rotation=270, labelpad=15)
        
        #cbar_grass = plt.colorbar(c_grass, ax=ax[0,p])
        #cbar_shrub = plt.colorbar(c_shrub, ax=ax[1,p])
        #cbar_grass.ax.set_ylabel('cut-off exponent s$_0$', rotation=270, labelpad=15)
        #cbar_shrub.ax.set_ylabel('cut-off exponent s$_0$', rotation=270, labelpad=15)
    SCFC_grass_avg, SCFC_shrub_avg = np.mean(SCFC_grass[p]), np.mean(SCFC_shrub[p])
    slope_grass_avg, slope_shrub_avg = np.mean([np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])
    scale_grass_avg, scale_shrub_avg = np.mean([np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)])
    corr_grass, corr_shrub = round(pearsonr(SCFC_grass[p], [np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)])[0],2), round(pearsonr(SCFC_shrub[p], [np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])[0],2)
    _ = ax[0,p].scatter(SCFC_grass_avg, slope_grass_avg, marker='x', c=c_avg, s=30)
    _ = ax[1,p].scatter(SCFC_shrub_avg, slope_shrub_avg, marker='x', c=c_avg, s=30)
    _ = ax[0,p].set_xlim(xlim_left, xlim_right)
    _ = ax[0,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[1,p].set_xlim(xlim_left, xlim_right)
    _ = ax[1,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[0,p].plot([xlim_left, SCFC_grass_avg], [scale_grass_avg, scale_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([xlim_left, SCFC_shrub_avg], [scale_shrub_avg, scale_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].plot([SCFC_grass_avg, SCFC_grass_avg], [ylim_bottom, scale_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([SCFC_shrub_avg, SCFC_shrub_avg], [ylim_bottom, scale_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].set_title(title_grass+' p$_c$='+str(.2+p*.3)+', r='+str(corr_grass), loc='left')
    _ = ax[1,p].set_title(title_shrub+' p$_c$='+str(.2+p*.3)+', r='+str(corr_shrub), loc='left')
    _ = ax[0,p].set_xlabel('SC-FC correlation')
    _ = ax[1,p].set_xlabel('SC-FC correlation')

ax[0,0].text(x=.25, y=0.0055, s='<$\\alpha$>=1.39', size=11)
ax[0,1].text(x=.26, y=0.0055, s='<$\\alpha$>=1.37', size=11)
ax[0,2].text(x=.26, y=0.0055, s='<$\\alpha$>=1.38', size=11)

ax[1,0].text(x=.52, y=0.048, s='<$\\alpha$>=1.56', size=11)
ax[1,1].text(x=.52, y=0.048, s='<$\\alpha$>=1.55', size=11)
ax[1,2].text(x=.52, y=0.048, s='<$\\alpha$>=1.57', size=11)
#1.3978
#1.5657
ax[0,0].set_ylabel('grassland\ncut-off parameter s$_0$')
ax[1,0].set_ylabel('shrubland\ncut-off parameter s$_0$')
#ax[0,0].set_ylabel('grassland\npower law slope $\\alpha$')
#ax[1,0].set_ylabel('shrubland\npower law slope $\\alpha$')
fig1.show()

with open('data2025/SCFC_grass.txt','wb') as f:
    pickle.dump(SCFC_grass,f)

with open('data2025/SCFC_shrub.txt','wb') as f:
    pickle.dump(SCFC_shrub,f)




### SLOPE (y-axis) VERSUS CUTOFF/SCALE (x-axis)
xlim_left, xlim_right = .001, .055
ylim_bottom, ylim_top = 1.11, 1.83
fig0, ax = plt.subplots(2,3, figsize=(8.5,5), constrained_layout=True)
for p in range(3):
    title_grass = '(a)'*(p==0) + '(b)'*(p==1) + '(c)'*(p==2)
    title_shrub = '(d)'*(p==0) + '(e)'*(p==1) + '(f)'*(p==2)
    _ = ax[0,p].scatter([np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)], [np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)], marker='.', c='green', alpha=.7)
    _ = ax[1,p].scatter([np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)], [np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)], marker='.', c='darkorange', alpha=.7)
    scale_grass_avg, scale_shrub_avg = np.mean([np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)])
    slope_grass_avg, slope_shrub_avg = np.mean([np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)]), np.mean([np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])
    corr_grass, corr_shrub = round(pearsonr([np.mean(truncpl_scale_grass[xmin,p,t]) for t in range(100)], [np.mean(truncpl_slope_grass[xmin,p,t]) for t in range(100)])[0],2), round(pearsonr([np.mean(truncpl_scale_shrub[xmin,p,t]) for t in range(100)], [np.mean(truncpl_slope_shrub[xmin,p,t]) for t in range(100)])[0],2)
    _ = ax[0,p].scatter(scale_grass_avg, slope_grass_avg, marker='x', c=c_avg, s=30)
    _ = ax[1,p].scatter(scale_shrub_avg, slope_shrub_avg, marker='x', c=c_avg, s=30)
    _ = ax[0,p].set_xlim(xlim_left, xlim_right)
    _ = ax[0,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[1,p].set_xlim(xlim_left, xlim_right)
    _ = ax[1,p].set_ylim(ylim_bottom, ylim_top)
    _ = ax[0,p].plot([xlim_left, scale_grass_avg], [slope_grass_avg, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([xlim_left, scale_shrub_avg], [slope_shrub_avg, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].plot([scale_grass_avg, scale_grass_avg], [ylim_bottom, slope_grass_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[1,p].plot([scale_shrub_avg, scale_shrub_avg], [ylim_bottom, slope_shrub_avg], lw=1.5, c=c_avg, linestyle='--', dashes=(4,5))
    _ = ax[0,p].set_title(title_grass+' p$_c$='+str(.2+p*.3)+', r='+str(corr_grass), loc='left')
    _ = ax[1,p].set_title(title_shrub+' p$_c$='+str(.2+p*.3)+', r='+str(corr_shrub), loc='left')
    _ = ax[0,p].set_xlabel('cut-off parameter s$_0$')
    _ = ax[1,p].set_xlabel('cut-off parameter s$_0$')

ax[0,0].set_ylabel('grassland\npower law slope $\\alpha$')
ax[1,0].set_ylabel('shrubland\npower law slope $\\alpha$')
fig0.show()



###

r0, r1 = powerlaw.Fit(test0), powerlaw.Fit(test1)
r0, r1 = powerlaw.Fit(test0, xmin=1), powerlaw.Fit(test1, xmin=1)

print(r0.truncated_power_law.alpha, r1.truncated_power_law.alpha)
print(r0.truncated_power_law.xmin, r1.truncated_power_law.xmin)

fig = plt.figure(constrained_layout=True)
r0.plot_pdf(color='darkgreen', linewidth=2, label='grass pdf')
r0.power_law.plot_pdf(color='darkgreen', linestyle='--', ax=fig)
r0.plot_ccdf(color='lightgreen', linewidth=2, ax=fig, label='grass cumul')
r0.power_law.plot_ccdf(color='lightgreen', linestyle='--', ax=fig)
r1.plot_pdf(color='darkorange', linewidth=2, label='shrub pdf')
r1.power_law.plot_pdf(color='darkorange', linestyle='--', ax=fig)
r1.plot_ccdf(color='orange', linewidth=2, ax=fig, label='shrub cumul')
r1.power_law.plot_ccdf(color='orange', linestyle='--', ax=fig)
plt.legend()
plt.show()

#r0, r1 = powerlaw.Fit(sizes0, xmin=xmin), powerlaw.Fit(sizes1, xmin=xmin)
r0, r1 = r0_xmin2, r0_xmin2

fig = r0.plot_pdf(color='darkgreen', linewidth=2, label='grass pdf')
r0.truncated_power_law.plot_pdf(color='darkgreen', linestyle='--', ax=fig)
r0.plot_ccdf(color='lightgreen', linewidth=2, ax=fig)
r0.truncated_power_law.plot_ccdf(color='lightgreen', linestyle='--', ax=fig)
r1.plot_pdf(color='darkorange', linewidth=2, label='shrub pdf')
r1.truncated_power_law.plot_pdf(color='darkorange', linestyle='--', ax=fig)
r1.plot_ccdf(color='orange', linewidth=2, ax=fig)
r1.truncated_power_law.plot_ccdf(color='orange', linestyle='--', ax=fig)
plt.legend()
plt.show()

sizes0, sizes1 = sum(Sizes_grass,[]), sum(Sizes_shrub,[])
bx0, by0 = log_binning_with_input(Counter(sizes0), binning)
bx1, by1 = log_binning_with_input(Counter(sizes1), binning)
bx0, by0 = bx0[~np.isnan(bx0)], by0[~np.isnan(by0)]
bx1, by1 = bx1[~np.isnan(bx1)], by1[~np.isnan(by1)]

r0, r1 = r0_xmin1, r1_xmin1
r0, r1 = r0_xmin2, r1_xmin2

#ax0 = r0.plot_pdf(color='darkgreen', linewidth=2, label='grass pdf')
fig, ax = plt.subplots(1, 2, figsize=(8.8,3.3), constrained_layout=True)
ax[0].set_xlabel('time step (10$^3$)')
ax[0].set_ylabel('correlation')
ax[0].set_title('(a)')
r0.truncated_power_law.plot_pdf(color='lightgreen', linestyle='--', ax=ax[1], alpha=.9)
ax[1].scatter(bx0, by0/sum(by0), marker='.', color='forestgreen', alpha=.8, s=60, label='grassland')
#ax1 = r1.plot_pdf(color='darkorange', linewidth=2)
r1.truncated_power_law.plot_pdf(color='sandybrown', linestyle='--', ax=ax[1], alpha=.9)
ax[1].scatter(bx1, by1/sum(by1), marker='.', color='darkorange', alpha=.8, s=60, label='shrubland')
#ax[1].legend(title='truncated power law fit', alignment='left', loc='lower left')
ax[1].legend(title='observed distribution', alignment='left', loc='lower left')
ax[1].set_xlabel('avalanche size')
ax[1].set_ylabel('frequency')
ax[1].set_title('(b)', loc='left')
ax[1].text(x=9, y=.0002, s='α = '+str(round(r0_xmin2.truncated_power_law.parameter1,2)), color='forestgreen', size=10)
ax[1].text(x=9, y=.000066, s='s  = '+str(round(r0_xmin2.truncated_power_law.parameter2*1e3))+' 10', color='forestgreen', size=10)
ax[1].text(x=35, y=.1, s='α = '+str(round(r1_xmin2.truncated_power_law.parameter1,2)), color='darkorange', size=10)
ax[1].text(x=35, y=.033, s='s  = '+str(round(r1_xmin2.truncated_power_law.parameter2*1e3))+' 10', color='darkorange', size=10)
plt.show()

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

fig, ax = plt.subplots(1, 2, figsize=(8.8,3.3), constrained_layout=True)
ax[0].fill_between(range(101), avgSCEFC0-stdSCEFC0, avgSCEFC0+stdSCEFC0, color='forestgreen', alpha=.8, label='grassland, SC$_E$ - FC(t)')
ax[0].fill_between(range(101), avgSCPDP0-stdSCPDP0, avgSCPDP0+stdSCPDP0, color='lightgreen', alpha=.8, label='grassland, SC$_P$ - DP(t)')
ax[0].fill_between(range(101), avgSCEFC1-stdSCEFC1, avgSCEFC1+stdSCEFC1, color='darkorange', alpha=.8, label='shrubland, SC$_E$ - FC(t)')
ax[0].fill_between(range(101), avgSCPDP1-stdSCPDP1, avgSCPDP1+stdSCPDP1, color='orange', alpha=.8, label='shrubland, SC$_P$ - DP(t)')
ax[0].plot(range(101), avgSCEFC0, lw=1, c='lime')
ax[0].plot(range(101), avgSCPDP0, lw=1, c='honeydew')
ax[0].plot(range(101), avgSCEFC1, lw=1, c='bisque')
ax[0].plot(range(101), avgSCPDP1, lw=1, c='yellow')
ax[0].set_ylim(bottom=.15)
ax[0].legend()
ax[0].set_xticks([0]+[20*k for k in range(1,6)], labels=[1]+[100*k for k in range(1,6)])
ax[0].set_xlabel('time step (10$^3$)')
ax[0].set_ylabel('correlation')
ax[0].set_title('(a)', loc='left')

ax[1].plot(bx_0, by_0, marker='.', color='darkgreen', label='grassland')
ax[1].plot(bx_1, by_1, marker='.', color='darkorange', label='shrubland')
ax[1].set_xlabel('avalanche size')
ax[1].set_ylabel('count')
ax[1].loglog([])
ax[1].legend(loc='lower left')
ax[1].set_title('(b)', loc='left')
plt.show()


for d in list(r0.supported_distributions.keys()):
    d, r0.distribution_compare('stretched_exponential', d)[0], r1.distribution_compare('stretched_exponential', d)[0]




{'power_law': <class 'powerlaw.Power_Law'>, 'lognormal': <class 'powerlaw.Lognormal'>, 'exponential': <class 'powerlaw.Exponential'>, 'truncated_power_law': <class 'powerlaw.Truncated_Power_Law'>, 'stretched_exponential': <class 'powerlaw.Stretched_Exponential'>, 'lognormal_positive': <class 'powerlaw.Lognormal_Positive'>}

p, t = 0, 0
data = pickle.load(np.array(sum(open('data2025/Coupling_grass_p'+str(p)+'_t'+str(t)+'.txt','rb'),[]))-1)
results = powerlaw.Fit(data)
results.truncated_power_law.parameter1
results.truncated_power_law.parameter2


Sizes_grass, Sizes_shrub = pickle.load(open('data2025/Sizes_grassland_natural.txt','rb')), pickle.load(open('data2025/Sizes_shrubland_natural.txt','rb'))

data0 = sum(Sizes_grass,[])
fit0 = powerlaw.Fit(data0)
a0 = fit0.truncated_power_law.parameter1
b0 = fit0.truncated_power_law.parameter2

data1 = sum(Sizes_shrub,[])
fit1 = powerlaw.Fit(data1)
a1 = fit1.truncated_power_law.parameter1
b1 = fit1.truncated_power_law.parameter2

bx_1, by_1 = log_binning_with_input(Counter(data1), binning)
bx_1, by_1 = bx_1[~np.isnan(bx_1)], by_1[~np.isnan(by_1)]
bx_0, by_0 = log_binning_with_input(Counter(data0), binning)
bx_0, by_0 = bx_0[~np.isnan(bx_0)], by_0[~np.isnan(by_0)]

plt.plot(bx_0, trunc(bx_0, a0, b0))
plt.loglog([])
plt.show()

plt.plot(bx_0, trunc(bx_0))

trunc_fit0 = fit0.truncated_power_law

fig = fit0.plot_pdf(color='b', linewidth=2, label='Empirical Data')
trunc_fit0.plot_pdf(ax=fig, color='r', linestyle='--', label='Truncated Power Law Fit')

wb2 = Fit_Weibull_2P(failures=data0)
wb3 = Fit_Weibull_3P(failures=data0)
plt.figure()
u1, u2 = np.unique(data0, return_counts=True)
plt.scatter(u1, u2, marker='o', s=20)
Y = plot_weibull(u1, wb.alpha, wb.beta)
plt.plot(u1, Y*sum(u2))
plt.loglog()
plt.title(str(grid))
plt.show()

from scipy.optimize import curve_fit
def truncated_power_law(x, alpha, lambd):
    return x**(-alpha) * np.exp(-lambd * x)

counts, bin_edges = np.histogram(sum(Sizes[0],[]), bins='auto', density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
initial_params = [1.5, 0.01]
params, covariance = curve_fit(truncated_power_law, bin_centers, counts, p0=initial_params)
alpha_fit, lambda_fit = params

from math import log10
from collections import Counter
from scipy.optimize import curve_fit

def truncated_power_law(x, alpha, lambd, beta):
    #return C * (x ** -alpha) * np.exp(-lambd * x)
    return (x ** -alpha) * np.exp(-(lambd * x)**beta)


binning = np.logspace(0, log10(max(t1)), 30)
bx_1, by_1 = log_binning_with_input(Counter(t1), binning)
bx_1, by_1 = drop_zeros(bx_1), drop_zeros(by_1)
bx_0, by_0 = log_binning_with_input(Counter(t0), binning)
bx_0, by_0 = bx_0[~np.isnan(bx_0)], by_0[~np.isnan(by_0)]

popt0,_=curve_fit(truncated_power_law, bx_0[3:], by_0[3:], p0=[2,.1,bx_0[3]])
popt1,_=curve_fit(truncated_power_law, bx_1, by_1, p0=[2,.1,bx_1[0]])

def trunc_func(x, alpha, x_min, x_max):
    return (x**(-alpha))*(x>=x_min)*(x<=x_max)

def trunc_func(x, alpha, x_min, x_max):
    return (x**(-alpha))*(x>=x_min)*(x<=x_max)

trunc_y = trunc_func(bx_1, 

plt.scatter(bx_0, by_0, label='grassland', marker='.')
plt.scatter(bx_1, by_1, label='shrubland', marker='.')
plt.loglog([])
plt.show()

import pymc as pm


x_min = 1.
with pm.Model() as model:
    alpha = pm.Uniform("alpha", lower=1.0, upper=3.0)  # Power-law exponent
    lambda_ = pm.Uniform("lam", lower=0.01, upper=1.0)
    def power_law_exp_cutoff_logp(x, alpha, lambda_, x_min):
        log_p = -alpha * pm.math.log(x) - lambda_ * x
        log_p = pm.math.switch(x >= x_min, log_p, -np.inf)  # Truncate below x_min
        return log_p
    likelihood = pm.DensityDist("likelihood", logp=power_law_exp_cutoff_logp, observed=data1, alpha=alpha, lambda_=lambda_, x_min=x_min)
    trace = pm.sample(2000, tune=1000, chains=2, target_accept=0.9)



>>> popt0
array([3.19553298e-01, 2.35079067e-01, 2.85656119e+06])
>>> by_0[:3]
array([2232083., 1542883.,  920747.])
>>> popt0,_=curve_fit(truncated_power_law, bx_0, by_0, p0=[2,.1, by_0[0]])
>>> popt
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'popt' is not defined. Did you mean: 'popt0'?
>>> popt0
array([3.19553636e-01, 2.35078903e-01, 2.85656075e+06])
>>> plt.scatter(bx_0, by_0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'plt' is not defined
>>> import matplotlib.pyplot as plt
>>> plt.scatter(bx_0, by_0)
<matplotlib.collections.PathCollection object at 0x000001817F8FF260>
>>> plt.loglog([])
[<matplotlib.lines.Line2D object at 0x0000018119213B00>]
>>> plt.plot(bx_0, truncated_power_law(bx_0,*popt0))
[<matplotlib.lines.Line2D object at 0x0000018119213E00>]
>>> plt.show()