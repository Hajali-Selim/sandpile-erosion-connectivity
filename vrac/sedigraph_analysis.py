from inits import *
import pandas as pd
# SEDIGRAPH

x1, x2 = 20,30
G, not_sinks = tgrid(x1,x2)

# merging 28.07.06, in names1[[14,15,21]]
extra = [14,15,16]

# STEP1: EXTRACT EVENTS

A = pd.ExcelFile('experimental_data_hydrograph/only_sed.xls')
nb_events = 22
data_tables2 = np.zeros(nb_events, np.ndarray)
for i in range(nb_events):
    data_tables2[i] = A.parse(sheet_name=i).values

names0, names00 = np.array(A.sheet_names), np.array([data_tables2[i][0,0] for i in range(nb_events)])
yields_monitored = {names00[i]:sum(all_sed[i]) for i in range(nb_events)}

new_ind = [14,15,20]
all_rain, all_sed = np.zeros(22, list), np.zeros(22, list)
all_rain[:21], all_sed[:21] = pickle.load(open('experimental_data_hydrograph/timeseries.txt','rb'))
all_rain[21], all_sed[21] = np.hstack(all_rain[new_ind]), np.hstack(all_sed[new_ind])

# let's readjust all_sed to all_rain
# by average
eventlen = [range(len(all_rain[i])) for i in select1+select4]
exp_sed, exp_rain = np.zeros(8, np.ndarray), np.zeros(8, np.ndarray)
for i in range(8):
    s_, r_ = all_sed[(select1+select4)[i]], all_rain[(select1+select4)[i]]
    exp_sed[i], exp_rain[i] = s_*max(r_)/max(s_), r_
    r_avg, s_avg = np.mean(r_[np.where(r_)]), np.mean(s_[np.where(s_)])
    exp_sed[i] = s_*r_avg/s_avg

fig, ax = plt.subplots(2,5,figsize=(19,7),tight_layout=True)
_ = [ax[i%2,i//2].plot(event_len[i], all_rain_p4[i]) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_rain_p4[0], label='rainfall amount')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4[i]) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4[0], label='sediment flow')
_ = [ax[i%2,i//2].set_title(names0[P4_loc0][i]) for i in range(10)]
_ = [ax[i%2,i//2].set_ylim(top=135) for i in range(10)]
ax[0,0].legend()
ax[1,2].set_xlabel('time')
plt.show()

# Merging 3 time-series 14+15+20, into 21


# total experimental yields
yields = {'P1 07.09.05':4113.7, 'P1 29.08.06':2048.12, 'P1 07.09.06':900, 'P4 29.09.05':443.33, 'P4 05.07.06':5835.96, 'P4 28.07.06m':4801.25, 'P4 23.08.06':2592.12, 'P4 29.08.06':0}
yields = {'P4 29.09.05':443.33, 'P4 05.07.06':5835.96, 'P4 28.07.06m':4801.25, 'P4 23.08.06':2592.12, 'P4 29.08.06':0}

# merging 280706a, 280706b, 280706m
extra_sed, extra_rain = np.zeros(3, np.ndarray), np.zeros(3, np.ndarray)
for i in range(3):
    s_, r_ = all_sed[extra[i]], all_rain[extra[i]]
    extra_sed[i], extra_rain[i] = s_*max(r_)/max(s_), r_
    r_avg, s_avg = np.mean(r_[np.where(r_)]), np.mean(s_[np.where(s_)])
    extra_sed[i] = s_*r_avg/s_avg

[sum(extra_sed[i]) for i in range(3)]
yields['P4 28.07.06'] = sum(sum(extra_sed[i]) for i in range(3))



# VISUALIZE MAHLERAN ON TOP OF EXPERIMENTAL SEDIGRAPHS
#select1, select4 = [0, 1, 2], [12, 13, 14, 19, 20]
new_ind = [14,15,20]
select1, select4 = [0, 1, 2], [12, 13, 19, 20, 21]

names1 = np.array(['P1 07.09.05', 'P1 29.08.06', 'P1 07.09.06', 'P2 07.09.05', 'P2 15.08.06', 'P2 07.09.06', 'P3 31.07.06', 'P3 01.08.06', 'P3 11.08.06', 'P3 29.08.06', 'P3 07.09.06', 'P4 07.09.05', 'P4 29.09.05', 'P4 05.07.06', 'P4 28.07.06a', 'P4 28.07.06b', 'P4 31.07.06', 'P4 01.08.06', 'P4 11.08.06', 'P4 23.08.06', 'P4 29.08.06', 'P4 28.07.06m', 'P4 29.07.06'])
folder = ['P1/plot 1 '+i[3:5]+i[6:8]+i[9:] for i in names1[select1]]+['P4/plot 4 '+i[3:5]+i[6:8]+i[9:] for i in names1[select4]]
# folder = ['plot 1 '+i[3:5]+i[6:8]+i[9:] for i in names1[select1]]+['plot 4 '+i[3:5]+i[6:8]+i[9:] for i in names1[select4]]
# folder[1] += '/uniform phi'
# folder[4] += ' detachment threshold and veg'
# folder[5] += ' calibrated'
# folder[6] += '/HZ = 0.004'

mahl_sed, mahl_hyd = np.zeros(len(select1+select4), np.ndarray), np.zeros(len(select1+select4), np.ndarray)
for i in range(len(select1+select4)):
    k = (names1+select4)[i]
    event = pd.read_excel('experimental_data_hydrograph/only_sed.xls', sheet_name='EVENT_'+str(k))
    sedtr, hydro = np.loadtxt('mahleran_outputs/'+folder[i]+'/sedtr001.dat', float), np.loadtxt('mahleran_outputs/'+folder[i]+'/hydro001.dat', float)
    sed = sedtr[:,1][:np.where(sedtr[:,1])[0][-1]]
    hyd = hydro[:len(sed)+1,1]
    names1[k], sum(sed), len(sed)
    mahl_sed[i], mahl_hyd[i] = np.array([np.mean(sed[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)]), np.array([np.mean(hyd[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)])

new_ind
sed_merge, hyd_merge = np.zeros(3, list), np.zeros(3, list)
for i in range(3):
    k = new_ind[i]
    name_i = names1[k]+'/hz = 0.004'
    folder_i = 'P4/plot 4 '+name_i[3:5]+name_i[6:8]+name_i[9:]
    #sed_merge[i], hyd_merge[i] = np.loadtxt('experimental_data_hydrograph/mahleran_outputs/'+folder_i+'/sedtr001.dat', float), np.loadtxt('experimental_data_hydrograph/mahleran_outputs/'+folder_i+'/hydro001.dat', float)[:,:12]
    sed_merge[i], hyd_merge[i] = np.loadtxt('mahleran_outputs/'+folder_i+'/sedtr001.dat', float, usecols=tuple(range(16))), np.loadtxt('mahleran_outputs/'+folder_i+'/hydro001.dat', float, usecols=tuple(range(12)))#[:,:12]
    i, sed_merge[i].shape, hyd_merge[i].shape

sed_merge, hyd_merge = np.vstack(list(sed_merge[i] for i in range(3))), np.vstack(list(hyd_merge[i] for i in range(3)))
np.savetxt('experimental_data_hydrograph/mahleran_outputs/P4/plot 4 280706m/sedtr001.dat',sed_merge)
np.savetxt('experimental_data_hydrograph/mahleran_outputs/P4/plot 4 280706m/hydro001.dat',hyd_merge)
i,k=-1,-1
s_, r_, mahl_ = all_sed[k], all_rain[k], mahl_sed[i]
sed_merge*max(r_)/max(mahl_)

all_sed_mahl, all_sed_aval = np.zeros(8, np.ndarray), np.zeros(8, np.ndarray)
for i in range(len(select1+select4)):
    k = (select1+select4)[i]
    s_, r_, mahl_ = all_sed[k], all_rain[k], mahl_sed[i]
    all_sed[i], all_rain[i], all_sed_mahl[i] = s_*max(r_)/max(s_), r_, mahl_*max(r_)/max(mahl_)

yields_mahl = np.array([sum(all_sed_mahl[i]) for i in range(8)])


# run sandpile sedigraph
# (1) steepest descent

Geff1, Geff4 = eff_lattice(G,EP[0]), eff_lattice(G,EP[3])
Gst1, Gst4 = eff_steep(Geff1, estp(Geff1,EP[0])), eff_steep(Geff4, estp(Geff4,EP[3]))

sed_aval, all_sed_aval = np.zeros(len(select1+select4), np.ndarray), np.zeros(len(select1+select4), np.ndarray)
for i in range(len(select1+select4)):
    k, sed_aval[i] = (select1+select4)[i], np.zeros(len(eventlen[i]), float)
    T_ratio = ceil(max(all_rain[k]))
    if i<3:
        G, V = Gst1, VP[0]
    else:
        G, V = Gst4, VP[3]
    S, current, branches, active, _, _ = simple_init(G)
    for run in range(30):
        run
        nb_rain = 0
        for T in eventlen[i]:
            #100*T/len(eventlen[i])
            nb_rain, exiting, p_rain = 0, 0, all_rain[k][T]/T_ratio
            while nb_rain < all_rain[k][T]:
                S, current, branches, active, nb_rain, exiting = simple_sandpilest_sedigraph_with_rain_intensity(G, V, S, current, branches, active, nb_rain, exiting, p_rain, not_sinks)
            sed_aval[i][T] += exiting
    sed_aval[i] /= 30

all_sed_mahl, all_sed_aval = np.zeros(8, np.ndarray), np.zeros(8, np.ndarray)
for i in range(len(select1+select4)):
    k = (select1+select4)[i]
    s_, r_, mahl_, aval_ = all_sed[k], all_rain[k], mahl_sed[i], sed_aval[i]
    r_avg, s_avg, mahl_avg, aval_avg = np.mean(r_[np.where(r_)]), np.mean(s_[np.where(s_)]), np.mean(mahl_[np.where(mahl_)]), np.mean(aval_[np.where(aval_)])
    all_sed[i], all_rain[i], all_sed_mahl[i], all_sed_aval[i] = s_*r_avg/s_avg, r_, mahl_*r_avg/mahl_avg, aval_*r_avg/aval_avg
    all_sed[i], all_rain[i], all_sed_mahl[i], all_sed_aval[i] = s_*max(r_)/max(s_), r_, mahl_*max(r_)/max(mahl_), aval_*max(r_)/max(aval_)

yields_exp, yields_mahl, yields_aval = np.array([yields[names1[(select1+select4)[i]]] for i in range(8)]), np.array([sum(all_sed_mahl[i]) for i in range(8)]), np.array([sum(all_sed_aval[i]) for i in range(8)])

names1[select1+select4]
yields_exp
yields_mahl

p1_aval = 






plt.plot([min(yields_exp[selected]),max(yields_exp[selected])], [min(yields_exp[selected]), max(yields_exp[selected])], linestyle='--', c='grey')
plt.scatter(yields_exp[selected], yields_exp[selected], c='grey', alpha=.8, s=40, label='empirical data')
plt.scatter(yields_exp[selected], yields_mahl[selected], alpha=.8, s=80, label='Mahleran')
plt.scatter(yields_exp[selected], yields_aval[selected], alpha=.8, s=80, label='avalanches')
for i in selected:
    _ = plt.text(yields_exp[i]-400, yields_exp[i]+100, s=names1[select1+select4][i], fontsize=12)

plt.xlabel('monitor sediment yield (g)')
plt.ylabel('modelled sediment yield (g)')
plt.legend()
c1,p1 = spearmanr(yields_exp[selected], yields_mahl[selected])
c2,p2 = spearmanr(yields_exp[selected], yields_aval[selected])
plt.title('Spearman correlation with Mahleran ($\\rho$='+str(round(c1,2))+', p-value='+str(round(p1,2))+') and Avalanches ($\\rho$='+str(round(c2,2))+', p-value='+str(round(p2,2))+')')
plt.show()

selected = [0,1,2,3,4,5,7]
[selected]

### 



# computing sediment yield ratios, mahleran over experiment, avalanches over experiment
point1_min, point1_max = min(yield1+yield1mahl+yield1aval), max(yield1+yield1mahl+yield1aval)
point4_min, point4_max = min(yield4+yield4mahl+yield4aval), max(yield4+yield4mahl+yield4aval)

x1_min, x1mahl_min, x1aval_min = min(yield1), min(yield1mahl), min(yield1aval)
x4_min, x4mahl_min, x4aval_min = min(yield4), min(yield4mahl), min(yield4aval)


plt.plot([point1_min,point1_min], [point1_min,point1_max], linestyle='--', c='grey')

plt.plot([min(yield1),max(yield1)], [min(yield1), max(yield1)], linestyle='--', c='grey')
plt.scatter(yield1, yield1, c='grey', label='monitored yield')
plt.scatter(yield1, yield1_mahl, label='modelled yield (Mahleran)')
plt.scatter(yield1, yield1_aval, label='modelled yield (avalanches)')
plt.ylabel('P1 sediment yield')
plt.title('Spearman correlations: '+str(round(spearmanr(yield1,yield1_mahl)[0],2))+' (with Mahleran), '+str(round(spearmanr(yield1,yield1_aval)[0],2))+' (with avalanches)')
plt.legend()
plt.figure()

# monitored yield on x-axis, modelled yield on y-axis

plt.plot([min(yield4),max(yield4)], [min(yield4), max(yield4)], linestyle='--', c='grey')
plt.scatter(yield4, yield4, c='grey', label='monitored yield')
plt.scatter(yield4, yield4_mahl, label='modelled yield (Mahleran)')
plt.scatter(yield4, yield4_aval, label='modelled yield (avalanches)')
plt.ylabel('P4 sediment yield')
plt.title('Spearman correlations: '+str(round(spearmanr(yield4,yield4_mahl)[0],2))+' (with Mahleran), '+str(round(spearmanr(yield4,yield4_aval)[0],2))+' (with avalanches)')
plt.legend()
plt.show()





fig, ax = plt.subplots(2,5,figsize=(19,7),tight_layout=True)
_ = [ax[i%2,i//2].plot(event_len[i], all_rain_p4[i], lw=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_rain_p4[0], lw=.7, label='rainfall amount')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4[i], lw=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4[0], lw=.7, label='sediment flow')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4mahl[i], linestyle='--', lw=2, alpha=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4mahl[0], linestyle='--', lw=2, alpha=.7, label='Mahleran sediment flow')

_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4aval[i], linestyle='-.', lw=2, alpha=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4aval[0], linestyle='-.', lw=2, alpha=.7, label='Avalanche sediment flow')

_ = [ax[i%2,i//2].set_title(names0[P4_loc0][i]) for i in range(10)]
_ = [ax[i%2,i//2].set_ylim(-5,135) for i in range(10)]
ax[0,0].legend()
ax[1,2].set_xlabel('time')
fig.suptitle('deterministic avalanches')
plt.show()

# (2) stochastic descent

sed_p4avaleff, all_sed_p4avaleff = np.zeros(10, np.ndarray), np.zeros(10, np.ndarray)
for i in range(10):
    k, sed_p4avaleff[i] = P4_loc1[i], np.zeros(len(event_len[i]), float)
    T_ratio = int(max(all_rain[k]))
    S, current, branches, active, _, _ = simple_init(Geff)
    for run in range(10):
        for T in event_len[i]:
            nb_rain, exiting, p_rain = 0, 0, all_rain[k][T]/T_ratio
            while nb_rain < all_rain[k][T]:
                S, current, branches, active, nb_rain, exiting = simple_sandpile_sedigraph_with_rain_intensity(Geff, VP[3], Ep, S, current, branches, active, nb_rain, exiting, p_rain, not_sinks)
            sed_p4avaleff[i][T] += exiting

for i in range(10):
    s_, r_, mahl_, aval_ = all_sed[P4_loc1][i], all_rain[P4_loc1][i], mahl_sed[i], sed_p4avaleff[i]
    r_avg, s_avg, mahl_avg, aval_avg = np.mean(r_[np.where(r_)]), np.mean(s_[np.where(s_)]), np.mean(mahl_[np.where(mahl_)]), np.mean(aval_[np.where(aval_)])
    all_sed_p4[i], all_rain_p4[i], all_sed_p4mahl[i], all_sed_p4avaleff[i] = s_*r_avg/s_avg, r_, mahl_*r_avg/mahl_avg, aval_*r_avg/aval_avg
    all_sed_p4[i], all_rain_p4[i], all_sed_p4mahl[i], all_sed_p4avaleff[i] = s_*max(r_)/max(s_), r_, mahl_*max(r_)/max(mahl_), aval_*max(r_)/max(aval_)

fig, ax = plt.subplots(2,5,figsize=(19,7),tight_layout=True)
_ = [ax[i%2,i//2].plot(event_len[i], all_rain_p4[i], lw=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_rain_p4[0], lw=.7, label='rainfall amount')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4[i], lw=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4[0], lw=.7, label='sediment flow')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4mahl[i], linestyle='--', lw=2, alpha=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4mahl[0], linestyle='--', lw=2, alpha=.7, label='Mahleran sediment flow')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4avaleff[i], linestyle='-.', lw=2, alpha=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4avaleff[0], linestyle='-.', lw=2, alpha=.7, label='Avalanche (stochastic)')
_ = [ax[i%2,i//2].set_title(names0[P4_loc0][i]) for i in range(10)]
_ = [ax[i%2,i//2].set_ylim(-5,135) for i in range(10)]
ax[0,0].legend()
ax[1,2].set_xlabel('time')
fig.suptitle('stochastic avalanches')
plt.show()


## sedigraph
fig, ax = plt.subplots(2,5,figsize=(19,7),tight_layout=True)
_ = [ax[i%2,i//2].plot(event_len[i], all_rain_p4[i]) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_rain_p4[0], label='rainfall amount')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4[i]) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4[0], label='sediment flow')
_ = [ax[i%2,i//2].plot(event_len[i], all_sed_p4mahl[i], linestyle='--', lw=2, alpha=.7) for i in range(1,10)]
_ = ax[0,0].plot(event_len[0], all_sed_p4mahl[0], linestyle='--', lw=2, alpha=.7, label='Mahleran sediment flow')
_ = [ax[i%2,i//2].set_title(names0[P4_loc0][i]) for i in range(10)]
_ = [ax[i%2,i//2].set_ylim(-5,235) for i in range(10)]
ax[0,0].legend()
ax[1,2].set_xlabel('time')
plt.show()



def simple_init(g):
    state, current_av, branches, new_active, size, dur = {}, [], {}, [], [], []
    for i in g:
        state[i] = 0
    return state, current_av, branches, new_active, size, dur


### SELECTING FROM plot_event_flow_data.xls

A = pd.ExcelFile('experimental_data_hydrograph/plot_event_flow_data.xls')
all_events = 40
data_tables = np.zeros(all_events, np.ndarray)
for i in range(all_events):
    data_tables[i] = A.parse(sheet_name=i).values

names = np.array([t[0,0] for t in data_tables])
select1, select4 = ['P1 07/09/2005', 'P1 29/08/2006', 'P1 07/09/2006'], ['P4 29/09/2005', 'P4 05/07/2006', 'P4 28/07/2006m', 'P4 23/08/2006', 'P4 29/08/2006']
names1, names4 = [int(np.where(names==e)[0][0]) for e in select1], [int(np.where(names==e)[0][0]) for e in select4]

rain_timeseries = np.zeros(10,list)
for i in range(10):
    i
    rain_col = data_tables[names[i]][1:,2].astype(str)
    names[i], len(rain_col[np.where(rain_col!='nan')[0][:-1]])
    rain_col = rain_col[np.where(rain_col!='nan')[0]]
    if i==0:
        rain1[i] = rain_col[:-2].astype(float)
    else:
        rain1[i] = rain_col[:-1].astype(float)

for i in range(4):
    i
    rain_col = data_tables[names4[i]][1:,2].astype(str)
    names4[i], select4[i], rain_col[np.where(rain_col!='nan')[0]]
    rain_col = rain_col[np.where(rain_col!='nan')[0]]
    if i==3:
        rain4[i] = rain_col[:-2].astype(float)
    else:
        rain4[i] = rain_col[:-1].astype(float)

rain = np.hstack((rain1,rain4))
with open('monitored_rainfall_data.txt','wb') as ll:
    pickle.dump(rain,ll)

rain = pickle.load(open('monitored_rainfall_data.txt','rb'))
yields_monitored = {'P1 07.09.05': 2709, 'P1 29.08.06': 1349, 'P1 07.09.06': 662, 'P2 07.09.05': 1197, 'P2 15.08.06': 36, 'P2 07.09.06': 608, 'P3 31.07.06': 1461, 'P3 01.08.06': 2086, 'P3 11.08.06': 503, 'P3 29.08.06': 294, 'P3 07.09.06': 2238, 'P4 07.09.05': 55, 'P4 29.09.05': 546, 'P4 05.07.06': 2114, 'P4 28.07.06a': 2405, 'P4 28.07.06b': 101, 'P4 31.07.06': 475, 'P4 01.08.06': 449, 'P4 11.08.06': 93, 'P4 23.08.06': 944, 'P4 29.08.06': 317, 'P4 28.07.06m': 328} #from only_sed => not good
# from hydrology2010
yields_monitored = {'P1 07/09/2005': 4.100, 'P1 29/08/2006': 2, 'P1 07/09/2006':.9, 'P4 05/07/2006':5.84, 'P4 28/07/2006':4.8, 'P4 23/08/2006':2.6}
yield_modelled = {'P1 07/09/2005': 4100, 'P1 29/08/2006', 'P1 07/09/2006'], ['P4 29/09/2005', 'P4 05/07/2006', 'P4 28/07/2006m', 'P4 23/08/2006', 'P4 29/08/2006'}

new_names = ['P1 07/09/2005', 'P1 05/07/2006', 'P1 01/08/2006', 'P1 29/08/2006', 'P1 07/09/2006', 'P4 07/09/2005', 'P4 05/07/2006', 'P4 28/07/2006m', 'P4 11/08/2006', 'P4 23/08/2006']
'P1 29/08/2006', 'P1 07/09/2006':.9, 'P4 05/07/2006':5.84, 'P4 28/07/2006':4.8, 'P4 23/08/2006':2.6}

eventlen = [len(rain[i]) for i in range(7)]


Geff1, Geff4 = eff_lattice(G,EP[0]), eff_lattice(G,EP[3])
Gst1, Gst4 = eff_steep(Geff1, estp(Geff1,EP[0])), eff_steep(Geff4, estp(Geff4,EP[3]))

sed_aval = np.zeros(7, np.ndarray)
for i in range(7):
    sed_aval[i] = np.zeros(eventlen[i], float)
    T_ratio = ceil(max(rain[i]))
    if i<3:
        G, V = Gst1, VP[0]
    else:
        G, V = Gst4, VP[3]
    S, current, branches, active, _, _ = simple_init(G)
    for run in range(50):
        run
        nb_rain = 0
        for T in range(eventlen[i]):
            nb_rain, exiting, p_rain = 0, 0, rain[i][T]/T_ratio
            while nb_rain < rain[i][T]:
                S, current, branches, active, nb_rain, exiting = simple_sandpilest_sedigraph_with_rain_intensity(G, V, S, current, branches, active, nb_rain, exiting, p_rain, not_sinks)
            sed_aval[i][T] += exiting
    sed_aval[i] /= 50

[float(round(sum(sed_aval[i]),2)) for i in range(7)]


all_sed_aval = np.zeros(7, np.ndarray)
for i in range(len(select1+select4)):
    k = (select1+select4)[i]
    s_, r_, mahl_, aval_ = all_sed[k], all_rain[k], mahl_sed[i], sed_aval[i]
    r_avg, s_avg, mahl_avg, aval_avg = np.mean(r_[np.where(r_)]), np.mean(s_[np.where(s_)]), np.mean(mahl_[np.where(mahl_)]), np.mean(aval_[np.where(aval_)])
    all_sed[i], all_rain[i], all_sed_mahl[i], all_sed_aval[i] = s_*r_avg/s_avg, r_, mahl_*r_avg/mahl_avg, aval_*r_avg/aval_avg
    all_sed[i], all_rain[i], all_sed_mahl[i], all_sed_aval[i] = s_*max(r_)/max(s_), r_, mahl_*max(r_)/max(mahl_), aval_*max(r_)/max(aval_)

yields_exp, yields_mahl, yields_aval = np.array([yields[names1[(select1+select4)[i]]] for i in range(8)]), np.array([sum(all_sed_mahl[i]) for i in range(8)]), np.array([sum(all_sed_aval[i]) for i in range(8)])


mahl_sed, mahl_hyd = np.zeros(len(select1+select4), np.ndarray), np.zeros(len(select1+select4), np.ndarray)
for i in range(len(names1+select4)):
    k = (names1+select4)[i]
    sedtr, hydro = np.loadtxt('mahleran_outputs/'+folder[i]+'/sedtr001.dat', float), np.loadtxt('mahleran_outputs/'+folder[i]+'/hydro001.dat', float)
    sed = sedtr[:,1][:np.where(sedtr[:,1])[0][-1]]
    hyd = hydro[:len(sed)+1,1]
    (select1+select4)[k], sum(sed), len(sed)
    #mahl_sed[i], mahl_hyd[i] = np.array([np.mean(sed[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)]), np.array([np.mean(hyd[t*60:(t+1)*60]) for t in range(int(len(sed)/60)+1)])



