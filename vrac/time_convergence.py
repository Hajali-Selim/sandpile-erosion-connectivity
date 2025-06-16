from avalanche_packages import *
from inits import *
from scipy.stats import zscore

G, not_sinks = tgrid(20,30)
DSCHG, SEDTR = np.zeros(4, dtype=np.ndarray), np.zeros(4, dtype=np.ndarray)
for k in range(4):
    DSCHG[k], SEDTR[k] = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'_rainA_highsm_dschg.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'_rainA_highsm_sedtr.asc', skiprows=6)[1:-1,1:-1]

vP, eP, VP, EP, DSCHGP, SEDTRP = np.zeros(4, np.ndarray), np.zeros(4, np.ndarray), np.zeros(4, dict), np.zeros(4, dict), np.zeros(4, dict), np.zeros(4, dict)
for k in range(4):
    vP[k], eP[k], VP[k], EP[k], DSCHGP[k], SEDTRP[k] = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'vegcover.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'dem.asc', skiprows=6)[1:-1,1:-1], {}, {}, {}, {}
    for i in range(60):
        for j in range(20):
            VP[k][(j,59-i)], EP[k][(j,59-i)], DSCHGP[k][(j,59-i)], SEDTRP[k][(j,59-i)] = vP[k][i,j]/100, eP[k][i,j], DSCHG[k][i,j], SEDTR[k][i,j]

Geff = [eff_lattice(G,EP[0])]+[0,0]+[eff_lattice(G,EP[3])]
Gst = [eff_steep(Geff[0],estp(Geff[0],EP[0]))]+[0,0]+[eff_steep(Geff[3],estp(Geff[3],EP[3]))]
Ceff = [init(Geff[0])[1]]+[0,0]+[init(Geff[3])[1]]
Cst = [init(Gst[0])[1]]+[0,0]+[init(Gst[3])[1]]
Epm_eff = [eprob_map(Geff[0],EP[0])]+[0,0]+[eprob_map(Geff[3],EP[3])]
#Epm_st = [eprob_map(Gst[k],EP[k]) for k in range(4)]
Epm_st = [{(i,j):1 for i,j in Gst[0].edges}]+[0,0]+[{(i,j):1 for i,j in Gst[3].edges}]
SCE_eff = [SC_EDGE(Geff[0],Epm_eff[0])]+[0,0]+[SC_EDGE(Geff[3],Epm_eff[3])]
SCE_st = [SC_EDGE(Gst[0],Epm_st[0])]+[0,0]+[SC_EDGE(Gst[3],Epm_st[3])]

Ep = [eprob(Geff[0],EP[0])]+[0,0]+[eprob(Geff[3],EP[3])]


Ceff,Cst,SC_st,SC_eff,Ep,Epm_eff,Epm_st = pickle.load(open('time_convergence_data.txt', 'rb'))
with open('time_convergence_data.txt','wb') as ll:
    pickle.dump((Ceff,Cst,SC_st,SC_eff,Ep,Epm_eff,Epm_st),ll)



plt.figure(figsize=(4,11))
nx.draw(G, dict(zip(G,G)), node_size=150, node_color=list(sc0.values()), edgelist=[])
plt.show()

with open('time_convergence_data.txt','wb') as ll:
    pickle.dump((Ceff,Cst,SC_st,SC_eff,Ep,Epm_eff,Epm_st),ll)

Ceff,Cst,SC_st,SC_eff,Ep,Epm_eff,Epm_st = pickle.load(open('time_convergence_data.txt','rb'))

Ep, Epm_eff, Epm_st = np.array([eprob(Geff[k],EP[k]) for k in range(4)]), np.array([eprob_map(Geff[k],EP[k]) for k in range(4)]), np.array([eprob_map(Gst[k],EP[k]) for k in range(4)])
#SCE_eff, SCE_st = np.array([SC_EDGE(Geff[k],Ceff[k],Epm_eff[k]) for k in range(4)]), np.array([SC_EDGE(Gst[k],Cst[k],Epm_st[k]) for k in range(4)])
SCE_eff = np.array([SC_EDGE(Geff[k],Epm_eff[k]) for k in range(4)])
#SCP_eff, SCP_st = [SC_PATH(Geff[k],Ceff[k],Epm_eff[k]) for k in range(3)]+[SCP], [SC_PATH(Gst[k],Cst[k],Epm_st[k]) for k in range(4)]

SCP_eff = pickle.load(open('time_convergence_data_previous.txt','rb'))[3][[0,3]]


cnew_from_init = init(gnew)
pnew = eprob(gnew, eP[3])


nb_runs, nb_steps, step_slice = 100, 200000, 1000
Sizes_eff, SCpFCl_eff, SCeFCv_eff = np.zeros(2, np.ndarray), np.zeros(2, np.ndarray), np.zeros(2, np.ndarray)
SCeFCl_eff, SCpFCv_eff = np.zeros(2, np.ndarray), np.zeros(2, np.ndarray)
FCl_eff, FCv_eff = np.zeros(2, np.ndarray), np.zeros(2, np.ndarray)
for k in range(1,2):
    print('LANDSCAPE N°',3*k+1,'...')
    SCeFCv_eff[k] = np.zeros((nb_runs,int(nb_steps/step_slice)), float)
    Sizes_eff[k], FCv_eff[k] = np.zeros(nb_runs, list), {node:0 for node in G}
    
    for run in range(1,nb_runs):
        print('.. run n°',run)
        #S, C, current, branches, active, Sizes_eff[k][run] = zeros(Ceff[k*3])
        S, C, current, branches, active, _ = zeros(Ctest)
        for step in range(nb_steps+1):
            S, C, current, branches, active, _ = sandpile_new_form(Geff[3], VP[3], eprob_test, S, C, current, branches, active, _,not_sinks)
            if step % step_slice == 0 and step>0:
                FCl, FCv = FC_path(C,'len'), FC_path(C,'var')
                SCpFCl_eff[k][run, int(step/step_slice-1)], SCeFCv_eff[k][run, int(step/step_slice-1)] = SCFC(SCP_eff[k],FCl,not_sinks), SCFC(SCE_eff[k],FCv,not_sinks)
                SCeFCl_eff[k][run, int(step/step_slice-1)], SCpFCv_eff[k][run, int(step/step_slice-1)] = SCFC(SCE_eff[k],FCl,not_sinks), SCFC(SCP_eff[k],FCv,not_sinks)
        for node in G:
           FCl_eff[k][node] += FCl[node]
           FCv_eff[k][node] += FCv[node]

SCeFCv_eff_new = np.zeros(2, np.ndarray)
SCeFCv_eff_new[0] = np.zeros((nb_runs,int(nb_steps/step_slice)), float)
SCeFCv_eff_new[1] = np.zeros((nb_runs,int(nb_steps/step_slice)), float)
for run in range(nb_runs):
    for ti in range(100):
        SCeFCv_eff_new[0][run,ti] = SCFC(SCE_eff[k],FCv,not_sinks)
        SCFC(SCE_eff[0],FCt[0,run,ti],not_sinks)

k=0
plt.fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCv_eff[k][:,i])-np.std(SCeFCv_eff[k][:,i]) for i in range(100)], [np.mean(SCeFCv_eff[k][:,i])+np.std(SCeFCv_eff[k][:,i]) for i in range(100)], alpha=.7, label='SCe-FCv')
plt.fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCl_eff[k][:,i])-np.std(SCpFCl_eff[k][:,i]) for i in range(100)], [np.mean(SCpFCl_eff[k][:,i])+np.std(SCpFCl_eff[k][:,i]) for i in range(100)], alpha=.7, label='SCp-FCl')
plt.legend()
plt.figure()

plt.fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCv_eff[k][:,i])-np.std(SCpFCv_eff[k][:,i]) for i in range(100)], [np.mean(SCpFCv_eff[k][:,i])+np.std(SCpFCv_eff[k][:,i]) for i in range(100)], alpha=.7, label='SCp-FCv')
plt.fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCl_eff[k][:,i])-np.std(SCeFCl_eff[k][:,i]) for i in range(100)], [np.mean(SCeFCl_eff[k][:,i])+np.std(SCeFCl_eff[k][:,i]) for i in range(100)], alpha=.7, label='SCe-FCl')
plt.legend()
plt.show()

## FCv (FC) figure
k=0
fig,ax=plt.subplots(1,2,figsize=(10,4),tight_layout=True)
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCv_eff[k][:,i])-np.std(SCeFCv_eff[k][:,i]) for i in range(100)], [np.mean(SCeFCv_eff[k][:,i])+np.std(SCeFCv_eff[k][:,i]) for i in range(100)], alpha=.8, label='grassland, SC$_E$-FC(t)', color='forestgreen')
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCl_eff[k][:,i])-np.std(SCpFCl_eff[k][:,i]) for i in range(100)], [np.mean(SCpFCl_eff[k][:,i])+np.std(SCpFCl_eff[k][:,i]) for i in range(100)], alpha=.8, label='grassland, SC$_P$-DP(t)', color='lightgreen')
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCv_eff[k][:,i]) for i in range(100)], lw=2, c='lime', alpha=.6)
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCl_eff[k][:,i]) for i in range(100)], lw=2, c='honeydew', alpha=.6)
k=1
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCv_eff[k][:,i])-np.std(SCeFCv_eff[k][:,i]) for i in range(100)], [np.mean(SCeFCv_eff[k][:,i])+np.std(SCeFCv_eff[k][:,i]) for i in range(100)], alpha=.8, label='shrubland, SC$_E$-FC(t)', color='darkorange')
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCl_eff[k][:,i])-np.std(SCpFCl_eff[k][:,i]) for i in range(100)], [np.mean(SCpFCl_eff[k][:,i])+np.std(SCpFCl_eff[k][:,i]) for i in range(100)], alpha=.8, label='shrubland, SC$_P$-DP(t)', color='orange')
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCeFCv_eff[k][:,i]) for i in range(100)], lw=2, c='bisque', alpha=.6)
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCpFCl_eff[k][:,i]) for i in range(100)], lw=2, c='yellow', alpha=.6)
ax[0].legend(fontsize=12)
ax[0].set_xlabel('step (10$^3$)', size=12)
ax[0].set_ylabel('correlation', size=12)
ax[0].set_xticks([5000,100000,200000,300000,400000,500000], labels=[5,100,200,300,400,500])

ax[1].scatter(u1_k0, u2_k0, marker='x', s=30, c='forestgreen', label='grassland')
ax[1].scatter(u1_k3, u2_k3, marker='x', s=30, c='darkorange', label='shrubland')
ax[1].loglog([])
ax[1].set_xlabel('avalanche size', size=12)
ax[1].set_ylabel('count', size=12)
ax[1].legend(fontsize=12)
ax[1].set_title(' (b)', loc='left')
plt.show()

with open('time_convergence_stochastic_dyn.txt','wb') as ll:
    pickle.dump((FCl_eff,FCv_eff, SCpFCl_eff,SCeFCv_eff,SCeFCl_eff,SCpFCv_eff, Sizes_eff),ll)

u1_k0, u2_k0 = np.unique(sum(Sizes[0],[]), return_counts=True)
u1_k3, u2_k3 = np.unique(sum(Sizes[1],[]), return_counts=True)

fig, ax = plt.subplots(1,2,tight_layout=True,figsize=(10,4))
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCFCv[0][:,i])-np.std(SCFCv[3][:,i]) for i in range(100)], [np.mean(SCFCv[3][:,i])+np.std(SCFCv[3][:,i]) for i in range(100)], alpha=.7, label='shrubland')
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCFCv[3][:,i]) for i in range(100)], c='cyan', lw=2)
ax[0].fill_between(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCFCv[0][:,i])-np.std(SCFCv[0][:,i]) for i in range(100)], [np.mean(SCFCv[0][:,i])+np.std(SCFCv[0][:,i]) for i in range(100)], alpha=.7, label='grassland')
ax[0].plot(np.arange(step_slice, nb_steps+1, step_slice), [np.mean(SCFCv[0][:,i]) for i in range(100)], c='bisque', lw=2)
ax[0].legend(fontsize=12)
ax[0].set_xticks(ticks=[step_slice]+list(np.arange(step_slice*20, nb_steps+1, step_slice*20)), labels=[2,40,80,120,160,200])
ax[0].set_yticks(ticks=[.2,.3,.4,.5,.6,.7,.8,.9])
ax[0].set_xlabel('step (10$^3$)', size=12)
ax[0].set_ylabel('SC$_E$-FC correlation', size=12)
ax[0].set_title(' (a)', loc='left')
ax[1].scatter(u1_k3, u2_k3, marker='x', s=30, c='tab:blue', label='shrubland')
ax[1].scatter(u1_k0, u2_k0, marker='x', s=30, c='tab:orange', label='grassland')
ax[1].loglog([])
ax[1].set_xlabel('avalanche size', size=12)
ax[1].set_ylabel('count', size=12)
ax[1].legend(fontsize=12)
ax[1].set_title(' (b)', loc='left')
plt.show()

### ### ### ### TIME-CONVERGENCE OF DETERMINISTIC DYNAMICS ON P4

nb_runs, nb_steps, step_slice = 100, 10000, 100
# SC vs Mahleran FC # SC vs sandpile FC # sandpile FC vs Mahleran FC
# SCFCh, SCFCsd → not a function of time, not direct numbers

for k in range(4):
    print('P_'+str(k)+':')
    r1, r2 = SCFC(SC_st[k],DSCHGP[k],not_sinks), SCFC(SC_st[k],SEDTRP[k],not_sinks)
    print('Hydrological connectivity correlates at:',round(r1,2))
    print('Sediment connectivity correlates at:',round(r2,2))

A, B = DSCHG[k], np.array([[SC_st[k][(j,59-i)] for j in range(20)] for i in range(60)])
fig, ax = plt.subplots(1,2,tight_layout=True)
ax[0].imshow(A, cmap='rainbow')
ax[1].imshow(B, cmap='rainbow')
_ = [ax[i].set_xticks([]) for i in range(2)]
_ = [ax[i].set_yticks([]) for i in range(2)]
plt.show()


nb_runs, nb_steps, step_slice = 100, 10000, 200
SCFCt, FCtFCwa, FCtFCsd, FCv_st = np.zeros(4, np.ndarray), np.zeros(4, np.ndarray), np.zeros(4, np.ndarray), np.zeros(4, dict)
for k in [3]:#range(4):
    SCFCt[k], FCtFCwa[k], FCtFCsd[k], FCv_st[k] = np.zeros((nb_runs,int(nb_steps/step_slice)), float), np.zeros((nb_runs,int(nb_steps/step_slice)), float), np.zeros((nb_runs,int(nb_steps/step_slice)), float), {node:0 for node in G}
    for run in range(nb_runs):
        print('run n°',run)
        S, C, current, branches, active, _ = zeros(Cst[k])
        for step in range(nb_steps+1):
            S, C, current, branches, active, _ = sandpilest(Gst[k], VP[k], S, C, current, branches, active, _, not_sinks)
            if step % step_slice == 0 and step>0:
                FCv = FC_path(C,'var')
                SCFCt[k][run, int(step/step_slice-1)], FCtFCwa[k][run, int(step/step_slice-1)], FCtFCsd[k][run, int(step/step_slice-1)] = SCFC(SC_st[k],FCv,not_sinks), SCFC(DSCHGP[k],FCv,not_sinks), SCFC(SEDTRP[k],FCv,not_sinks)
        for node in G:
            FCv_st[k][node] += FCv[node]

#k=0 running for 60000/500 steps

SCEtest_eff = [SC_EDGE(Geff[k],Ceff[k],Epm_eff[k]) for k in range(4)]
SCEtest_st = [SC_EDGE(Gst[k],Cst[k],Epm_st[k]) for k in range(4)]

SCFCt, FCtFCwa, FCtFCsd, FCt_st = pickle.load(open('time_convergence_deterministic_dyn.txt','rb'))

k = 3
plt.figure(tight_layout=True)
SCFCt_avg, FCtFCwa_avg, FCtFCsd_avg = np.array([np.mean(SCFCt[k][:,t]) for t in range(len(SCFCt[k][0]))]), np.array([np.mean(FCtFCwa[k][:,t]) for t in range(len(SCFCt[k][0]))]), np.array([np.mean(FCtFCsd[k][:,t]) for t in range(len(SCFCt[k][0]))])
SCFCt_std, FCtFCwa_std, FCtFCsd_std = np.array([np.std(SCFCt[k][:,t]) for t in range(len(SCFCt[k][0]))]), np.array([np.std(FCtFCwa[k][:,t]) for t in range(len(SCFCt[k][0]))]), np.array([np.std(FCtFCsd[k][:,t]) for t in range(len(SCFCt[k][0]))])
plt.fill_between(np.arange(step_slice,nb_steps+1,step_slice), SCFCt_avg-SCFCt_std, SCFCt_avg+SCFCt_std, label='FC-SC$_L$: sandpile FC and link-based SC', color='tab:blue', alpha=.6, zorder=2)
#plt.fill_between(np.arange(step_slice,nb_steps+1,step_slice), FCtFCwa_avg-FCtFCwa_std, FCtFCwa_avg+FCtFCwa_std, label='FC-FC$_w$: sandpile FC and MAHLERAN water connectivity', color='tab:orange', alpha=.7)
plt.fill_between(np.arange(step_slice,nb_steps+1,step_slice), FCtFCsd_avg-FCtFCsd_std, FCtFCsd_avg+FCtFCsd_std, label='FC-FC$_s$: sandpile FC and MAHLERAN sediment connectivity', color='darkorange', alpha=.8)
#plt.plot(np.arange(step_slice,nb_steps+1,step_slice), FCtFCwa_avg, c='bisque', lw=2)
plt.plot(np.arange(step_slice,nb_steps+1,step_slice), SCFCt_avg, c='cyan', lw=2)
plt.plot(np.arange(step_slice,nb_steps+1,step_slice), FCtFCsd_avg, c='bisque', lw=2)
plt.xlabel('step (10$^3$)', size=12)
plt.ylabel('correlation', size=12)
plt.xticks([200,2000,4000,6000,8000,10000], labels=[.2,2,4,6,8,10])
plt.legend(loc='lower right', fontsize=12)
plt.xlim(0, 10200)
#plt.xlim(right=10000)
plt.show()

fig, ax = plt.subplots(1,5,figsize=(10,6),tight_layout=True)
p = [ax[0].imshow(vP[3]/100, cmap='Greens'), ax[1].imshow(nx_to_mat(SCE_st[3]), cmap='terrain_r'), ax[2].imshow(nx_to_mat(FCt_st[3]), cmap='terrain_r'), ax[3].imshow(DSCHG[3], cmap='terrain_r'), ax[4].imshow(SEDTR[3], cmap='terrain_r')]
_ = [fig.colorbar(p[i],ax=ax[i],fraction=.015,pad=.008,orientation='horizontal') for i in range(5)]
_ = [ax[i].set_xticks([]) for i in range(5)]
_ = [ax[i].set_yticks([]) for i in range(5)]
ax[0].set_title(' (a)', loc='left', size=12)
ax[1].set_title(' (b)', loc='left', size=12)
ax[2].set_title(' (c)', loc='left', size=12)
ax[3].set_title(' (d)', loc='left', size=12)
ax[4].set_title(' (e)', loc='left', size=12)
plt.subplot_tool()
plt.show()

fig, ax = plt.subplots(1,5,figsize=(10,6),tight_layout=True)
#p = [ax[0].imshow(vP[0]/100, cmap='Greens'), ax[1].imshow(nx_to_mat(SCE_st[0]), cmap='terrain_r'), ax[2].imshow(nx_to_mat(FCt_st[0]), cmap='terrain_r'), ax[3].imshow(DSCHG[0], cmap='terrain_r'), ax[4].imshow(SEDTR[0], cmap='terrain_r')]
p = [ax[0].imshow(vP[0]/100, cmap='Greens'), ax[1].imshow(nx_to_mat(sc0), cmap='terrain_r'), ax[2].imshow(nx_to_mat(FCtest_st), cmap='terrain_r'), ax[3].imshow(DSCHG[0], cmap='terrain_r'), ax[4].imshow(SEDTR[0], cmap='terrain_r')]
_ = [fig.colorbar(p[i],ax=ax[i],fraction=.015,pad=.008,orientation='horizontal') for i in range(5)]
_ = [ax[i].set_xticks([]) for i in range(5)]
_ = [ax[i].set_yticks([]) for i in range(5)]
ax[0].set_title(' (a)', loc='left', size=12)
ax[1].set_title(' (b)', loc='left', size=12)
ax[2].set_title(' (c)', loc='left', size=12)
ax[3].set_title(' (d)', loc='left', size=12)
ax[4].set_title(' (e)', loc='left', size=12)
plt.subplot_tool()
plt.show()

plt.imshow(np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'_rainA_highsm_dschg.asc', skiprows=6), cmap='terrain_r')


plt.figure(figsize=(3,9))
nx.draw(Gst[0],pos=dict(zip(G,G)), edgelist=[], node_size=80,node_color=list(SCE_st[0].values()))

ax[4].fill_between(np.arange(step_slice,nb_steps+1,step_slice), SCFCt_avg-SCFCt_std, SCFCt_avg+SCFCt_std, label='SC-FC(t)', color='grey', alpha=.8)
ax[4].fill_between(np.arange(step_slice,nb_steps+1,step_slice), FCtFCwa_avg-FCtFCwa_std, FCtFCwa_avg+FCtFCwa_std, label='FC(t)-FC$^*$: MAHLERAN\nwater connectivity', color='tab:blue', alpha=.8)
ax[4].fill_between(np.arange(step_slice,nb_steps+1,step_slice), FCtFCsd_avg-FCtFCsd_std, FCtFCsd_avg+FCtFCsd_std, label='FC(t)-FC$^*$: MAHLERAN\nsediment connectivity', color='tab:orange', alpha=.8)
ax[4].plot(np.arange(step_slice,nb_steps+1,step_slice), SCFCt_avg, c='silver', lw=2)
ax[4].plot(np.arange(step_slice,nb_steps+1,step_slice), FCtFCwa_avg, c='cyan', lw=2)
ax[4].plot(np.arange(step_slice,nb_steps+1,step_slice), FCtFCsd_avg, c='bisque', lw=2)
ax[4].set_xlabel('step (10$^3$)', size=12)
ax[4].set_ylabel('correlation', size=12)
ax[4].set_xticks([1000,10000,20000,30000,40000,50000], labels=[1,10,20,30,40,50])
ax[4].legend(loc='lower right')
_ = [ax[i].set_xticks([]) for i in range(4)]
_ = [ax[i].set_yticks([]) for i in range(4)]
plt.show()


# stochastic descent p4

SCFCt, SCwFCt, FCt = np.zeros(4, np.ndarray), np.zeros(4, np.ndarray), np.zeros(4, dict)
with open('time_convergence_stochastic_dyn.txt','wb') as ll:
    pickle.dump((SCFCv, Sizes, FC_final),ll)

SCFCv, Sizes, FC_final = pickle.load(open('time_convergence_stochastic_dyn.txt','rb'))

fig,ax=plt.subplots(1,2,figsize=(6,7),tight_layout=True)
p=ax[0].imshow(np.log2(1+nx_to_mat(FC_final[k])), cmap=cmap)
fig.colorbar(p,ax=ax[0],shrink=.8)
p=ax[1].imshow(np.log2(1+nx_to_mat(sce_eff)), cmap=cmap)
fig.colorbar(p,ax=ax[1],shrink=.8)
plt.show()

fig,ax=plt.subplots(1,2,figsize=(6,7),tight_layout=True)
p=ax[0].imshow(nx_to_mat(FC_final[k]), cmap=cmap)
fig.colorbar(p,ax=ax[0],shrink=.8)
p=ax[1].imshow(nx_to_mat(sce_eff), cmap=cmap)
fig.colorbar(p,ax=ax[1],shrink=.8)
plt.show()

nb_runs, nb_steps, step_slice = 100, 200000, 2000
for k in [3]:#range(4):
    print('P'+str(k)+':')
    SCFCt[k], SCwFCt[k], FCt[k] = np.zeros((nb_runs,int(nb_steps/step_slice)), float), np.zeros((nb_runs,int(nb_steps/step_slice)), float), {node:0 for node in G}
    for run in range(nb_runs):
        print('run n°',run)
        S, C, current, branches, active, _ = zeros(Ceff[k])
        for step in range(nb_steps+1):
            S, C, current, branches, active, _, pexit = sandpile(Geff[k], VP[k], Ep[k], S, C, current, branches, active, [], 0, not_sinks)
            if step % step_slice == 0 and step>0:
                fc = FC_path(C,'var')
                SCFCt[k][run, int(step/step_slice-1)], SCwFCt[k][run, int(step/step_slice-1)] = SCFC(SC_eff[k],fc,not_sinks), SCFC(SCw_eff[k],fc,not_sinks)
        for node in G:
            FCt[k][node] += fc[node]

with open('time_convergence_stochastic_dyn.txt','wb') as ll:
    pickle.dump((SCFCt, SCwFCt),ll)

### VARYING DISSIPATION COEFFICIENT
D, SCFCD, ExitD = np.arange(.1,1.01,.1), np.zeros((10,nb_runs), np.ndarray), np.zeros((10,nb_runs), np.ndarray)
nb_steps, k = 10000, 3
for d_idx in range(len(D)):
    d=D[d_idx]
    d
    for run in range(nb_runs):
        print('run n°',run)
        S, C, current, branches, active, _ = zeros(Cst[k])
        ExitD[d_idx,run] = 0
        for step in range(nb_steps):
            S, C, current, branches, active, _, ExitD[d_idx,run] = sandpilest_vary_dissipation(Gst[k], VP[k], S, C, current, branches, active, _, ExitD[d_idx,run], not_sinks, d)
        FC = FC_path(C,'var')
        SCFCD[d_idx,run] = SCFC(SC_st[k],FC,not_sinks)

SCFCD_avg, SCFCD_std = np.array([np.mean(SCFCD[i]) for i in range(10)]), np.array([np.std(SCFCD[i]) for i in range(10)])

plt.figure(tight_layout=True)
plt.fill_between(D, SCFCD_avg-SCFCD_std, SCFCD_avg+SCFCD_std, alpha=.8)
plt.plot(D, SCFCD_avg, lw=2, color='cyan')
plt.xlabel('dissipation coefficient')
plt.ylabel('SC-FC correlation')
plt.show()

ExitD_avg, ExitD_std = np.array([np.mean(ExitD[i]) for i in range(10)]), np.array([np.std(ExitD[i]) for i in range(10)])
plt.fill_between(D, ExitD_avg-ExitD_std, ExitD_avg+ExitD_std, alpha=.8)
plt.plot(D, ExitD_avg, lw=2, color='cyan')
plt.show()

with open('SCFC_dissipation_p4.txt','wb') as ll:
    pickle.dump(SCFCD, ll)

### RECONSTRUCTION SCATTERPLOT

fig,ax=plt.subplots(1,2,figsize=(10.5,4),tight_layout=True)
ax[0].scatter(land_scores_0, MTmat[0], marker='x', s=25, label='grassland')
ax[0].plot([0,5], [beta_0,alpha_0*5+beta_0], lw=2, c='r')
ax[0].text(4.45,.013,s=str(round(alpha_0,4)*1e3)+' 10$^{-3}$', size=12, c='r', weight='bold')
ax[0].set_xlabel('surrounding vegetation density',size=12)
ax[0].set_ylabel('microtopography',size=12)
ax[0].legend(loc='lower right', fontsize=12,framealpha=1)
ax[1].scatter(land_scores_3, MTmat[3], marker='x', s=25, label='shrubland')
ax[1].plot([0,5], [beta_3,alpha_3*5+beta_3], lw=2, c='r')
ax[1].text(4.2,-.003,s=str(round(alpha_3,4)*1e3)+' 10$^{-3}$', size=12, c='r', weight='bold')
ax[1].set_xlabel('surrounding vegetation density',size=12)
ax[1].set_ylabel('microtopography',size=12)
ax[0].set_ylim(-.07,.07)
ax[1].set_ylim(-.07,.07)
ax[0].set_title(' (a)', loc='left')
ax[1].set_title(' (b)', loc='left')
ax[1].legend(loc='lower right', fontsize=12,framealpha=1)
plt.show()

### 




