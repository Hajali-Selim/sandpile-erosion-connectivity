import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from operator import mul
from functools import reduce
from scipy.stats import pearsonr
import pickle, random, itertools, sys

### SECTION 1

def tgrid(x, y):
    ## Generate grid of width x and length y/2 (by default, set x=20, y=30)
    g = nx.grid_2d_graph(x,int(2*y))
    g.remove_edges_from([((i1,0),(i1+1,0)) for i1 in range(x-1)])
    g = g.to_directed()
    g.remove_edges_from([((i1,i2),(i1,i2+1)) for i1,i2 in itertools.product(range(x), range(int(2*y)-1))])
    g.add_edges_from([((0,i2),(x-1,i2)) for i2 in range(1,int(2*y))]+[((x-1,i2),(0,i2)) for i2 in range(1,int(2*y))])
    s_ = [(i, 0) for i in range(x)]
    return g, [i for i in g if i not in s_]

## example
x1, x2 = 20, 30
G, not_sinks = tgrid(x1,x2) # not_sinks are all the nodes minus the lowest layer, it is useful to other functions

### SECTION 2

vi = 1
def tclustering_square(g, cover, pclust, vi):
    ## Populate grid with vegetation of value 1 (see th line to change this to randomly select v_i)
    # cover: ratio of nodes with v_i != 0
    # pclust: vegetation clustering probability
    v = {}
    for i in g:
        v[i] = 0
    free_nodes, nb_veg_final, nb_veg = set(list(g)), int(len(g)*cover), 0
    node = random.choice(list(free_nodes))
    while nb_veg < nb_veg_final:
        v[node] = vi#random.random()
        free_nodes.remove(node)
        nb_veg += 1
        if random.random() < pclust:
            i,j = node
            free_neighbors = {(i,j+1), (i,j-1), (i-1,j), (i+1,j), (i-1,j+1), (i+1,j+1), (i+1,j-1), (i-1,j+1)} & free_nodes
            if free_neighbors:
                node = random.choice(list(free_neighbors))
            else:
                node = random.choice(list(free_nodes))
        else:
            node = random.choice(list(free_nodes))
    return v

## example
vcover, pc = .3, .5
V = tclustering_square(G, vcover, pc, vi)
plt.figure(figsize=(3,9))
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(V.values()), cmap='Greens')
plt.show()

### SECTION 3: SAMPLING SHRUBLAND

## Importing data from natural landscapes
# Set: k=0 (grassland), k=1 (grass-shrubland), k=2 (shrub-grassland), k=3 (shrubland)
# the term 'int(k==3)' (below) is added because the file k=3 doesn't have the same dimensions as the others
k = 3
VD, ED, Vp, Ep, MT = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'vegcover.asc', skiprows=6)[1+int(k==3):-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'dem.asc', skiprows=6)[1:-1,1:-1], {}, {}, {}
for i,j in G:
    Vp[i,j], Ep[i,j] = VD[59-j,i], ED[59-j,i]
    MT[i,j] = Ep[i,j] - np.mean(ED[59-j])

## Computing scores from chosen natural landscape, the score of a node i = v_i + sum(v_j) where j are the 4 neighbors of i
land_scores = np.zeros((60,20), dtype=float)
for j in range(1,19):
    for i in range(1,59):
        land_scores[i,j] = 1*VD[i,j] + VD[i-1,j] + VD[i+1,j] + VD[i,j-1]+ VD[i,j+1]
        land_scores[i,0] = 2*VD[i,0] + VD[i-1,0] + VD[i+1,0] + VD[i,1]
        land_scores[i,19] = 2*VD[i,19] + VD[i-1,19] + VD[i+1,19] + VD[i,18]
    land_scores[0,j] = 2*VD[0,j] + VD[0,j-1] + VD[0,j+1] + VD[1,j]
    land_scores[59,j] = 2*VD[59,j] + VD[59,j-1] + VD[59,j+1] + VD[58,j]

land_scores[0,0] = 3*VD[0,0] + VD[0,1] + VD[1,0]
land_scores[59,0] = 3*VD[59,0] + VD[59,1] + VD[58,0]
land_scores[0,19] = 3*VD[0,19] + VD[0,18] + VD[1,19]
land_scores[59,19] = 3*VD[59,19] + VD[58,19] + VD[59,18]
land_scores /= 100

land1 = land_scores
land4 = land_scores

### SECTION 4: GENERATING SHRUBLAND

## Linear regression between measured scores vs microtopography of each node
X, Y = [land_scores[59-j,i] for i,j in G], [Ep[i,j]-np.mean([Ep[l,j] for l in range(20)]) for i,j in G]
alpha, beta = np.polyfit(X,Y,1)
Y_line = [i*alpha+beta for i in X]
alpha_std = np.std([Y[k]-Y_line[k] for k in range(len(G))])

def generate_elevation_extended(g, v, ep, a, a_std, b):
    v_score, mt, e, g_x, g_y, plane = {}, {}, {}, list(g)[-1][0], list(g)[-1][1], np.sum(ep,axis=1)/20
    print(plane)
    extended_plane = sum([list(plane+k*(1-.8*plane[-1])) for k in reversed(range(int(g_y/60)+1))],[])
    print(extended_plane)
    for i in range(1,g_x):
        for j in range(1,g_y):
            v_score[i,j] = v[i,j] + v[i-1,j]+ v[i+1,j] + v[i,j-1] + v[i,j+1]
            v_score[0,j] = 2*v[0,j] + v[1,j] + v[0,j-1] + v[0,j+1]
            v_score[g_x,j] = 2*v[g_x,j] + v[g_x-1,j] + v[g_x,j-1] + v[g_x,j+1]
        v_score[i,0] = 2*v[i,0] + v[i-1,0] + v[i+1,0] + v[i,1]
        v_score[i,g_y] = 2*v[i,g_y] + v[i-1,g_y] + v[i+1,g_y] + v[i,g_y-1]
    v_score[0,g_y] = 3*v[0,g_y] + v[1,g_y] + v[0,g_y-1]
    v_score[g_x,g_y] = 3*v[g_x,g_y] + v[g_x-1,g_y] + v[g_x,g_y-1]
    v_score[0,0] = 3*v[0,0] + v[1,0] + v[0,1]
    v_score[g_x,0] = 3*v[g_x,0] + v[g_x-1,0] + v[g_x,1]
    for i,j in g:
        mt[i,j] = v_score[i,j]*a + (2*random.random()-1)*a_std/2 + b
        e[i,j] = mt[i,j] + extended_plane[j]
    return e, mt

def generate_elevation(g, v, ep, a, a_std, b):
    v_score, mt, e, g_x, g_y = {}, {}, {}, list(g)[-1][0], list(g)[-1][1]
    for i in range(1,19):
        for j in range(1,59):
            v_score[i,j] = v[i,j] + v[i-1,j]+ v[i+1,j] + v[i,j-1] + v[i,j+1]
            v_score[0,j] = 2*v[0,j] + v[1,j] + v[0,j-1] + v[0,j+1]
            v_score[19,j] = 2*v[19,j] + v[18,j] + v[19,j-1] + v[19,j+1]
        v_score[i,0] = 2*v[i,0] + v[i-1,0] + v[i+1,0] + v[i,1]
        v_score[i,59] = 2*v[i,59] + v[i-1,59] + v[i+1,59] + v[i,58]
    v_score[0,59] = 3*v[0,59] + v[1,59] + v[0,58]
    v_score[19,59] = 3*v[19,59] + v[18,59] + v[19,58]
    v_score[0,0] = 3*v[0,0] + v[1,0] + v[0,1]
    v_score[19,0] = 3*v[19,0] + v[18,0] + v[19,1]
    for i,j in g:
        mt[i,j] = v_score[i,j]*a + (2*random.random()-1)*a_std + b
        e[i,j] = mt[i,j] + np.mean([ep[k,j] for k in range(20)])
    return e, mt

def generate_elevation_new(g, not_sinks_, v, ep, a, a_std, b):
    v_score, mt, e = {}, {}, {}
    for i in range(1,19):
        for j in range(1,59):
            v_score[i,j] = v[i,j] + v[i-1,j]+ v[i+1,j] + v[i,j-1] + v[i,j+1]
            v_score[0,j] = 2*v[0,j] + v[1,j] + v[0,j-1] + v[0,j+1]
            v_score[19,j] = 2*v[19,j] + v[18,j] + v[19,j-1] + v[19,j+1]
        v_score[i,0] = 2*v[i,0] + v[i-1,0] + v[i+1,0] + v[i,1]
        v_score[i,59] = 2*v[i,59] + v[i-1,59] + v[i+1,59] + v[i,58]
    v_score[0,59] = 3*v[0,59] + v[1,59] + v[0,58]
    v_score[19,59] = 3*v[19,59] + v[18,59] + v[19,58]
    v_score[0,0] = 3*v[0,0] + v[1,0] + v[0,1]
    v_score[19,0] = 3*v[19,0] + v[18,0] + v[19,1]
    for i,j in g:
        mt[i,j] = v_score[i,j]*a + (2*random.random()-1)*a_std + b
        e[i,j] = mt[i,j] + np.mean([ep[k,j] for k in range(20)])
    c = 1
    geff = eff_lattice(g, e)
    zero_nodes = [i for i in not_sinks_ if geff.out_degree(i)==0]
    while len(zero_nodes):
        #print(len(zero_nodes))
        for i,j in zero_nodes:
            next_nodes = list(g.successors((i,j)))
            closest_node = np.array(next_nodes)[np.argmin([e[k] for k in next_nodes])]
            #print(closest_node)
            mt[i,j] = mt[closest_node[0], closest_node[1]] + 1e-2
            #mt[i,j] = mt[i,j-1] + 1e-3
            e[i,j] = mt[i,j] + np.mean([ep[k,j] for k in range(20)])
        geff = eff_lattice(g, e)
        zero_nodes = [i for i in not_sinks_ if geff.out_degree(i)==0]
        c += 1
    #print(c)
    return e, mt

    for j in reversed(range(1,60)):
        for i in range(1,19):
            if e[i,j]<e[i-1,j] and e[i,j]<e[i+1,j] and e[i,j]<e[i,j-1]:
                e[i,j] = e[i,j-1]+1e-3
                mt[i,j] = e[i,j] - np.mean([ep[k,j] for k in range(20)])
        if e[0,j]<e[19,j] and e[0,j]<e[1,j] and e[0,j]<e[0,j-1]:
            e[0,j] = e[0,j-1]+1e-3
            mt[0,j] = e[0,j] - np.mean([ep[k,j] for k in range(20)])
        if e[19,j]<e[18,j] and e[19,j]<e[0,j] and e[19,j]<e[19,j-1]:
            e[19,j] = e[19,j-1]+1e-3
            mt[19,j] = e[19,j] - np.mean([ep[k,j] for k in range(20)])
    if e[0,59]<e[0,58] and e[0,59]<e[19,59] and e[0,59]<e[1,59]:
        e[0,59] = e[0,58]+1e-3
        mt[0,59] = e[0,59] - np.mean([ep[k,59] for k in range(20)])
    if e[19,59]<e[19,58] and e[19,59]<e[0,59] and e[19,59]<e[18,59]:
        e[19,59] = e[19,58]+1e-3
        mt[19,59] = e[19,59] - np.mean([ep[k,59] for k in range(20)])

## example 1: reconstructing the elevation profiles of the chosen landscape

V, E = VP[k], EP[k]
V = Vs[p][t]

E, MTmat = generate_elevation(G, VP[3], EP[3], alpha, alpha_std, beta)
plt.imshow(nx_to_mat(MTmat),cmap='rainbow')
plt.show()

#MTmat = eP[3]-np.array([[np.mean(eP[3][i]) for k in range(20)] for i in range(60)])
fig, ax = plt.  subplots(1,3,figsize=(6,6),tight_layout=True)
p0 = ax[0].imshow(vP[3], cmap='Greens')
p1 = ax[1].imshow(nx_to_mat(MT), cmap='rainbow')
p2 = ax[2].imshow(nx_to_mat(MTmat), cmap='rainbow')
p = p0, p1, p2
_ = [ax[i].set_xticks([]) for i in range(3)]
_ = [ax[i].set_yticks([]) for i in range(3)]
_ = [ax[i].set_title(' (a)'*int(i==0)+' (b)'*int(i==1)+' (c)'*int(i==2), loc='left', size=12) for i in range(3)]
_ = [fig.colorbar(p[i], ax=ax[i], orientation='horizontal', fraction=.02, pad=.02, shrink=.82) for i in range(3)]
plt.show()

fig, ax = plt.subplots(1,3,figsize=(6,6),tight_layout=True)
_ = [ax[i].imshow(nx_to_mat(Vs[i][t]), cmap='Greens') for i in range(3)]
_ = [ax[i].set_xticks([]) for i in range(3)]
_ = [ax[i].set_yticks([]) for i in range(3)]
_ = [ax[i].set_title(' (a)'*int(i==0)+' (b)'*int(i==1)+' (c)'*int(i==2), loc='left', size=12) for i in range(3)]
plt.show()

Geff = eff_lattice(G,E2)

def rerouting(g, not_sinks, ep):
    geff, k = eff_lattice(g,ep), 0
    wells = [i for i in geff if geff.out_degree(i)==0 and i in not_sinks]
    for i in wells:
        ep[list(g.successors(i))[np.argmin([ep[k] for k in g.successors(i)])]] = ep[i] - 1e-3
    geff = eff_lattice(g,ep)
    wells = [i for i in geff if geff.out_degree(i)==0 and i in not_sinks]
    while len(wells):
        k += 1
        for i in wells:
            ep[(i[0],i[1]-1)] = ep[i] - 1e-4
        geff = eff_lattice(g,ep)
        wells = [i for i in geff if geff.out_degree(i)==0 and i in not_sinks]
    print('required',k,'rerouting(s)')
    return ep

fig, ax = plt.subplots(1,3, figsize=(9,9), tight_layout=True)
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(V.values()), cmap='Greens', ax=ax[0])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(E.values()), cmap='rainbow', ax=ax[1])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MT.values()), cmap='rainbow', ax=ax[2])
plt.show()

## example 2: generating elevation for a randomly generated landscape
V = tclustering_square(G, vcover, pc)
E3, MT3 = generate_elevation(G, V, Ep, alpha, alpha_std, beta)

fig, ax = plt.subplots(1,2,figsize=(6,9))
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(V.values()), cmap='Greens', ax=ax[0])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MT3.values()), cmap='rainbow', ax=ax[1])
plt.show()

### SECTION 5: SAMPLING GRASSLAND

k = 1
VD, ED, Vp, Ep, MTp1 = np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'vegcover.asc', skiprows=6)[1+int(k==3):-1,1:-1], np.loadtxt('MAHLERAN_FC/p'+str(k+1)+'dem.asc', skiprows=6)[1:-1,1:-1], {}, {}, {}

MTarr = np.zeros((60,20), dtype=float)
for i,j in G:
    Vp[i,j], Ep[i,j] = VD[59-j,i], ED[59-j,i]
    MTp1[i,j] = Ep[i,j]-np.mean(ED[59-j,:])
    MTarr[59-j,i] = MTp1[i,j]

##
VD1, VD2, VD3, VD4 = np.loadtxt('MAHLERAN_FC/p1vegcover.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p2vegcover.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p3vegcover.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p4vegcover.asc', skiprows=6)[2:-1,1:-1]
ED1, ED2, ED3, ED4 = np.loadtxt('MAHLERAN_FC/p1dem.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p2dem.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p3dem.asc', skiprows=6)[1:-1,1:-1], np.loadtxt('MAHLERAN_FC/p4dem.asc', skiprows=6)[1:-1,1:-1]
MT1, MT2, MT3, MT4 = np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float), np.zeros((60,20), dtype=float)
for j in range(20):
    for i in range(60):
        MT1[i,j], MT2[i,j], MT3[i,j], MT4[i,j] = ED1[i,j]-np.mean(ED1[i]), ED2[i,j]-np.mean(ED2[i]), ED3[i,j]-np.mean(ED3[i]), ED4[i,j]-np.mean(ED4[i])

fig, ax = plt.subplots(1,8,figsize=(16,6),tight_layout=True)
ax[0].imshow(VD1, cmap='Greens')
ax[1].imshow(MT1, cmap='rainbow')
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[0].set_title('p1 (grassland)')
ax[2].imshow(VD2, cmap='Greens')
ax[3].imshow(MT2, cmap='rainbow')
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[3].set_xticks([])
ax[3].set_yticks([])
ax[2].set_title('p2 (grass-shrubland)')
ax[4].imshow(VD3, cmap='Greens')
ax[5].imshow(MT3, cmap='rainbow')
ax[4].set_xticks([])
ax[4].set_yticks([])
ax[5].set_xticks([])
ax[5].set_yticks([])
ax[4].set_title('p3 (shrub-grassland)')
ax[6].imshow(VD4, cmap='Greens')
ax[7].imshow(MT4, cmap='rainbow')
ax[6].set_xticks([])
ax[6].set_yticks([])
ax[7].set_xticks([])
ax[7].set_yticks([])
ax[6].set_title('p4 (shrubland)')
plt.show()


## Computing scores from chosen natural landscape, the score of a node i = v_i + sum(v_j) where j are the 4 neighbors of i
land_scores = np.zeros((60,20), dtype=float)
for j in range(1,19):
    for i in range(1,59):
    #    land_scores[i,j] = 1*VD[i,j] + VD[i-1,j] + VD[i+1,j] + VD[i,j-1]+ VD[i,j+1]
    #    land_scores[i,0] = 2*VD[i,0] + VD[i-1,0] + VD[i+1,0] + VD[i,1]
    #    land_scores[i,19] = 2*VD[i,19] + VD[i-1,19] + VD[i+1,19] + VD[i,18]
    #land_scores[0,j] = 2*VD[0,j] + VD[0,j-1] + VD[0,j+1] + VD[1,j]
    #land_scores[59,j] = 2*VD[59,j] + VD[59,j-1] + VD[59,j+1] + VD[58,j]
        
        land_scores[i,j] = 1*VD[i,j] - 0*(VD[i-1,j] + VD[i+1,j]) + VD[i,j-1]+ VD[i,j+1]
        land_scores[i,0] = 2*VD[i,0] - 0*(VD[i-1,0] + VD[i+1,0]) + VD[i,1]
        land_scores[i,19] = 2*VD[i,19] - 0*(VD[i-1,19] + VD[i+1,19]) + VD[i,18]
    land_scores[0,j] = VD[0,j] + VD[0,j-1] + VD[0,j+1] - 0*VD[1,j]
    land_scores[59,j] = VD[59,j] + VD[59,j-1] + VD[59,j+1] - 0*VD[58,j]

land_scores[0,0] = 2*VD[0,0] + VD[0,1] - 0*VD[1,0]
land_scores[59,0] = 2*VD[59,0] + VD[59,1] - 0*VD[58,0]
land_scores[0,19] = 2*VD[0,19] + VD[0,18] - 0*VD[1,19]
land_scores[59,19] = 2*VD[59,19] - 0*VD[58,19] + VD[59,18]
land_scores /= 100

MTarr = np.zeros((60,20), dtype=float)
for j in range(20):
    for i in range(60):
        MTarr[i,j] = MTp1[j,59-i]

pearsonr(MTarr.flatten(), land_scores.flatten())
plt.scatter(MTarr.flatten(), land_scores.flatten(), marker='x', s=50)
plt.show()

fig, ax = plt.subplots(1,2,figsize=(6,9))
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(Vp.values()), cmap='Greens', ax=ax[0])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MTp1.values()), cmap='rainbow', ax=ax[1])
plt.show()

### SECTION 4: GENERATING SHRUBLAND

## Linear regression between measured scores vs microtopography of each node
X, Y = [land_scores[59-j,i] for i,j in G], [Ep[i,j]-np.mean([Ep[l,j] for l in range(20)]) for i,j in G]
alpha, beta = np.polyfit(X,Y,1)
Y_line = [i*alpha+beta for i in X]
alpha_std = np.std([Y[k]-Y_line[k] for k in range(len(G))])

## example 1: reconstructing the elevation profiles of the chosen landscape
E2, MT2 = generate_elevation(G, Vp, Ep, alpha, alpha_std, beta)

fig, ax = plt.subplots(1,3, figsize=(9,9), tight_layout=True)
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(Vp.values()), cmap='Greens', ax=ax[0])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MTp4.values()), cmap='rainbow', ax=ax[1])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MT2.values()), cmap='rainbow', ax=ax[2])

## example 2: generating elevation for a randomly generated landscape
V = tclustering_square(G, vcover, pc)
E3, MT3 = generate_elevation(G, V, Ep, alpha, alpha_std, beta)

fig, ax = plt.subplots(1,2,figsize=(6,9))
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(V.values()), cmap='Greens', ax=ax[0])
nx.draw(G, dict(zip(G,G)), edgelist=[], node_size=80, node_color=list(MT3.values()), cmap='rainbow', ax=ax[1])
plt.show()

#### SECTION 5: WIDELY SAMPLED SCORE (2nd neighborhood) → calculations currently messed up
land_scores = np.zeros((60,20), float)
for j in range(1,19):
    for i in range(1,59):
        land_scores[i,j] = 1*VD[i,j] + VD[i-1,j] + VD[i+1,j] + VD[i,j-1]+ VD[i,j+1]
        land_scores[i,j] += .5*(VD[i-1,j-1] + VD[i+1,i-1] + VD[i-1,i+1] + VD[i+1,i+1])
        land_scores[i,0] = 2*VD[i,0] + VD[i-1,0] + VD[i+1,0] + VD[i,1]
        land_scores[i,0] += 0.5*(VD[i-1,0] + VD[i+1,0] + VD[i-1,1] + VD[i+1,1])
        land_scores[i,19] = 2*VD[i,19] + VD[i-1,19] + VD[i+1,19] + VD[i,18]
        land_scores[i,19] += 0.5*(VD[i-1,18] + VD[i-1,18] + VD[i-1,19] + VD[i+1,19])
    land_scores[0,j] = 2*VD[0,j] + VD[0,j-1] + VD[0,j+1] + VD[1,j]
    land_scores[0,j] += 0.5*(VD[0,j-1] + VD[0,j+1] + VD[1,j-1] + VD[1,j+1])
    land_scores[59,j] = 2*VD[59,j] + VD[59,j-1] + VD[59,j+1] + VD[58,j]
    land_scores[59,j] += 0.5*(VD[59,j-1] + VD[59,j+1] + VD[58,j-1] + VD[58,j+1])

for j in range(2,18):
    for i in range(2,58):
        land_scores[i,j] += 0.5*(VD[i-1,j] + VD[i+1,j] + VD[i,j-1] + VD[i,j+1])

land_scores[0,0] = 3*VD[0,0] + VD[0,1] + VD[1,0]
land_scores[0,0] += 0.5*(3*VD[0,0] + VD[0,1] + VD[1,0]+ VD[0,2] + VD[2,0])
land_scores[59,0] = 3*VD[59,0] + VD[59,1] + VD[58,0]
land_scores[59,0] += 0.5*(3*VD[59,0] + VD[59,1] + VD[58,0] + VD[59,2] + VD[57,0])
land_scores[0,19] = 3*VD[0,19] + VD[0,18] + VD[1,19]
land_scores[0,19] += 0.5*(3*VD[0,19] + VD[0,18] + VD[1,19] + VD[0,17] + VD[2,19])
land_scores[59,19] = 3*VD[59,19] + VD[58,19] + VD[59,18]
land_scores[59,19] += 0.5*(3*VD[59,19] + VD[58,19] + VD[59,18] + VD[59,17] + VD[57,19])

scores2_2 = land_scores/100

plt.scatter(land1.flatten(), nx_to_mat(MT1).flatten(), marker='x')
plt.scatter(land4.flatten(), nx_to_mat(MT4).flatten(), marker='x')
plt.plot([0,5], [beta1, 5*alpha1+beta1], c='r')
plt.plot([0,5], [beta4, 5*alpha4+beta4], c='r')
plt.show()


X1, Y1 = [land1[59-j,i] for i,j in G], [Ep1[i,j]-np.mean([Ep1[l,j] for l in range(20)]) for i,j in G]
alpha1, beta1 = np.polyfit(X1,Y1,1)
Y1_line = [i*alpha1+beta1 for i in X1]
alpha1_std = np.std([Y1[k]-Y1_line[k] for k in range(len(G))])

X4, Y4 = [land4[59-j,i] for i,j in G], [Ep4[i,j]-np.mean([Ep4[l,j] for l in range(20)]) for i,j in G]
alpha4, beta4 = np.polyfit(X4,Y4,1)
Y4_line = [i*alpha1+beta1 for i in X4]
alpha4_std = np.std([Y4[k]-Y4_line[k] for k in range(len(G))])

fig = plt.figure(figsize=(10.3,6), tight_layout=True)
gs = gridspec.GridSpec(nrows=2, ncols=4, height_ratios=[1,1], width_ratios=[2.2,1,1,1])
ax_scatter = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[1,0])
ax_maps = fig.add_subplot(gs[:,1]), fig.add_subplot(gs[:,2]), fig.add_subplot(gs[:,3])
for ax in ax_maps:
    _ = ax.set_xticks([])
    _ = ax.set_yticks([])

ax_scatter[0].scatter(land1.flatten(), nx_to_mat(MT1).flatten(), marker='x', s=16, label='grassland')
ax_scatter[0].plot([0,5], [beta1, 5*alpha1+beta1], c='r')
ax_scatter[0].text(x=4.37, y=.01, s=str(round(1e3*alpha1,1))+' 10$^{-3}$', c='r', weight='bold', size=9)
ax_scatter[0].set_xlabel('surrounding vegetation density')
ax_scatter[0].set_ylabel('microtopography')
ax_scatter[0].legend()
ax_scatter[1].scatter(land4.flatten(), nx_to_mat(MT4).flatten(), marker='x', s=16, label='shrubland')
ax_scatter[1].plot([0,5], [beta4, 5*alpha4+beta4], c='r')
ax_scatter[1].text(x=4.16, y=.015, s=str(round(1e3*alpha4,1))+' 10$^{-3}$', c='r', weight='bold', size=9)
ax_scatter[1].set_xlabel('surrounding vegetation density')
ax_scatter[1].set_ylabel('microtopography')
ax_scatter[1].legend()
ax_scatter[0].set_title(' (a)', loc='left')
ax_scatter[1].set_title(' (b)', loc='left')
ax_maps[0].set_title(' (c)', loc='left')
ax_maps[1].set_title(' (d)', loc='left')
ax_maps[2].set_title(' (e)', loc='left')
p = ax_maps[0].imshow(vP[3]/100, cmap='Greens'), ax_maps[1].imshow(nx_to_mat(MT4), cmap='rainbow'), ax_maps[2].imshow(nx_to_mat(MT2), cmap='rainbow')
divider = make_axes_locatable(ax_maps[0]), make_axes_locatable(ax_maps[1]), make_axes_locatable(ax_maps[2])
cax = [d.append_axes('bottom', size='2%', pad=0.1) for d in divider]
_ = [fig.colorbar(p[i], cax=cax[i], orientation='horizontal') for i in range(3)]
plt.show()

from mpl_toolkits.axes_grid1 import make_axes_locatable







