import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from operator import mul
from functools import reduce
from scipy.stats import pearsonr
import pickle, random, itertools, sys
from math import log
from collections import Counter

def vgrid(x, y):
    g = nx.DiGraph()
    for j in range(y+1):
        for i in range(x):
            g.add_edge( (2*i-1, 2*(y-j)+1), (2*i,2*(y-j)) )
            g.add_edge( (2*i+1, 2*(y-j)+1), (2*i,2*(y-j)) )
            g.add_edge( (2*i, 2*(y-j)), (2*i-1, 2*(y-j)-1) )
            g.add_edge( (2*i, 2*(y-j)), (2*i+1, 2*(y-j)-1) )
        g.add_edge( (1, 2*(y-j)+1), (0, 2*(y-j)) )
    for j in range(1, y+1):
        g.add_edge( (2*x-1, 2*(y-j)+1), (0, 2*(y-j)) )
        g.add_edge( (0, 2*(y-j)), (2*x-1, 2*(y-j)-1) )
    for i,j in list(g):
        if i<0 or j<0 or i>=2*x or j>=2*y:
            g.remove_node((i,j))
    s_ = [(2*i, 0) for i in range(x)]
    return g, [i for i in g if i not in s_]

def tgrid(x, y):
    g = nx.grid_2d_graph(x,int(2*y))
    g.remove_edges_from([((i1,0),(i1+1,0)) for i1 in range(x-1)])
    g = g.to_directed()
    g.remove_edges_from([((i1,i2),(i1,i2+1)) for i1,i2 in itertools.product(range(x), range(int(2*y)-1))])
    g.add_edges_from([((0,i2),(x-1,i2)) for i2 in range(1,int(2*y))]+[((x-1,i2),(0,i2)) for i2 in range(1,int(2*y))])
    s_ = [(i, 0) for i in range(x)]
    return g, [i for i in g if i not in s_]

def hgrid(x, y):
    g = nx.grid_2d_graph(x,int(2*y))
    g = g.to_directed()
    g.remove_edges_from(list(g.edges))
    g.add_edges_from([((i,j+1),(i+1,j)) for i,j in itertools.product(range(x-1),range(2*y-1))])
    g.add_edges_from([((i+1,j+1),(i,j)) for i,j in itertools.product(range(x-1),range(2*y-1))])
    g.add_edges_from([((0,j+1),(x-1,j)) for j in range(2*y-1)] + [((x-1,j+1),(0,j)) for j in range(2*y-1)]) # pbc
    g.add_edges_from([((i,j),(i,j-1)) for i,j in itertools.product(range(x),range(1,2*y))])
    s_ = [(i, 0) for i in range(x)]
    return g, [i for i in g if i not in s_]

def init(g):
    state, coupling, current_av, branches, new_active, size, ordlist = {}, {}, [], {}, [], [], sorted(list(g), key=lambda a:a[1], reverse=True)
    for i in g:
        state[i], coupling[i] = 0, {}
    for i,j in itertools.combinations(ordlist, 2):
        if nx.has_path(g,i,j):
            coupling[i][j] = []
        if i[1] == j[1]:		# for tgrid
            if nx.has_path(g,j,i):
                coupling[j][i] = []
    return state, coupling, current_av, branches, new_active, size

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

def SC_path_thres(g, inter_dist, lim):
    sc, gconn = {}, deepcopy(g)
    gconn.remove_edges_from(list(g.edges))
    for i in inter_dist:
        for j in inter_dist[i]:
            if inter_dist[i][j] <= lim:
                _ = gconn.add_edge(i,j,weight=inter_dist[i][j])
    for i in g:
        sc[i] = gconn.in_degree(i, weight='weight') - gconn.out_degree(i, weight='weight')
    return sc

def slope(l):
    xf = max(l)//5
    yf = log(sum((i//5 in [xf,xf+1]) for i in l))
    xi = 0
    yi = log(sum((i in [0,5]) for i in l))
    return (yf-yi)/(xf-xi)

def slope_line(x,y):
    return (y[-1]-y[0])/(x[-1]-x[0])

def sandpile(g, v, ep, state, coupling, current_av, branches, old_active, size, not_sinks):# exit_,
    new_active, unstable = [], [node for node in old_active if state[node] >= g.out_degree(node)]
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                spreading_scheme = Counter(random.choices(list(g.successors(node1)), weights=ep[(node1)], k=nb_partcl))
                for node2 in spreading_scheme:
                    state[node2] += np.sum(np.random.rand(spreading_scheme[node2]) > v[node2]/2)
                    current_av[node2] = current_av[node1]+1
                    if node2 not in branches:
                        branches[node2] = [node1]
                    elif node1 not in branches[node2]:
                        branches[node2].append(node1)
                    # record all the previous sources of the new incident node (at least, from a new path)
                    origin, all_sources, sources, counted = list(current_av.keys())[0], [node1], [node1], []
                    while origin not in sources:
                        new_sources = []
                        for source in sources:
                            new_sources += branches[source]
                        new_sources = list(set(new_sources))
                        all_sources += new_sources
                        sources = new_sources
                    # record couplings between new incident node and all its previous sources
                    for source in all_sources: # collapse 'all_sources' to avoid redundance of nodes
                        coupling[source][node2].append(current_av[node2]-current_av[source])
                    if node2 not in new_active:
                        new_active.append(node2)
            #else:
            #    exit_ += nb_partcl
    else:
        if len(current_av):
            size.append(len(current_av)-1)
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        #print(node2, v[node2], type(node2), type(v[node2]))
        #if v[node2] == 0:
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, coupling, current_av, branches, new_active, size#, exit_

def sandpile_only_scfc(g, v, ep, state, coupling, current_av, branches, old_active, not_sinks):
    new_active, unstable = [], [node for node in old_active if state[node] >= g.out_degree(node)]
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                spreading_scheme = Counter(random.choices(list(g.successors(node1)), weights=ep[(node1)], k=nb_partcl))
                for node2 in spreading_scheme:
                    state[node2] += np.sum(np.random.rand(spreading_scheme[node2]) > v[node2]/2)
                    current_av[node2] = current_av[node1]+1
                    if node2 not in branches:
                        branches[node2] = [node1]
                    elif node1 not in branches[node2]:
                        branches[node2].append(node1)
                    # record all the previous sources of the new incident node (at least, from a new path)
                    origin, all_sources, sources, counted = list(current_av.keys())[0], [node1], [node1], []
                    while origin not in sources:
                        new_sources = []
                        for source in sources:
                            new_sources += branches[source]
                        new_sources = list(set(new_sources))
                        all_sources += new_sources
                        sources = new_sources
                    # record couplings between new incident node and all its previous sources
                    for source in all_sources: # collapse 'all_sources' to avoid redundance of nodes
                        coupling[source][node2].append(current_av[node2]-current_av[source])
                    if node2 not in new_active:
                        new_active.append(node2)
    else:
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, coupling, current_av, branches, new_active

def sandpile_only_avalanches(g, v, ep, state, current_av, branches, old_active, size, dura, not_sinks):# exit_,
    new_active, unstable = [], [node for node in old_active if state[node] >= g.out_degree(node)]
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                spreading_scheme = Counter(random.choices(list(g.successors(node1)), weights=ep[(node1)], k=nb_partcl))
                for node2 in spreading_scheme:
                    state[node2] += np.sum(np.random.rand(spreading_scheme[node2]) > v[node2]/2)
                    current_av[node2] = current_av[node1]+1
                    if node2 not in branches:
                        branches[node2] = [node1]
                    elif node1 not in branches[node2]:
                        branches[node2].append(node1)
                    # record all the previous sources of the new incident node (at least, from a new path)
                    origin, all_sources, sources, counted = list(current_av.keys())[0], [node1], [node1], []
                    while origin not in sources:
                        new_sources = []
                        for source in sources:
                            new_sources += branches[source]
                        new_sources = list(set(new_sources))
                        all_sources += new_sources
                        sources = new_sources
                    # record couplings between new incident node and all its previous sources
                    if node2 not in new_active:
                        new_active.append(node2)
    else:
        if len(current_av):
            size.append(len(current_av)-1)
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, current_av, branches, new_active, size, dura

def sandpilest(g, v, state, coupling, current_av, branches, old_active, size, exit_, not_sinks):#duration, 
    unstable, new_active = [], []
    for node in old_active:
        if state[node] > g.out_degree(node):
            unstable.append(node)
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                node2 = list(g.successors(node1))[0]
                #if node2 in not_sinks:
                #    state[node2] += np.sum(np.random.rand(nb_partcl) > v[node2]/2)
                state[node2] += np.sum(np.random.rand(nb_partcl) > v[node2]/2)
                current_av[node2] = current_av[node1]+1
                if node2 not in branches:
                    branches[node2] = [node1]
                elif node1 not in branches[node2]:
                    branches[node2].append(node1)
                # record all the previous sources of the new incident node (at least, from a new path)
                origin, all_sources, sources, counted = list(current_av.keys())[0], [node1], [node1], []
                while origin not in sources:
                    new_sources = []
                    for source in sources:
                        new_sources += branches[source]
                    new_sources = list(set(new_sources))
                    all_sources += new_sources
                    sources = new_sources
                # record couplings between new incident node and all its previous sources
                for source in all_sources: # collapse 'all_sources' to avoid redundance of nodes
                    coupling[source][node2].append(current_av[node2]-current_av[source])
                if node2 not in new_active:
                    new_active.append(node2)
            else:
                exit_ += nb_partcl
    else:
        if len(current_av):
            size.append(len(current_av))
            #duration.append(max(list(current_av.values())))
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, coupling, current_av, branches, new_active, size, exit_

def simple_sandpilest_sedigraph_with_rain_intensity(g, v, state, current_av, branches, old_active, nb_rain_, exit_, p_rain_, not_sinks):#duration, 
    new_active, unstable = [], [node for node in old_active if state[node] >= g.out_degree(node)]
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                node2 = list(g.successors(node1))[0]
                state[node2] += np.sum(np.random.rand(nb_partcl) > v[node2]/2)
                current_av[node2] = current_av[node1]+1
                if node2 not in branches:
                    branches[node2] = [node1]
                elif node1 not in branches[node2]:
                    branches[node2].append(node1)
                if node2 not in new_active:
                    new_active.append(node2)
            else:
                exit_ += nb_partcl
    elif random.random() < p_rain_:
        current_av, branches, nb_rain_ = {}, {}, nb_rain_+1
        node2 = random.choice(not_sinks)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, current_av, branches, new_active, nb_rain_, exit_

def sandpilest_vary_dissipation(g, v, state, coupling, current_av, branches, old_active, size, exit_, not_sinks, d_coeff):#duration, 
    unstable, new_active = [], []
    for node in old_active:
        if state[node] > g.out_degree(node):
            unstable.append(node)
    if len(unstable):
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if g.out_degree(node1):
                node2 = list(g.successors(node1))[0]
                #if node2 in not_sinks:
                #    state[node2] += np.sum(np.random.rand(nb_partcl) > v[node2]/2)
                state[node2] += np.sum(np.random.rand(nb_partcl) > d_coeff*v[node2])
                current_av[node2] = current_av[node1]+1
                if node2 not in branches:
                    branches[node2] = [node1]
                elif node1 not in branches[node2]:
                    branches[node2].append(node1)
                # record all the previous sources of the new incident node (at least, from a new path)
                origin, all_sources, sources, counted = list(current_av.keys())[0], [node1], [node1], []
                while origin not in sources:
                    new_sources = []
                    for source in sources:
                        new_sources += branches[source]
                    new_sources = list(set(new_sources))
                    all_sources += new_sources
                    sources = new_sources
                # record couplings between new incident node and all its previous sources
                for source in all_sources: # collapse 'all_sources' to avoid redundance of nodes
                    coupling[source][node2].append(current_av[node2]-current_av[source])
                if node2 not in new_active:
                    new_active.append(node2)
            else:
                exit_ += nb_partcl
    else:
        if len(current_av):
            size.append(len(current_av))
            #duration.append(max(list(current_av.values())))
        current_av, branches = {}, {}
        node2 = random.choice(not_sinks)
        if random.random() > d_coeff*v[node2]:
            state[node2] += 1
            new_active.append(node2)
    return state, coupling, current_av, branches, new_active, size, exit_

def SC_area(g, c):
    sc, gconn = {}, deepcopy(g)
    gconn.remove_edges_from(list(g.edges))
    for i in c:
        for j in c[i]:
            _ = gconn.add_edge(i,j)
    for i in g:
        sc[i] = gconn.in_degree(i)
    return sc

def SCw_area(g, c, epmap):
    sc, i, gs = {}, 0, {}
    for v1 in c:
        gs[v1] = nx.DiGraph()
    for v1 in c:
        print('NODE N°',i,':',v1,'HAS',len(c[v1]),' possible targets...')
        i+=1
        for v2 in c[v1]:
            nb_steps = nx.shortest_path_length(g,v1,v2)
            if nb_steps<=30:
                paths = list(nx.all_shortest_paths(g,v1,v2))
                for p in paths:
                    gs[v2].add_weighted_edges_from([(p[k],p[k+1],epmap[(p[k],p[k+1])]) for k in range(len(p)-1)])
    for v1 in c:
        print('For',v1,', edgelist=',list(gs[v1].edges))
        sc[v1] = gs[v1].size(weight='weight')
    return gs, sc

def SC_EDGE_before(g,epmap):
    gs = {}
    for i in g:
        gs[i] = nx.DiGraph()
    for j in range(60):
        for i in range(20):
            v1 = (i,59-j)
            source_edges = list((v2,v1) for v2 in g.predecessors(v1))
            while source_edges:
                next_source_edges = source_edges
                for s,t in source_edges:
                    gs[v1].add_edge(s,t,weight=epmap[(s,t)])
                    next_source_edges.remove((s,t))
                    if g.predecessors(s):
                        for node in g.predecessors(s):
                            if node[1]==59-j:
                                next_source_edges += [(node,s) for node in list(g.predecessors(s))]
                            else:
                                gs[v1].add_weighted_edges_from([(i,j,gs[node][i][j]['weight']) for i,j in gs[node].edges])
                source_edges = next_source_edges
    return gs
    return {node:sum(gs[node].edges(data=True)) for node in gs}, {node:len(gs[node].edges) for node in gs}, {node:len(gs[node]) for node in gs}

# now use this one
def SC_EDGE(g,epmap):
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

def SC_PATH(g,c,epmap):
    gconn = deepcopy(g)
    gconn.remove_edges_from(list(g.edges()))
    for j in range(60):
        next_layer_successors, same_layer_successors = {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]-1] for v1 in [(i,j) for i in range(20)]}, {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]] for v1 in [(i,j) for i in range(20)]}
        for i in range(20):
            v1 = (i,j)
            gconn.add_weighted_edges_from([(v1,v2,epmap[(v1,v2)]) for v2 in next_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,epmap[(v1,v2)]) for v3 in gconn.successors(v2)] for v2 in next_layer_successors[v1]],[]))
            for v2 in next_layer_successors[v1]:
                for v3 in gconn.successors(v2):
                    gconn[v1][v3]['weight'] += gconn[v2][v3]['weight']*epmap[(v1,v2)]
        for v1 in sort_layer(same_layer_successors):
            gconn.add_weighted_edges_from([(v1,v2,epmap[(v1,v2)]) for v2 in same_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,gconn[v2][v3]['weight']) for v3 in gconn.successors(v2)] for v2 in same_layer_successors[v1]],[]))
            for v2 in same_layer_successors[v1]:
                for v3 in gconn.successors(v2):
                    gconn[v1][v3]['weight'] += gconn[v2][v3]['weight']*epmap[(v1,v2)]
    return {node:gconn.in_degree(node, weight='weight') for node in g}

def sort_layer(nodeset):
    gl = nx.DiGraph()
    gl.add_nodes_from(nodeset.keys())
    gl.add_edges_from([(i,nodeset[i][0]) for i in nodeset if len(nodeset[i])])
    return list(nx.topological_sort(gl))[::-1]

def SC_PATH_before(g,c,epmap):
    sc, i, gconn, nlist, counted = {}, 0, deepcopy(g), sorted(list(g), key=lambda x:x[1]), []
    gconn.remove_edges_from(list(g.edges))
    for v1 in nlist:
        print('NODE N°',i,':',v1,'HAS',len(c[v1]),' possible targets...')
        i+=1
        for v2 in c[v1]:
            if v2 not in counted:
                if nx.shortest_path_length(g,v1,v2)<=30:
                    paths = list(nx.all_shortest_paths(g,v1,v2))
                    prod = np.sum([np.prod([epmap[(p[i],p[i+1])] for i in range(len(p)-1)]) for p in paths])
                    gconn.add_edge(v1,v2,weight=prod)
            else:
                gconn.add_edge(v1,v2,weight=gconn.in_degree(v2,weight='weight')+epmap[(v1,v2)])
    for v1 in c:
        sc[v1] = gconn.in_degree(v1,weight='weight') 
    return sc

def add_veg_to_area(gs, sc, v):
    sc_new = {}
    for v1 in gs:
        g = gs[v1]
        sc_new[v1] = sc[v1]*sum(v[v2] for v2 in g)
    return sc_new

# keep this one
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

# fix increments on paths
def SCw_path_quick20(g,c,epmap):
    sc, i, gconn, nlist, counted = {}, 0, deepcopy(g), sorted(list(g), key=lambda x:x[1]), []
    gconn.remove_edges_from(list(g.edges))
    #gconn_uw.remove_edges_from(list(g.edges))
    for v1 in nlist:
        print('NODE N°',i,':',v1,'HAS',len(c[v1]),' possible targets...')
        i+=1
        for v2 in c[v1]:
            if v2 not in counted:
                if nx.shortest_path_length(g,v1,v2)<=30:
                    paths = list(nx.all_shortest_paths(g,v1,v2))
                    prod = np.sum([np.prod([epmap[(p[i],p[i+1])] for i in range(len(p)-1)]) for p in paths])
                    gconn.add_edge(v1,v2, weight=prod)
                    #gconn_uw.add_edge(v1,v2, weight=len(paths))
            else:
                gconn.add_edge(v1,v2, weight=gconn.in_degree(v2,weight='weight')+epmap[(v1,v2)])
                #gconn_uw.add_edge(v1,v2, weight=gconn_uw.in_degree(v2,weight='weight')+1)
    for v1 in c:
        sc[v1] = gconn.in_degree(v1,weight='weight')
    return sc, sc_uw

def SCw_path_quick3(g,c,epmap):
    i, nlist, counted = 0, sorted(list(g), key=lambda x:x[1]), {v1:[] for v1 in g}
    for v1 in nlist:
        print('NODE N°',i,':',v1,'HAS',len(c[v1]),' possible targets...')
        i+=1
        for v2 in c[v1]:
            if counted[v2]:
                v1_to_v2 = list(nx.all_shortest_paths(g,v1,v2))
                print('from',v1,'to',v2,', there are',len(v1_to_v2),'paths')
                p = v1_to_v2
                path_prod = np.prod([epmap[(p[i],p[i+1])] for i in range(len(p)-1)])
                counted[v1] += [path_prod*i for i in counted[v2]]
            elif nx.shortest_path_length(g,v1,v2)<=30:
                paths = list(nx.all_shortest_paths(g,v1,v2))
                counted[v1] += [np.prod([epmap[(p[i],p[i+1])] for i in range(len(p)-1)]) for p in paths]
    return {v1:sum(counted[v1]) for v1 in g}, {v1:len(counted[v1]) for v1 in g}

def FC_path(c, tt):
    fc, fci = {}, {}
    for i in c:
        fci[i] = []
    for i in c:
        for j in c[i]:
            fci[j] += c[i][j]
    if tt == 'len':
        for i in c:
            fc[i] = 0
            if len(fci[i]):
                fc[i] += len(fci[i])
    elif tt == 'sum':
        for i in c:
            fc[i] = 0
            if len(fci[i]):
                fc[i] += sum(fci[i])
    elif tt == 'var':
        for i in c:
            fc[i] = 0
            if len(fci[i]):
                fc[i] += np.var(fci[i])
    return fc

def FC_transport(c):
    fc = {i:0 for i in c}
    for i in c:
        for j in c[i]:
            if c[i][j]:
                fc[j] += np.var(c[i][j])
    return fc

def FC_deposit(c):
    fc, fci, fco = {i:0 for i in c}, {i:0 for i in c}, {i:0 for i in c}
    for i in c:
        for j in c[i]:
            if len(c[i][j]):
                fci[j] += len(c[i][j])
                #fco[i] += len(c[i][j])
    for i in c:
        fc[i] += fci[i]
    #    fc[i] -= fco[i]
    return fc

def FC_nb(c):
    fc, fci, fco = {i:0 for i in c}, {i:0 for i in c}, {i:0 for i in c}
    for i in c:
        for j in c[i]:
            if len(c[i][j]):
                fci[j] += len(c[i][j])
                fco[i] += len(c[i][j])
    for i in c:
        fc[i] += fci[i]
        fc[i] -= fco[i]
    return fc, fci, fco

def nx_to_mat(x):
    mat = np.zeros((60,20), float)
    for i in range(60):
        for j in range(20):
            mat[i,j] = x[(j,59-i)]
    return mat

def nx_to_mat2(x):
    g_x, g_y = list(x.keys())[-1]
    mat = np.zeros((g_y+1,g_x+1), float)
    for i in range(g_y+1):
        for j in range(g_x+1):
            mat[i,j] = x[(j,g_y-i)]
    return mat

def zeros(coupling):
    state, coupling_empty, current_av, branches, new_active, size = {}, {}, [], {}, [], []
    for i in coupling:
        state[i], coupling_empty[i] = 0, {}
        for j in coupling[i]:
            coupling_empty[i][j] = []
    return state, coupling_empty, current_av, branches, new_active, size

def zeros_no_avalanches(c):
    state, c_empty, current_av, branches, new_active = {}, {}, [], {}, []
    for i in c:
        state[i], c_empty[i] = 0, {}
        for j in c[i]:
            c_empty[i][j] = []
    return state, c_empty, current_av, branches, new_active

def eff_lattice(g, e):
    geff = deepcopy(g)
    for i in g:
        for j in g.successors(i):
            if e[i] <= e[j]:
                geff.remove_edge(i,j)
    for i in g:
        if geff.out_degree(i)==0:
            geff.add_edge(i, (i[0],i[1]-1))
    return geff

def eff_lattice(g, e): # this one (07.02.2025)
    geff = nx.DiGraph()
    geff.add_nodes_from(g.nodes)
    for i in g:
        for j in g.successors(i):
            #print(i,j,e[i],e[j])
            if e[i] >= e[j]:
                geff.add_edge(i,j)
    return geff

def eff_lattice(g, e): # last one = robust one (written on 02.06.2025)
    geff = deepcopy(g)
    for i in g:
        for j in g.successors(i):
            if e[i] <= e[j]:
                geff.remove_edge(i,j)
    for i,j in geff.edges:
        if (j,i) in geff.edges:
            if e[j] < e[i]:
                geff.remove_edge(j,i)
                if len(geff.successors(j)) == 0:
                    geff.add_edge(j,(j[0],j[1]-1))
            else:
                geff.remove_edge(i,j)
                if len(geff.successors(i)) == 0:
                    geff.add_edge(i,(i[0],i[1]-1))
    return geff

def eff_elev(g, e):
    geff = deepcopy(g)
    for i in g:
        for j in g.successors(i):
            if e[i] < e[j]:
                geff.remove_edge(i,j)
    return geff

def eff_steep(g, est):
    gst = deepcopy(g)
    gst.remove_edges_from(list(g.edges))
    for i in est:
        _ = gst.add_edge(i, est[i])
    return gst

def eprob(g, e):# eprob == propagation_probability
    ep = {}
    for i in g:
        ep[i] = []
        for j in g.successors(i):
            if e[i]>e[j]:
                ep[i].append(e[i]-e[j])
        if len(list(g.successors(i)))==1:
            ep[i] = [1]
    return ep

def eprob_map(g, e):# eprob_map == edge_slope
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

def estp(g, e):
    est = {}
    for i in g:
        next_nodes = list(g.successors(i))
        if len(next_nodes):
            est[i] = next_nodes[np.argmin(next_elev)]
    return est

def estp(g, e):
    est = {}
    for i in g:
        next_nodes = list(g.successors(i))
        if len(next_nodes):
            elev_nodes = np.array([e[k] for k in next_nodes])
            min_elev = np.where(elev_nodes==np.min(elev_nodes))[0]
            if len(min_elev)>1:
                est[i] = (i[0],i[1]-1)
            else:
                est[i] = next_nodes[min_elev[0]]
            #est[i] = next_nodes[np.argmin(next_elev)]
    return est

def SCFC(sc, fc, not_sinks):
    return pearsonr(list(sc[i] for i in not_sinks), list(fc[i] for i in not_sinks))[0]

def SCFC2(sc, fc, not_sinks):
    return spearmanr(list(sc[i] for i in not_sinks), list(fc[i] for i in not_sinks))[0]

def veg_vclust(g, cover, not_sinks_, clust):
    v, nb_added_veg = {}, 0
    for node in g:
        v[node] = 0
    empty_nodes, veg_nodes = deepcopy(not_sinks_), []
    target_nb_veg = int(len(not_sinks_)*cover)
    node = random.choice(empty_nodes)
    v[node] = random.random()
    _ = empty_nodes.remove(node)
    while len(veg_nodes) < target_nb_veg:
        n1, n2 = node
        if random.random() < clust:
            all_neighbs = list(g.predecessors(node)) + list(g.successors(node)) + [int((n1-1,n2) in g)*(n1-1,n2), int((n1+1,n2) in g)*(n1+1,n2)]
            free_neighbs = [j for j in all_neighbs if j in empty_nodes]
            while len(free_neighbs) == 0:
                node = random.choice(veg_nodes)
                n1, n2 = node
                all_neighbs = list(g.predecessors(node)) + list(g.successors(node)) + [int((n1-1,n2) in g)*(n1-1,n2), int((n1+1,n2) in g)*(n1+1,n2)]
                free_neighbs = [j for j in all_neighbs if j in empty_nodes]
            node = random.choice(free_neighbs)
        else:
            all_neighbs = list(g.predecessors(node)) + list(g.successors(node)) + [int((n1-1,n2) in g)*(n1-1,n2), int((n1+1,n2) in g)*(n1+1,n2)]
            node = random.choice([j for j in empty_nodes if j not in all_neighbs])
        v[node] = random.random()
        veg_nodes.append(node)
        _ = empty_nodes.remove(node)
    return v

def veg_tclust(g, x1_, x2_, cover, not_sinks_, clust):
    v = {}
    for node in g:
        v[node] = 0
    empty_nodes = deepcopy(not_sinks_)
    nb_veg, nb_added_veg = int(g.order()*cover), 1
    node = random.choice(empty_nodes)
    v[node] = random.random()
    _ = empty_nodes.remove(node)
    while nb_added_veg < nb_veg:
        n1, n2 = node
        if n1 == 0:
            all_neighbs = [(x1_-1,n2), (n1+1,n2)]
        elif n1 == x1_-1:
            all_neighbs = [(n1-1,n2), (0,n2)]
        else:
            all_neighbs = [(n1-1,n2), (n1+1,n2)]
        if n2 == int(2*x2_)-1:
            all_neighbs.append((n1,n2-1))
        elif n2 > 0:
            all_neighbs.append((n1,n2-1))
            all_neighbs.append((n1,n2+1))
        free_neighbs = [j for j in all_neighbs if j in empty_nodes]
        if len(free_neighbs) == 0:
            node = random.choice(all_neighbs)
        else:
            if random.random() < clust:
                node = random.choice(free_neighbs)
            else:
                node = random.choice([j for j in empty_nodes if j not in all_neighbs])
            v[node] = random.random()
            nb_added_veg += 1
            _ = empty_nodes.remove(node)
    return v

def tclustering(g, cover, pclust):
    v = {}
    for i in g:
        v[i] = 0
    free_nodes, nb_veg_final, nb_veg = list(g), int(len(g)*cover), 0
    node = random.choice(free_nodes)
    while nb_veg < nb_veg_final:
        v[node] = random.random()
        free_nodes.remove(node)
        nb_veg += 1
        if random.random() < pclust:
            neighbors = set(G.predecessors(node))|set(G.successors(node))
            free_neighbors = neighbors & set(free_nodes)
            if free_neighbors:
                node = random.choice(list(free_neighbors))
            else:
                node = random.choice(list(free_nodes))
        else:
            node = random.choice(list(free_nodes))
    return v

def tclustering_square(g, cover, pclust):
    v = {}
    for i in g:
        v[i] = 0
    free_nodes, nb_veg_final, nb_veg = set(list(g)), int(len(g)*cover), 0
    node = random.choice(list(free_nodes))
    while nb_veg < nb_veg_final:
        v[node] = 1#random.random()
        free_nodes.remove(node)
        nb_veg += 1
        if random.random() < pclust:
            i,j = node
            free_neighbors = {(i,j+1), (i,j-1), (i-1,j), (i+1,j), (i-1,j+1), (i+1,j+1), (i+1,j-1), (i-1,j+1)} & free_nodes
            if free_neighbors:
                node = random.choice(list(free_neighbors))
            else:
                node = random.choice(list(free_nodes))
                #free_neighbors = {(i-2,j+2), (i-1,j+2), (i,j+2), (i+1,j+2), (i+2,j+2), (i-2,j+1), (i+2,j+1), (i-2,j), (i+2,j), (i-2,j-1), (i+2,j-1), (i-2,j-2), (i-1,j-2), (i,j-2), (i+1,j-2), (i+2,j+2)} & free_nodes
                #if free_neighbors:
                #    node = random.choice(list(free_neighbors))
                #else:
                #    node = random.choice(list(free_nodes))
        else:
            node = random.choice(list(free_nodes))
    return v

def tclustering_square_from_px(g, vpx, pclust):
    v = {}
    for i in g:
        v[i] = 0
    free_nodes, vsum, vsum_final, vpx_list = set(list(g)), 0, sum(list(vpx.values())), list(i for i in vpx.values() if i>0)
    node = random.choice(list(free_nodes))
    while vsum < vsum_final:
        v[node] = random.choice(vpx_list)
        vsum += v[node]
        free_nodes.remove(node)
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

def tclustering_from_px(g, vx, pclust):
    v, veg_px = {}, list(i for i in vx.values() if i>0)
    for i in g:
        v[i] = 0
    random.shuffle(veg_px)
    free_nodes, nb_veg_final, nb_veg = list(g), sum(i>0 for i in veg_px), 0
    node = random.choice(free_nodes)
    while nb_veg < nb_veg_final:
        v[node] = veg_px[nb_veg]
        free_nodes.remove(node)
        nb_veg += 1
        if random.random() < pclust:
            neighbors = set(G.predecessors(node))|set(G.successors(node))
            free_neighbors = neighbors & set(free_nodes)
            if free_neighbors:
                node = random.choice(list(free_neighbors))
            else:
                node = random.choice(list(free_nodes))
        else:
            node = random.choice(list(free_nodes))
    return v

def estimate_weibull(sample, x1, x2):
    gamma = sp.symbols('gamma')
    GAMMA = sum(log(i)*i**(1/gamma) for i in sample)/sum(i**(1/gamma) for i in sample) - sum(log(i) for i in sample)/len(sample) - gamma
    dy = float(GAMMA.subs({gamma:x1})) - float(GAMMA.subs({gamma:x2}))
    dy_prev, k = dy+1, 0
    while dy < dy_prev:# is spacing around y=0 decreasing or not, if not, leave
        dy_prev = dy
        k += 1
        gamma_range = np.linspace(x1,x2+1e-3,10)
        Y = [float(GAMMA.subs({gamma:i})) for i in gamma_range]
        index_old = np.argmin(np.abs(Y))
        x_old, y_old = gamma_range[index_old], Y[index_old]
        if index_old in range(1,9):
            index_new = (y_old*Y[index_old+1]<0)*(index_old+1) + (y_old*Y[index_old-1]<0)*(index_old-1)
        elif index_old == 9:
            index_new = 8
        else:
            index_new = 1
        x_new, y_new = gamma_range[index_new], Y[index_new]
        slope = (y_old-y_new)/(x_old-x_new)
        x1 = gamma_range[index_old]
        x2 = x1 - Y[index_old]/slope
        dy = abs(y_new-y_old)
    c = 2/(x_old+x_new)
    x0 = (sum(i**c for i in sample)/len(sample))**(1/c)
    return c, x0, k

def plot_weibull(x, x0, c):
    return (c*x**(c-1))/(x0**c) *np.exp(-(x/x0)**c)

def plot_weibull3(x, mu, x0, c):
    return (c*(x-mu)**(c-1))/(x0**c) *np.exp(-((x-mu)/x0)**c)



