import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from copy import deepcopy
from itertools import product, combinations
from collections import Counter
from scipy.stats import pearsonr

def matrix_to_dict(matrix: np.ndarray):
    '''
    Convert matrix data into dictionary to be used to create a networkx 2D directed lattice for later use.
    '''
    d, height, width = {}, matrix.shape[0], matrix.shape[1]
    for i in range(height):
        for j in range(width):
            d[(j,59-i)] = matrix[i,j]
    return d

def generate_baseline_lattice(width: int, height: int):
    '''
    Create a baseline directed 2D lattice where all lateral and downwards directions are allowed.
    Inputs:
        width: number of nodes per row, type=int
        height: number of nodes per column, type=int
    Outputs
        g: 2D lattice, type=nx.diGraph
        bulk_: list of non-outflow nodes, type=list
    '''
    g = nx.grid_2d_graph(width,height)
    g.remove_edges_from([((i1,0),(i1+1,0)) for i1 in range(width-1)])
    g = g.to_directed()
    g.remove_edges_from([((i1,i2),(i1,i2+1)) for i1,i2 in product(range(width), range(height-1))])
    g.add_edges_from([((0,i2),(width-1,i2)) for i2 in range(1,height)]+[((width-1,i2),(0,i2)) for i2 in range(1,height)])
    bulk_ = [i for i in g if i not in [(i, 0) for i in range(width)]]
    return g, bulk_
    
def adapt_lattice_to_topography(g, t, bulk_nodes): # This seems to work, g_prob has the correct number of edges with their normalised slopes as edge weights
    '''
    Create lattices adapted to the topography for probabilistic or deterministic descent dynamics, and computes the propagation probabilities out of each node, and slope differences for each edge for future simulations.
    Inputs
        g: baseline lattice, type=nx.DiGraph
        t: landscape topography, type=dict
    Outputs
        g_prob: lattice for probabilistic dynamics
        propagation_prob: normalised propagation probability (value) from each node (key) to its successors, for probabilistic dynamics
        edge_slope: normalised propagation probability (value) of each edge (key), for probabilistic dynamics
        g_deter: lattice for deterministic dynamics
        propagation_prob_deter: normalised propagation probability (value) from each node (key) to its successors, for deterministic dynamics (trivial)
        edge_slope_deter: normalised propagation probability (value) of each edge (key), for deterministic dynamics (trivial)
    '''
    g_prob, g_deter, propagation_prob, edge_slope = nx.DiGraph(), nx.DiGraph(), {}, {}
    g_prob.add_nodes_from(g.nodes)
    g_deter.add_nodes_from(g.nodes)
    for i in bulk_nodes:
        next_nodes = [j for j in g.successors(i) if t[i]>t[j]]
        propagation_prob[i] = [t[i]-t[j] for j in next_nodes]
        propagation_prob[i] = [e/sum(propagation_prob[i]) for e in propagation_prob[i]]
        steepest_next = next_nodes[np.argmax(propagation_prob[i])]
        g_deter.add_edge(i, steepest_next)
        for k in range(len(next_nodes)):
            edge_slope[i, next_nodes[k]] = propagation_prob[i][k]
        if g.out_degree(i)==1:
            propagation_prob[i] == [1]
            edge_slope[i, list(g.successors(i))[0]] = 1
    for i,j in edge_slope:
        g_prob.add_edge(i, j, weight=edge_slope[i,j])
    propagation_prob_deter = {i:[1] for i in bulk_nodes}
    edge_slope_deter = {(i,j):1 for i,j in g_deter.edges}
    return g_prob, propagation_prob, edge_slope, g_deter, propagation_prob_deter, edge_slope_deter

def identify_steepest_direction(g_prob: nx.DiGraph, t: dict):
    '''
    Identify the single successor of each node according to the steepest descent direction.
    Inputs
        g_prob: landscape lattice assuming probabilistic dynamics, type=nx.DiGraph
        t: landscape topography, type=dict
    Output
        steep: associates each node (key) to its successor (value) according to the steepest descent direction, type=dict
    '''
    steep = {}
    for i in g_prob:
        next_nodes = list(g_prob.successors(i))
        if len(next_nodes):
            elev_nodes = np.array([t[k] for k in next_nodes])
            steep[i] = next_nodes[np.argmin(elev_nodes)]
    return steep

def probabilistic_to_deterministic(g_prob: nx.DiGraph, steep: dict):
    '''
    Adapting the landscape lattice from probabilistic to deterministic dynamics.
    Inputs
        g_prob: landscape lattice assuming probabilistic dynamics, type=nx.DiGraph
        steep: associates each node (key) to its successor (value) according to the steepest descent direction, type=dict
    Outputs
        g_deter: landscape lattice assuming deterministic dynamics, type=nx.DiGraph
    '''
    g_deter = deepcopy(g_prob)
    g_deter.remove_edges_from(list(g_prob.edges))
    for i in steep:
        _ = g_deter.add_edge(i, steep[i])
    return g_deter

### Simulation functions

def declare_simulation_variables(g: nx.DiGraph):
    '''
    Initialise all the variables required for a sandpile simulation.
    Input
        g: landscape lattice, type=nx.DiGraph
    Outputs
        state_: number of particles per node, type=dict
        coupling_: distances covered by particles exchanged between nodes, type=dict
        current_av_: nodes participating in the currently occuring avalanche, empty during accumulation phase, type=list
        branches_: length of branches described by the currently occuring avalanche, type=dict
        new_av_: nodes having received particles during the last time step, type=list
        size_: list of avalanche sizes having occured during the entire simulation, type=list
        exit_: number of particles having exited the landscae during the entire simulation, type=int
    '''
    state_, current_av_, branches_, new_active_, size_ = {i:0 for i in g}, {}, {}, [], []
    return state_, current_av_, branches_, new_active_, size_

def create_node_coupling_dictionary(g: nx.DiGraph):
    '''
    Identify the existence of paths allowing particle exchange between source (key) and target (value) nodes.
    This function is computationally expensive and only required for connectivity analysis, i.e. computing SC and FC.
    Input
        g: landscape lattice, type=nx.DiGraph
    Output
        coupling: tracks coupling between pairs of nodes (path existence), type=dict
    '''
    coupling, ordlist = {i:{} for i in g}, sorted(list(g), key=lambda a:a[1], reverse=True)
    for i,j in combinations(ordlist, 2):
        if nx.has_path(g,i,j):
            coupling[i][j] = []
        if i[1] == j[1]:		# for tgrid
            if nx.has_path(g,j,i):
                coupling[j][i] = []
    return coupling

def compute_slope_differences(g: nx.DiGraph, t: dict):
    '''
    Compute two variables required for probabilistic simulations (propagation_per_node) and the computation of structural connectivity (edge_slope)
    Inputs
        g: landscape lattice, type=nx.DiGraph
        t: landscape topography, type=dict
    Outputs
        propagation_probability: associate to each node (key) its list of propagation probabilities (value), type=dict
        edge_slope: associate to each edge (key) its elevation difference (value), type=dict
    '''
    propagation_probability, edge_slope = {i:[t[i]-t[j] for j in g.successors(i) if t[i]>t[j]] for i in g}, {(e1,e2):[0] for e1,e2 in g.edges}
    for i in g:
        edge_slope[i] = [edge_slope[i][j]/sum(edge_slope[i].values()) for j in edge_slope[i]]
        if g.out_degree(i)==1:
            propagation_probability[i] = [1]
            if list(edge_slope[i]) == 1:
                edge_slope[(i, list(g.successors(i))[0])] = 1
    return propagation_probability, edge_slope

def compute_SC(g: nx.DiGraph, edge_slope: dict):
    '''
    Compute structural connectivity
    Inputs
        g: landscape lattice, type=nx.DiGraph
        edge_slope: associate to each edge (key) its elevation difference (value), type=dict
    Outputs
        sc: landscape structural connectivity, type=dict
    '''
    gconn = deepcopy(g)
    gconn.remove_edges_from(list(g.edges()))
    for j in range(60):
        next_layer_successors, same_layer_successors = {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]-1] for v1 in [(i,j) for i in range(20)]}, {v1:[v2 for v2 in g.successors(v1) if v2[1]==v1[1]] for v1 in [(i,j) for i in range(20)]}
        for i in range(20):
            v1 = (i,j)
            gconn.add_weighted_edges_from([(v1,v2,edge_slope[(v1,v2)]) for v2 in next_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,gconn[v2][v3]['weight']) for v3 in gconn.successors(v2)] for v2 in next_layer_successors[v1]],[]))
        for v1 in sort_layer(same_layer_successors):
            gconn.add_weighted_edges_from([(v1,v2,edge_slope[(v1,v2)]) for v2 in same_layer_successors[v1]])
            gconn.add_weighted_edges_from(sum([[(v1,v3,gconn[v2][v3]['weight']) for v3 in gconn.successors(v2)] for v2 in same_layer_successors[v1]],[]))
    sc = {node:gconn.in_degree(node, weight='weight') for node in g}
    return sc

def sort_layer(nodeset):
    '''
    Sort subset of nodes to create a more computational efficient loop.
    Input
        nodeset: associate a node and its successors (subset of nodes), type=dict
    Output
        sorted_nodeset: sorted input, type=list
    '''
    gl = nx.DiGraph()
    gl.add_nodes_from(nodeset.keys())
    gl.add_edges_from([(i,nodeset[i][0]) for i in nodeset if len(nodeset[i])])
    sorted_nodeset = list(nx.topological_sort(gl))[::-1]
    return sorted_nodeset

def compute_FC(coupling):
    '''
    Compute FC
    Input
        coupling: tracked particle exchange (list of integers) between a source (key) and target nodes (value), type=dict
    '''
    fc = {}
    for i in coupling:
        fci = sum(coupling[i].values(),[])
        if len(fci):
            fc[i] = np.var(fci)
        else:
            fc[i] = 0
    return fc

def dict_to_mat(x):
    '''
    Convert a (60,20) dictionary of nodes into a matrix for a better visualisation of the landscape map.
    Input
        x: input dictionary of nodes of a (60,20) lattice, dtype=dict
    Output
        mat: output matrix, dtype=np.ndarray of shape (60,20)
    '''
    mat = np.zeros((60,20), float)
    for i in range(60):
        for j in range(20):
            mat[i,j] = x[(j,59-i)]
    return mat

# Data generation functions

def generate_vegetation(g: nx.DiGraph, ratio_v: float, pclust: float, v_source: dict):
    '''
    Generate vegetation distribution.
    Inputs
        g: landscape lattice, type=nx.Graph
        ratio_v: ratio of vegetated nodes, type=float between 0 and 1
        pclust: clustering probability, type=float between 0 and 1
        v_source: empirical vegetation to sample vegetation densities from, type=dict
    Output
        v: vegetation density (value) per node (key), type=dict
    '''
    v = {i:0 for i in g}
    free_nodes, nb_veg_final, nb_veg, v_sample = set(list(g)), int(len(g)*ratio_v), 0, list(i for i in v_source.values() if i>0)
    node = random.choice(list(free_nodes))
    while nb_veg < nb_veg_final:
        v[node] = random.choice(v_sample)
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

def compute_landscape_vegetation_score(vmat: np.ndarray):
    '''
    Compute the vegetation scores (vegetation density neighborhood) of each node of an empirical landscape.
    Input
        vmat: vegetation matrix of empirical landscape (original imported field data), type=np.ndarray
    Output
        land_scores: matrix of vegetation scores
    '''
    land_scores = np.zeros((60,20), dtype=float)
    for j in range(1,19):
        for i in range(1,59):
            land_scores[i,j] = 1*vmat[i,j] + vmat[i-1,j] + vmat[i+1,j] + vmat[i,j-1]+ vmat[i,j+1]
            land_scores[i,0] = 2*vmat[i,0] + vmat[i-1,0] + vmat[i+1,0] + vmat[i,1]
            land_scores[i,19] = 2*vmat[i,19] + vmat[i-1,19] + vmat[i+1,19] + vmat[i,18]
        land_scores[0,j] = 2*vmat[0,j] + vmat[0,j-1] + vmat[0,j+1] + vmat[1,j]
        land_scores[59,j] = 2*vmat[59,j] + vmat[59,j-1] + vmat[59,j+1] + vmat[58,j]
    land_scores[0,0] = 3*vmat[0,0] + vmat[0,1] + vmat[1,0]
    land_scores[59,0] = 3*vmat[59,0] + vmat[59,1] + vmat[58,0]
    land_scores[0,19] = 3*vmat[0,19] + vmat[0,18] + vmat[1,19]
    land_scores[59,19] = 3*vmat[59,19] + vmat[58,19] + vmat[59,18]
    land_scores /= 100
    return land_scores

def generate_topography(g: nx.DiGraph, tmat: dict, v: dict, a: float, a_std: float, b: float):
    '''
    Generate a topography
    Inputs
        g: baseline lattice, type=nx.DiGraph
        tmax: empirical landscape topography matrix to be used to generate planar elevation, type=dict
        v: generated vegetation, type=dict
        a: slope of the linear regression between vegetation density scores and microtopography, type=float
        a_std: standard deviation of the slope a, type=float
        b: intercept of the linear regression between vegetation density scores and microtopography, type=float
    Output
        t_gen: generated topography, type=dict
    '''
    v_score, g_x, g_y, plane = {}, list(g)[-1][0], list(g)[-1][1], np.sum(tmat,axis=1)/20
    extended_plane = sum([list(plane+k*(1-.8*plane[-1])) for k in reversed(range(int(g_y/60)+1))],[])
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
    t_gen = {(i,j): v_score[i,j]*a + (2*random.random()-1)*a_std/2 + b + extended_plane[j] for i,j in g}
    return t_gen

# Simulation functions

def sandpile(g, v, propagation_probability, state, coupling, current_av, branches, old_active, size, bulk_nodes):
    unstable, new_active = [node for node in old_active if state[node] >= g.out_degree(node)], []
    if len(unstable):
        # Relaxation phase
        for node1 in unstable:
            nb_partcl, state[node1] = state[node1], 0
            if len(current_av) == 0:
                current_av[node1] = 0
            if len(list(g.successors(node1))):
                spreading_scheme = Counter(random.choices(list(g.successors((node1))), weights=propagation_probability[(node1)], k=nb_partcl))
                for node2 in spreading_scheme:
                    if g.out_degree(node2):
                        state[node2] += np.sum(np.random.rand(spreading_scheme[node2]) > v[node2]/2)
                        current_av[node2] = current_av[node1]+1
                        if node2 not in branches:
                            branches[node2] = [node1]
                        elif node1 not in branches[node2]:
                            branches[node2].append(node1)
                        # Record all the previous sources of the new incident node (at least, from a new path)
                        origin, all_sources, sources = list(current_av.keys())[0], [node1], [node1]
                        while origin not in sources:
                            new_sources = []
                            for source in sources:
                                new_sources += branches[source]
                            all_sources += new_sources
                            sources = new_sources
                        # Record couplings between new incident node and all its previous sources
                        for source in list(set(all_sources)): # collapse 'all_sources' to avoid redundance of nodes
                            coupling[source][node2].append(current_av[node2]-current_av[source])
                    if node2 not in new_active:
                        new_active.append(node2)
    else:
        # Accumulation phase
        if len(current_av):
            size.append(len(current_av))
        current_av, branches = {}, {}
        node2 = random.choice(bulk_nodes)
        if random.random() > v[node2]/2:
            state[node2] += 1
            new_active.append(node2)
    return state, coupling, current_av, branches, new_active, size

def SCFC(sc, fc, bulk_):
    return pearsonr(list(sc[i] for i in bulk_), list(fc[i] for i in bulk_))[0]

