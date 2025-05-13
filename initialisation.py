import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from deepcopy import deepcopy
from itertools import product, combinations

def matrix_to_dict(matrix: np.ndarray):
    '''
    Convert matrix data into dictionary to be used to create a networkx 2D directed lattice for later use.
    '''
    #print(matrix.shape)
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
        bulk_nodes: list of non-outflow nodes, type=list
    '''
    g = nx.grid_2d_graph(width,height)
    g.remove_edges_from([((i1,0),(i1+1,0)) for i1 in range(width-1)])
    g = g.to_directed()
    g.remove_edges_from([((i1,i2),(i1,i2+1)) for i1,i2 in product(range(width), range(height-1))])
    g.add_edges_from([((0,i2),(width-1,i2)) for i2 in range(1,height)]+[((width-1,i2),(0,i2)) for i2 in range(1,height)])
    bulk_nodes = [i for i in g if i not in [(i, 0) for i in range(width)]]
    return g, bulk_nodes

def adapt_lattice_to_topography(g: nx.DiGraph, t: dict, dyn: str):
    '''
    Adapt lattice to topography for probabilistic or deterministic descent sandpile simulations.
    Inputs
        g: baseline lattice, type=nx.diGraph
        t: landscape topography, type=dict
        dyn: probabilistic or deterministic dynamics selection, type=str
    Outputs
        g_prob: lattice adapted to topography
    '''
    g_prob = nx.DiGraph()
    g_prob.add_nodes_from(g.nodes)
    for i in g:
        for j in g.successors(i):
            if e[i] >= e[j]:
                g_prob.add_edge(i,j)
    if dyn == 'probabilistic':
        return g_prob
    elif dyn == 'deterministic':
        steep = identify_steepest_direction(g_prob, t)
        return probabilistic_to_deterministic(g_prob, steep)
    else:
        print('Error: Wrong type of dynamics, please enter either \'probabilistic\' or \'deterministic\' (third input variable).')

def identify_steepest_direction(g_prob, t):
    steep = {}
    for i in g_prob:
        next_nodes = list(g_prob.successors(i))
        if len(next_nodes):
            elev_nodes = np.array([t[k] for k in next_nodes])
            min_elev = np.where(elev_nodes==np.min(elev_nodes))[0]
            if len(min_elev)>1:
                steep[i] = (i[0],i[1]-1)
            else:
                steep[i] = next_nodes[min_elev[0]]
    return steep

def probabilistic_to_deterministic(g_prob, steep):
    g_deter = deepcopy(g_prob)
    g_deter.remove_edges_from(list(g_prob.edges))
    for i in steep:
        _ = g_deter.add_edge(i, steep[i])
    return g_deter

### Simulation functions
def declare_simulation_variables(g):
    '''
    Initialise all the variables required for a sandpile simulation.
    Input
        g: landscape lattice, type=nx.diGraph
    Outputs
        state_: number of particles per node, type=dict
        coupling_: distances covered by particles exchanged between nodes, type=dict
        current_av_: nodes participating in the currently occuring avalanche, empty during accumulation phase, type=list
        branches_: length of branches described by the currently occuring avalanche, type=dict
        new_av_: nodes having received particles during the last time step, type=list
        size_: list of avalanche sizes having occured during the entire simulation, type=list
    '''
    state_, coupling_, current_av_, branches_, new_, size_, ordlist = {}, {}, [], {}, [], [], sorted(list(g), key=lambda a:a[1], reverse=True)
    for i in g:
        state_[i], coupling_[i] = 0, {}
    for i,j in combinations(ordlist, 2):
        if nx.has_path(g,i,j):
            coupling_[i][j] = []
        if i[1] == j[1]:
            if nx.has_path(g,j,i):
                coupling_[j][i] = []
    return state_, coupling_, current_av_, branches_, new_active_, size_

def declare_simulation_variables_for_criticality_only(g):
    '''
    Initialise all the variables required for a sandpile simulation.
    Input
        g: landscape lattice, type=nx.diGraph
    Outputs
        state_: number of particles per node, type=dict
        coupling_: distances covered by particles exchanged between nodes, type=dict
        current_av_: nodes participating in the currently occuring avalanche, empty during accumulation phase, type=list
        branches_: length of branches described by the currently occuring avalanche, type=dict
        new_av_: nodes having received particles during the last time step, type=list
        size_: list of avalanche sizes having occured during the entire simulation, type=list
    '''
    state_, coupling_, current_av_, branches_, new_active_, size_ = {}, {}, [], {}, [], [], sorted(list(g), key=lambda a:a[1], reverse=True)
    for i in g:
        state_[i], coupling_[i] = 0, {}
    return state_, current_av_, branches_, new_active_, size_

def compute_node_connectivity(g):
    '''
    Compute dictionary of node connectivity, to track particle exchange between node.
    This function is computationally expensive and only necessary for connectivity analysis, i.e. computing SC and FC.
    '''



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



