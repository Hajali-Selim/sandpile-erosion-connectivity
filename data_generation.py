import random

def generate_vegetation(g: nx.DiGraph, cover: float, pclust: float, v_source: dict):
    '''
    Generate vegetation distribution.
    Inputs
        g: landscape lattice, type=nx.Graph
        cover: ratio of vegetated nodes, type=float between 0 and 1
        pclust: clustering probability, type=float between 0 and 1
        v_source: empirical vegetation to sample vegetation densities from, type=dict
    Output
        v: vegetation density (value) per node (key), type=dict
    '''
    v = {i:0 for i in g}
    free_nodes, nb_veg_final, nb_veg, v_sample = set(list(g)), int(len(g)*cover), 0, list(i for i in v_source.values() if i>0)
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
    return v

def compute_landscape_vegetation_score(vmat):
    '''
    Compute the vegetation scores (vegetation density neighborhood) of each node of an empirical landscape.
    Input
        vmat: vegetation matrix of empirical landscape (original imported field data), type=np.ndarray
    Output
        land_score: matrix of vegetation density scores
    '''
    land_score = np.zeros((60,20), dtype=float)
    for j in range(1,19):
        for i in range(1,59):
            land_score[i,j] = 1*vmat[i,j] + vmat[i-1,j] + vmat[i+1,j] + vmat[i,j-1]+ vmat[i,j+1]
            land_score[i,0] = 2*vmat[i,0] + vmat[i-1,0] + vmat[i+1,0] + vmat[i,1]
            land_score[i,19] = 2*vmat[i,19] + vmat[i-1,19] + vmat[i+1,19] + vmat[i,18]
        land_score[0,j] = 2*vmat[0,j] + vmat[0,j-1] + vmat[0,j+1] + vmat[1,j]
        land_score[59,j] = 2*vmat[59,j] + vmat[59,j-1] + vmat[59,j+1] + vmat[58,j]
    land_score[0,0] = 3*vmat[0,0] + vmat[0,1] + vmat[1,0]
    land_score[59,0] = 3*vmat[59,0] + vmat[59,1] + vmat[58,0]
    land_score[0,19] = 3*vmat[0,19] + vmat[0,18] + vmat[1,19]
    land_score[59,19] = 3*vmat[59,19] + vmat[58,19] + vmat[59,18]
    return land_score/100

