#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import combinations, permutations
import logging

import networkx as nx
import numpy as np
import pandas as pd


# +
# generate a random adjacency matrix
# traces: Number or Domino Traces
# If traces>1 the output will be a data frame of list
# nodes: Number of nodes 
# parent_max: number max of possible parents per node

def adjacency_generator(traces,nodes,parents_max):
        
    def connected_graph(adjacency_matrix):
         g = nx.DiGraph(adjacency_matrix)
         connect = nx.is_weakly_connected(g)
         return connect
    
    data = []
    
    for k in range(traces):
        
        finished = False
        
        while not finished:
        
            permutation = np.random.permutation(range(0,nodes)) #creating permutation matrix
            idx = np.empty_like(permutation)
            idx[permutation] = np.arange(len(permutation))       
            adjacency_matrix = np.zeros((nodes, nodes),dtype=int)

            for j in range(1, nodes):
                nb_parents = np.random.randint(0, min([parents_max, j])+1) # selecting number of parents for each node min 1
                for i in np.random.choice(range(0,j), nb_parents, replace=True):# selecting randomly connections betwwen nodes
                    adjacency_matrix[i, j] = 1 

            adjacency_matrix[:, idx]
            adjacency_matrix[:] = adjacency_matrix[:, idx]
            adjacency_matrix[:] = adjacency_matrix[idx,:]
            
            finished = connected_graph(adjacency_matrix)
        data.append(adjacency_matrix) # generating nested list of adjacency matrix

    return data


# +
# Add causal relation between parents and childs nodes
# adjacency_matrix : a adjacency matrix it must be a numpy array
# nodes: Number of nodes 
# parent_max: number max of possible parents per node
# N : Number of experimental cases (sample) per domino trace
# p : probability that the initial root fall or not in the Domino trace
# eps: Quantity of noise level in the data

def theta_generator_multi(adjacency_matrix,nodes,parents_max, N, p, eps):
    

    data = [] # creating empty data sets
    num_rows = len(adjacency_matrix)

    for b in range(num_rows):
        
        matrix = adjacency_matrix[b]
        nodes = len(matrix)
        X = np.zeros(shape=(nodes, N), dtype=int) 
        
        for t in range (2):
            for i in range(0, nodes):
                if not sum(matrix[:, i]):
                    X[i,:] = np.random.binomial(1, p, size=N)
                    for k in range(0, nodes):
                        if sum(matrix[:, k]):
                            parents = np.where(matrix[:,k] == 1)[0]
                            X[k, :] = np.ones_like(X[k, :])
                            for a in parents:
                                noise = np.random.binomial(1, eps, size=N)
                                X[k, :] = X[k,:]*[(1-noise[j])*X[a,j]+noise[j]*(1-X[a,j]) for j in range(N)]

                    theta = X.sum(axis=1)/N                    
                   
        data.append({"Trace": matrix, "Theta": theta, "Matrix_theta": X}) # generating nested list of adjacency matrix
        
    df = pd.DataFrame(data=data).sample(frac=1).reset_index(drop=True)

    
    return df


# +
# Metrics of recall and precision of the skeleton
# between ground truth and predicted graph

def diff_DAG(dag, skel):
    
    dag_edges = set(list(dag.edges()))
    skel_edges = set(list(skel.edges()))
    
    dag_edges_inv = {(j, i) for i, j in dag_edges }
    edges_skel_inv = {(j, i) for i, j in skel_edges }
    
    additions = skel_edges - dag_edges - dag_edges_inv
    deletions = dag_edges - skel_edges - edges_skel_inv
    
    diff = len(additions) + len(deletions)
    
    true_positives = len(dag_edges) - len(deletions) 
    false_negatives = len(deletions) 
    false_positives = len(additions)
    
    if false_positives + true_positives != 0 and true_positives + false_negatives != 0 :
        precision = true_positives/(false_positives + true_positives)
        recall = true_positives/(true_positives + false_negatives)
    
    else : 
        precision = 0
        recall = 0         
    
    return precision, recall, len(dag_edges)


# +
# Metrics of recall and precision of the direction of edges
# between ground truth and predicted graph

def rec_directions(dag, pred_dag):
    
    dag_edges = set(list(dag.edges()))
    pred_dag_edges = set(list(pred_dag.edges())) 
    
    dag_edges_inv = {(j, i) for i, j in dag_edges }
    edges_pred_dag_inv = {(j, i) for i, j in pred_dag_edges }
    pred_dag_inv_diff = edges_pred_dag_inv - pred_dag_edges 
    
    additions = pred_dag_edges - dag_edges - pred_dag_inv_diff
    deletions = dag_edges - pred_dag_edges - dag_edges_inv
    
    true_positives = len(pred_dag_edges) - len(additions)
    false_positives = len(additions) 
    false_negatives = len(deletions)
    
    if false_positives + true_positives != 0 and true_positives + false_negatives != 0 :
        precision = true_positives/(false_positives + true_positives)
        recall = true_positives/(true_positives + false_negatives)
        
    else : 
        precision = 0
        recall = 0 
    
    return precision, recall


# -

def estimate_skeleton_mulalpha(indep_test_func, data_matrix, alpha, **kwargs):
    l_pval = []
    def method_stable(kwargs):
        return ('method' in kwargs) and kwargs['method'] == "stable"

    node_ids = range(data_matrix.shape[1])
    node_size = data_matrix.shape[1]
    sep_set = [[set() for i in range(node_size)] for j in range(node_size)]
    if 'init_graph' in kwargs:
        g = kwargs['init_graph']
        if not isinstance(g, nx.Graph):
            raise ValueError
        elif not g.number_of_nodes() == len(node_ids):
            raise ValueError('init_graph not matching data_matrix shape')
        for (i, j) in combinations(node_ids, 2):
            if not g.has_edge(i, j):
                sep_set[i][j] = None
                sep_set[j][i] = None
    else:
        g = _create_complete_graph(node_ids)

    l = 0
    print("multi")
    while True:
        cont = False
        remove_edges = []
        for (i, j) in permutations(node_ids, 2):
            adj_i = list(g.neighbors(i))
            if j not in adj_i:
                continue
            else:
                adj_i.remove(j)
            if len(adj_i) >= l:
                _logger.debug('testing %s and %s' % (i,j))
                _logger.debug('neighbors of %s are %s' % (i, str(adj_i)))
                if len(adj_i) < l:
                    continue
                
                
                for k in combinations(adj_i, l):
                    
                    p_val = indep_test_func(data_matrix, i, j, set(k),
                                            **kwargs)

                    l_pval.append({"i":i,"j":j,"set_k":set(k),"p_val":p_val})
                   
                    if p_val > alpha:
                        if g.has_edge(i, j):
                            if method_stable(kwargs):
                                remove_edges.append((i, j))
                            else:
                                g.remove_edge(i, j)
                        sep_set[i][j] |= set(k)
                        sep_set[j][i] |= set(k)
                        break
                cont = True
        l += 1
        if method_stable(kwargs):
            g.remove_edges_from(remove_edges)
        if cont is False:
            break
        if ('max_reach' in kwargs) and (l > kwargs['max_reach']):
            break
    df_pval = pd.DataFrame(data=l_pval).sample(frac=1).reset_index(drop=True)
    return (g, sep_set,df_pval )

def estimate_skeleton_list(data_matrix, alpha, l_pval, **kwargs):
    
    def method_stable(kwargs):
        return ('method' in kwargs) and kwargs['method'] == "stable"

    node_ids = range(data_matrix.shape[1])
    node_size = data_matrix.shape[1]
    sep_set = [[set() for i in range(node_size)] for j in range(node_size)]
    if 'init_graph' in kwargs:
        g = kwargs['init_graph']
        if not isinstance(g, nx.Graph):
            raise ValueError
        elif not g.number_of_nodes() == len(node_ids):
            raise ValueError('init_graph not matching data_matrix shape')
        for (i, j) in combinations(node_ids, 2):
            if not g.has_edge(i, j):
                sep_set[i][j] = None
                sep_set[j][i] = None
    else:
        g = _create_complete_graph(node_ids)

    l = 0
    while True:
        cont = False
        remove_edges = []
        for (i, j) in permutations(node_ids, 2):
            adj_i = list(g.neighbors(i))
            if j not in adj_i:
                continue
            else:
                adj_i.remove(j)
            if len(adj_i) >= l:
                _logger.debug('testing %s and %s' % (i,j))
                _logger.debug('neighbors of %s are %s' % (i, str(adj_i)))
                if len(adj_i) < l:
                    continue
                
                
                for k in combinations(adj_i, l):
                    
                    #p_val = indep_test_func(data_matrix, i, j, set(k),
                    p_val = l_pval.p_val[(l_pval.i == i) & (l_pval.j == j) & (l_pval.set_k == set(k))].values.item(0)
                
                    if p_val > alpha:
                        if g.has_edge(i, j):
                            if method_stable(kwargs):
                                remove_edges.append((i, j))
                            else:
                                g.remove_edge(i, j)
                        sep_set[i][j] |= set(k)
                        sep_set[j][i] |= set(k)
                        break
                cont = True
        l += 1
        if method_stable(kwargs):
            g.remove_edges_from(remove_edges)
        if cont is False:
            break
        if ('max_reach' in kwargs) and (l > kwargs['max_reach']):
            break
            
    return (g, sep_set )
