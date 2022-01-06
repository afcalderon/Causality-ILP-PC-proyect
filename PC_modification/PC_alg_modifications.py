#!/usr/bin/env python
# -*- coding: utf-8 -*-
# +
# This functions are modifications of the PC algorithm implementation
# created by chufangao https://github.com/chufangao/pcalg
# -
from itertools import combinations, permutations
import logging

import networkx as nx
import numpy as np
import pandas as pd


# Estimate skeleton injecting the ground truth
def estimate_skeleton_true(adj, data_matrix, **kwargs):
    
    def method_stable(kwargs):
        return ('method' in kwargs) and kwargs['method'] == "stable"
    
    D = nx.DiGraph(adj, directed=True)
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
                if len(adj_i) < l:
                    continue
                for k in combinations(adj_i, l):
                    p_val = nx.d_separated(D, {i}, {j},set(k)) 
                    if p_val is True :
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

    return (g, sep_set)


# +
# Calculate CI test with different number of alpha  and save them in a dictionary

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


# -

# get the alpha in the dictionary create it in "estimate_skeleton_mulalpha"
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
