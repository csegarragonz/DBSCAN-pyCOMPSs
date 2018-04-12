# DBSCAN 4 PyCOMPSs
# This version is an attempt to divide merge_task_ps_0 in smaller tasks.
# The reason why it has been captured in a different script is that previous 
# version is sufficiently good to consider not ruining it, as done in previous
# cases.
# carlos.segarra @ bsc.es

# IMPORTED MODULES
    # General Imports
from collections import defaultdict
from ast import literal_eval
import itertools
import time
from collections import deque
import numpy as np
import sys
import os
    # PyCOMPSs Imports
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.api import compss_delete_object
    # DBSCAN Imports
from classes.Data import Data
from task_src.init_data import count_lines, orquestrate_init_data
from task_src.init_data import init_data, merge_task_init, neigh_squares_query
from task_src.partial_scan import orquestrate_scan_merge, partial_scan_merge
from task_src.partial_scan import merge_task_ps_0, merge_task_ps_1
from task_src.dict_wait_on import dict_compss_wait_on
from task_src.sync_clusters import orquestrate_sync_clusters, sync_clusters
from task_src.sync_clusters import merge_task_sync
from task_src.expand_cluster import unwrap_adj_mat, expand_cluster

def DBSCAN(epsilon, min_points, file_id, TH_1, TH_2):
    #   TODO: code from scratch the Disjoint Set
    #   TODO: comment the code apropriately
    #   TODO: separate the code in different modules

    # DBSCAN Algorithm
    initial_time = time.time()

    # Initial Definitions (necessary?)
    epsilon = float(epsilon)
    min_points = int(min_points)
    dataset_info = "dataset.txt"

    # Data inisialitation
    dimensions = []
    f = open(dataset_info, "r")
    for line in f:
        split_line = line.split()
        if split_line[0] == file_id:
            dimensions = literal_eval(split_line[1])
            break
    dimension_perms = [range(i) for i in dimensions]
    dataset = defaultdict()
    dataset_tmp = defaultdict()
    neigh_sq_coord = defaultdict()
    len_datasets = defaultdict()
    adj_mat = defaultdict()
    tmp_mat = defaultdict()
    border_points = defaultdict()
    for comb in itertools.product(*dimension_perms):
        # TODO: implement init_data as a class method maybe inside
        # the initialisation
        dataset[comb] = Data()
        dataset_tmp[comb] = Data()
        len_datasets[comb] = count_lines(comb, file_id)
        adj_mat[comb] = [0]
        tmp_mat[comb] = [0]
        border_points[comb] = defaultdict(list)
        fut_list = orquestrate_init_data(comb, file_id, len_datasets[comb], 1,
                                         0, [], TH_1)
        dataset_tmp[comb] = merge_task_init(*fut_list)
        neigh_sq_coord[comb] = neigh_squares_query(comb, epsilon,
                                                   dimensions)

    compss_barrier()
    print "Data Inisialisation Tasks Finished"
    di_time = time.time() - initial_time
    print "DI Lasted: "+str(di_time)

    # Partial Scan And Initial Cluster merging
    for comb in itertools.product(*dimension_perms):
        adj_mat[comb] = [0]
        tmp_mat[comb] = [0]
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset_tmp[coord])
        fut_list_0, fut_list_1 = orquestrate_scan_merge(dataset_tmp[comb],
                                                        epsilon, min_points,
                                                        len_datasets[coord],
                                                        1, 0, [[], []], TH_1,
                                                        *neigh_squares)
        dataset[comb],tmp_mat[comb],adj_mat[comb] = merge_task_ps_0(*fut_list_0)
        border_points[comb] = merge_task_ps_1(*fut_list_1)

    compss_barrier()
    print "Partial Scan Tasks Finished"
    ps_time = time.time() - initial_time
    print "PS Lasted: "+str(ps_time - di_time)

    # Cluster Synchronisation
    for comb in itertools.product(*dimension_perms):
        compss_delete_object(dataset_tmp[comb])
        neigh_squares = []
        neigh_squares_loc = []
        len_neighs = 0
        for coord in neigh_sq_coord[comb]:
            neigh_squares_loc.append(coord)
            neigh_squares.append(dataset[coord])
            len_neighs += len_datasets[coord]
        fut_list = orquestrate_sync_clusters(dataset[comb], adj_mat[comb],
                                             epsilon, comb, neigh_squares_loc,
                                             len_neighs, 1, 0, [], TH_2,
                                             *neigh_squares)
        adj_mat[comb] = merge_task_sync(adj_mat[comb], *fut_list)

    # Cluster list update
    border_points = dict_compss_wait_on(border_points, dimension_perms)
    tmp_mat = dict_compss_wait_on(tmp_mat, dimension_perms)
    adj_mat = dict_compss_wait_on(adj_mat, dimension_perms)

    print "Cluster Synchronisation Tasks Finished"
    cs_time = time.time() - initial_time
    print "CS Lasted: "+str(cs_time - ps_time)

    links_list = unwrap_adj_mat(tmp_mat, adj_mat, neigh_sq_coord,
                                dimension_perms)
    for comb in itertools.product(*dimension_perms):
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset[coord])
        expand_cluster(dataset[comb], epsilon, border_points[comb],
                       dimension_perms, links_list, comb, tmp_mat[comb],
                       file_id, *neigh_squares)
    print "Time elapsed: " + str(time.time()-initial_time)
    return 1

if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]),
            int(sys.argv[5]))
