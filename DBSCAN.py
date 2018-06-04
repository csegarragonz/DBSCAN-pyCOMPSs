# DBSCAN 4 PyCOMPSs
# carlos.segarra @ bsc.es

# Imports
from collections import defaultdict # General Imports
from ast import literal_eval
import itertools, time, sys, os, argparse
from collections import deque
import numpy as np
from pycompss.api.task import task # PyCOMPSs Imports
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.api import compss_delete_object
from classes.Data import Data # DBSCAN Imports
from classes import constants
from classes.DS import DisjointSet
from classes.square import Square

def sync_relations(cluster_rules):
    out = defaultdict(set)
    for comb in cluster_rules:
        for key in cluster_rules[comb]:
            out[key] |= cluster_rules[comb][key]
    mf_set = DisjointSet(out.keys())
    for key in out:
        tmp = list(out[key])
        for i in range(len(tmp)-1):
            mf_set.union(tmp[i], tmp[i+1])
    return mf_set.get()

# DBSCAN Algorithm
def DBSCAN(epsilon, min_points, datafile, is_mn, print_times, *args, **kwargs):
    initial_time = time.time()

    if is_mn:
        TH_1=11000
        TH_2=1100
    else:
        TH_1=10000
        TH_2=100

    # Initial Definitions (necessary?)
    dataset_info = "dataset.txt"
    count_tasks = 0
    count_tasks_2 = 0

    # Data inisialitation
    dimensions = []
    f = open(dataset_info, "r")
    for line in f:
        split_line = line.split()
        if int(split_line[0]) == datafile:
            dimensions = literal_eval(split_line[1])
            break
    dimension_perms = [range(i) for i in dimensions]
    dataset = defaultdict()
    links = defaultdict()

    for comb in itertools.product(*dimension_perms):
        dataset[comb] = Square(comb, epsilon, dimensions)
        count_tasks_2 += dataset[comb].init_data(datafile, is_mn, TH_1,count_tasks_2)
        count_tasks += dataset[comb].partial_scan(min_points, TH_1, count_tasks)

    if print_times:
        compss_barrier()
        print "Partial Scan Tasks Finished"
        ps_time = time.time() - initial_time
        print "PS Lasted: "+str(ps_time - di_time)
        cp_count = 0
        for comb in itertools.product(*dimension_perms):
            dataset[comb].core_points = compss_wait_on(dataset[comb].core_points)
            cp_count += dataset[comb].core_points.count(constants.CORE_POINT)
        print "Number of core points found: "+str(cp_count)

    # We loop again since the first loop initialises all the variables that
    # are used afterwards
    for comb in itertools.product(*dimension_perms):
        neigh_sq_id = dataset[comb].neigh_sq_id
        labels_versions = []
        for neigh_comb in neigh_sq_id:
            # We obtain the labels found for our points by our neighbours.
            labels_versions.append(dataset[neigh_comb].cluster_labels[comb])
        links[comb] = dataset[comb].sync_labels(*labels_versions)
        count_tasks += 1

    for comb in itertools.product(*dimension_perms):
        links[comb] = compss_wait_on(links[comb])
        count_tasks += 1

    updated_links = sync_relations(links)
    for comb in itertools.product(*dimension_perms):
        dataset[comb].update_labels(updated_links, is_mn, datafile)

    print "Total number of tasks scheduled: "+str(count_tasks)
    print "Number of clusters found: "+str(len(updated_links))
    print "Time elapsed: " + str(time.time()-initial_time)

    return 1
#        dataset_tmp[comb] = Data()
#        len_datasets[comb] = count_lines(comb, datafile, is_mn)
#        adj_mat[comb] = [0]
#        tmp_mat[comb] = [0]
#        border_points[comb] = defaultdict(list)
#        fut_list, count_tasks = orquestrate_init_data(comb, datafile,
#                                                      len_datasets[comb], 1,
#                                                      0, [], TH_1, is_mn,
#                                                      count_tasks)
#        count_tasks += 1
#        # dataset[comb].points[comb]
#        dataset_tmp[comb] = merge_task_init(*fut_list)
#        neigh_sq_coord[comb] = neigh_squares_query(comb, epsilon,
#                                                   dimensions)

    if print_times:
        compss_barrier()
        print "Data Inisialisation Tasks Finished"
        di_time = time.time() - initial_time
        print "DI Lasted: "+str(di_time)

    # Partial Scan And Initial Cluster merging
#    for comb in itertools.product(*dimension_perms):
#        count_tasks += dataset[comb].partial_scan(min_points, TH_1, count_tasks)
#        cluster_labels = compss_wait_on(dataset[comb].cluster_labels)
#        relations = compss_wait_on(dataset[comb].relations)
#        print cluster_labels
#        print relations
##        adj_mat[comb] = [0]
##        tmp_mat[comb] = [0]
##        neigh_squares = []
##        for coord in neigh_sq_coord[comb]:
##            neigh_squares.append(dataset_tmp[coord])
##        [fut_list_0,
##         fut_list_1,
##         count_tasks] = orquestrate_scan_merge(dataset_tmp[comb], epsilon,
##                                               min_points, len_datasets[coord],
##                                               1, 0, [[], []], TH_1, count_tasks,
##                                               *neigh_squares)
##        count_tasks += 2
##        dataset[comb],tmp_mat[comb],adj_mat[comb] = merge_task_ps_0(*fut_list_0)
##        border_points[comb] = merge_task_ps_1(*fut_list_1)
#    return


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
        [fut_list,
         count_tasks] = orquestrate_sync_clusters(dataset[comb], adj_mat[comb],
                                                  epsilon, comb,
                                                  neigh_squares_loc,
                                                  len_neighs, 1, 0, [], TH_2,
                                                  count_tasks, *neigh_squares)
        count_tasks += 1
        adj_mat[comb] = merge_task_sync(adj_mat[comb], *fut_list)

    # Cluster list update
    border_points = dict_compss_wait_on(border_points, dimension_perms)
    tmp_mat = dict_compss_wait_on(tmp_mat, dimension_perms)
    adj_mat = dict_compss_wait_on(adj_mat, dimension_perms)

    if print_times:
        print "Cluster Synchronisation Tasks Finished"
        cs_time = time.time() - initial_time
        print "CS Lasted: "+str(cs_time - ps_time)

    links_list = unwrap_adj_mat(tmp_mat, adj_mat, neigh_sq_coord,
                                dimension_perms)
    for comb in itertools.product(*dimension_perms):
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset[coord])
        count_tasks += 1
        expand_cluster(dataset[comb], epsilon, border_points[comb],
                       dimension_perms, links_list, comb, tmp_mat[comb],
                       datafile, is_mn, *neigh_squares)
    return 1

if __name__ == "__main__":
    time.sleep(5)
    parser = argparse.ArgumentParser(description='DBSCAN Clustering Algorithm implemented within the PyCOMPSs framework. For a detailed guide on the usage see the user guide provided.')
    parser.add_argument('epsilon', type=float, help='Radius that defines the maximum distance under which neighbors are looked for.')
    parser.add_argument('min_points', type=int, help='Minimum number of neighbors for a point to be considered core point.')
    parser.add_argument('datafile', type=int, help='Numeric identifier for the dataset to be used. For further information see the user guide provided.')
    parser.add_argument('--is_mn', action='store_true', help='If set to true, this tells the algorithm that you are running the code in the MN cluster, setting the correct paths to the data files and setting the correct parameters. Otherwise it assumes you are running the code locally.')
    parser.add_argument('--print_times', action='store_true', help='If set to true, the timing for each task will be printed through the standard output. NOTE THAT THIS WILL LEAD TO EXTRA BARRIERS INTRODUCED IN THE CODE. Otherwise only the total time elapsed is printed.')
    args = parser.parse_args()
    DBSCAN(**vars(args))
