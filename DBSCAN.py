# DBSCAN 4 PyCOMPSs
# Version with nesting
# carlos.segarra @ bsc.es

# Imports
from pycompss.api.task import task
from pycompss.api.parameter import INOUT
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_delete_object
from collections import defaultdict
from classes.DS import DisjointSet
from classes.Data import Data
from ast import literal_eval
import copy
import itertools
import time
import numpy as np
import sys
import os


def count_lines(tupla, file_id):
    # path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_id)
    path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(tupla[0])
    for num, j in enumerate(tupla):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += ".txt"
    with open(tmp_string) as infile:
        for i, line in enumerate(infile):
            pass
        return i+1


@task(data_pos=INOUT)
def init_data(data_pos, tupla, file_id):
    # path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_id)
    path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(tupla[0])
    for num, j in enumerate(tupla):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += ".txt"
    data_pos.value = np.loadtxt(tmp_string)
    if len(np.shape(data_pos.value)) == 1:
        data_pos.value = np.array([data_pos.value])
    tmp_vec = -2*np.ones(np.shape(data_pos.value)[0])
    data_pos.value = [data_pos.value, tmp_vec]
#    return data_pos


def neigh_squares_query(square, epsilon, dimensions):
    # Only multiples of 10 are supported as possible number
    # of squares per side.
    dim = len(square)
    neigh_squares = []
    border_squares = [int(min(max(epsilon*i, 1), i-1)) for i in dimensions]
#    border_squares = int(max(epsilon/0.1, 1))
    perm = []
    for i in range(dim):
        perm.append(range(-border_squares[i], border_squares[i] + 1))
    for comb in itertools.product(*perm):
        current = []
        for i in range(dim):
            if square[i] + comb[i] in range(dimensions[i]):
                current.append(square[i]+comb[i])
        if len(current) == dim and current != list(square):
            neigh_squares.append(tuple(current))
    neigh_squares.append(tuple(square))
    return tuple(neigh_squares)


@task(returns=(Data, defaultdict, list, list, list))
def partial_scan_merge(data, epsilon, min_points, *args):
    # Core point location in the data chunk
    data_copy = Data()
    data_copy.value = copy.deepcopy(data.value)
    tmp_unwrap = [i.value[0] for i in args]
    tmp_unwrap_2 = [i.value[1] for i in args]
    neigh_points = np.vstack(tmp_unwrap)
    neigh_points_clust = np.concatenate(tmp_unwrap_2)
    non_assigned = defaultdict(int)
    dist_mat = np.linalg.norm(neigh_points -
                              data_copy.value[0][:, np.newaxis, :],
                              axis=2)-epsilon < 0
    for num in range(dist_mat.shape[0]):
        neigh_count = np.sum(dist_mat[num])
        if neigh_count > min_points:
            data_copy.value[1][num] = -1
        elif neigh_count > 0:
            tmp = []
            for pos, proxim in enumerate(dist_mat[num]):
                if proxim:
                    # Adding three since later to detect core points we will
                    # require value > -1 and -0 is > 1 and not to confuse it
                    # with a noise point
                    if neigh_points_clust[pos] == -1:
                        data_copy.value[1][num] = -(pos+3)
                        break
                    else:
                        tmp.append(pos)
            else:
                non_assigned[num] = tmp
    # Cluster the core points found
    cluster_count = 0
    core_points = [[num, p] for num, p in enumerate(data_copy.value[0]) if
                   data_copy.value[1][num] == -1]
    for pos, point in core_points:
        # TODO: make a one-liner out of this
        tmp_vec = (np.linalg.norm(data_copy.value[0]-point, axis=1)-epsilon) < 0
        for num, poss_neigh in enumerate(tmp_vec):
            if poss_neigh and data_copy.value[1][num] > -1:
                data_copy.value[1][pos] = data_copy.value[1][num]
                break
        else:
            data_copy.value[1][pos] = cluster_count
            cluster_count += 1
    return [data_copy, non_assigned, [cluster_count], [cluster_count],
            dist_mat]
#    return data, non_assigned


def dict_compss_wait_on(dicc, dimension_perms):
    for comb in itertools.product(*dimension_perms):
        dicc[comb] = compss_wait_on(dicc[comb])
    return dicc


def orquestrate_sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc,
                              dist_mat, len_neighs, quocient, res,
                              fut_list, *args):
    # TODO: currently hardcoded
    THRESHOLD = 100
    if (len_neighs/quocient) > THRESHOLD:
        orquestrate_sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc,
                                  dist_mat, len_neighs, quocient*2,
                                  res*2 + 0, fut_list, *args)
        orquestrate_sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc,
                                  dist_mat, len_neighs, quocient*2,
                                  res*2 + 1, fut_list, *args)
    else:
        fut_list.append(sync_clusters(data, adj_mat, epsilon, coord,
                                      neigh_sq_loc, dist_mat, quocient, res,
                                      *args))
    return fut_list
    # TODO: picar el merge
    # TODO: fer len_neigh accessible


@task(returns=list)
def merge_task(adj_mat, *args):
    adj_mat_copy = [[] for _ in range(max(adj_mat[0], 1))]
    for args_i in args:
        for num, list_elem in enumerate(args_i):
            for elem in list_elem:
                if elem not in adj_mat_copy[num]:
                    adj_mat_copy[num].append(elem)
    return adj_mat_copy


@task(returns=list)
# @task(adj_mat=INOUT)
def sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc, dist_mat,
                  quocient, res, *args):
    # TODO: change *args
    adj_mat_copy = [[] for _ in range(max(adj_mat[0], 1))]
    tmp_unwrap_1 = [i.value[1] for i in args]
    tmp_unwrap_1 = np.array([i for num, i in
                             enumerate(np.concatenate(tmp_unwrap_1))
                             if ((num % quocient) == res)])
    tmp_unwrap_2 = [neigh_sq_loc[i] for i in range(len(neigh_sq_loc))
                    for j in range(len(args[i].value[1]))]
    tmp_unwrap_2 = [i for num, i in enumerate(tmp_unwrap_2) if
                    ((num % quocient) == res)]
    nice_args = [tmp_unwrap_1, tmp_unwrap_2]
    for num in range(dist_mat.shape[0]):
        current_clust_id = int(data.value[1][num])
        if current_clust_id > -1:
            # for tmp2, loc2 in nice_args:
            for num2 in range(len(nice_args[0])):
                if dist_mat[num][num2*quocient + res]:
                    clust_ind = int(nice_args[0][num2])
                    adj_mat_elem = [nice_args[1][num2], clust_ind]
                    if ((clust_ind > -1) and
                        (adj_mat_elem not in
                            adj_mat_copy[current_clust_id])):
                        adj_mat_copy[current_clust_id].append(adj_mat_elem)
    return adj_mat_copy


def unwrap_adj_mat(tmp_mat, adj_mat, neigh_sq_coord, dimension_perms):
    links_list = []
    for comb in itertools.product(*dimension_perms):
        for k in range(tmp_mat[comb][0]):
            links_list.append([comb, k])
    mf_set = DisjointSet(links_list)
    for comb in itertools.product(*dimension_perms):
            # For an implementation issue, elements in adj_mat are
            # wrapped with an extra pair of [], hence the j[0]
        for k in range(len(adj_mat[comb][0])-1):
            mf_set.union(adj_mat[comb][0][k], adj_mat[comb][0][k+1])
    out = mf_set.get()
    return out


@task()
def expand_cluster(data, epsilon, border_points, dimension_perms, links_list,
                   square, cardinal, file_id, *args):
    data_copy = Data()
    data_copy.value = copy.deepcopy(data.value)
    tmp_unwrap_2 = [i.value[1] for i in args]
    neigh_points_clust = np.concatenate(tmp_unwrap_2)
    for elem in border_points:
        for p in border_points[elem]:
            if neigh_points_clust[p] > -1:
                data_copy.value[1][elem] = neigh_points_clust[p]
                break
    for num, elem in enumerate(data_copy.value[1]):
        if elem < -2:
            clust_ind = -1*elem - 3
            data_copy.value[1][num] = neigh_points_clust[clust_ind]

    # Map all cluster labels to general numbering
    mappings = []
    for k, links in enumerate(links_list):
        for pair in links:
            if pair[0] == square:
                mappings.append([pair[1], k])
    for num, elem in enumerate(data_copy.value[1]):
        if elem > -1:
            for pair in mappings:
                if int(elem) == pair[0]:
                    data_copy.value[1][num] = pair[1]

    # Update all files (for the moment writing to another one)
    # path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_id)
    path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(square[0])
    for num, j in enumerate(square):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += "_OUT.txt"
    f_out = open(tmp_string, "w")
    for num, val in enumerate(data_copy.value[0]):
        f_out.write(str(data_copy.value[0][num])+" "
                    + str(int(data_copy.value[1][num])) + "\n")
    f_out.close()


def DBSCAN(epsilon, min_points, file_id):
    #   TODO: code from scratch the Disjoint Set
    #   TODO: comment the code apropriately

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
    dist_mat = defaultdict()
    neigh_sq_coord = defaultdict()
    len_datasets = defaultdict()
    for comb in itertools.product(*dimension_perms):
        # TODO: implement init_data as a class method maybe inside
        # the initialisation
        dataset[comb] = Data()
        dataset_tmp[comb] = Data()
        dist_mat[comb] = Data()
        len_datasets[comb] = count_lines(comb, file_id)
        init_data(dataset_tmp[comb], comb, file_id)
        neigh_sq_coord[comb] = neigh_squares_query(comb, epsilon,
                                                   dimensions)

    # Partial Scan And Initial Cluster merging
#    dataset_tmp = dict_compss_wait_on(dataset_tmp, dimension_perms)
    adj_mat = defaultdict()
    tmp_mat = defaultdict()
    border_points = defaultdict()
    for comb in itertools.product(*dimension_perms):
        adj_mat[comb] = [0]
        tmp_mat[comb] = [0]
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset_tmp[coord])
        # Adj_mat as inout?
        [dataset[comb],
         border_points[comb],
         tmp_mat[comb],
         adj_mat[comb],
         dist_mat[comb]] = partial_scan_merge(dataset_tmp[comb],
                                              epsilon, min_points,
                                              *neigh_squares)

    # Cluster Synchronisation
    for comb in itertools.product(*dimension_perms):
#        compss_delete_object(dataset_tmp[comb])
        neigh_squares = []
        neigh_squares_loc = []
        len_neighs = 0
        for coord in neigh_sq_coord[comb]:
            neigh_squares_loc.append(coord)
            neigh_squares.append(dataset[coord])
            len_neighs += len_datasets[coord]
        # TODO: make as INOUT instead of OUT, currently not working
        fut_list = orquestrate_sync_clusters(dataset[comb], adj_mat[comb],
                                             epsilon, comb, neigh_squares_loc,
                                             dist_mat[comb], len_neighs, 1, 0,
                                             [], *neigh_squares)
        adj_mat[comb] = merge_task(adj_mat[comb], *fut_list)

    # Cluster list update
    # TODO: join the three
    border_points = dict_compss_wait_on(border_points, dimension_perms)
    tmp_mat = dict_compss_wait_on(tmp_mat, dimension_perms)
    adj_mat = dict_compss_wait_on(adj_mat, dimension_perms)
    links_list = unwrap_adj_mat(tmp_mat, adj_mat, neigh_sq_coord,
                                dimension_perms)
    for comb in itertools.product(*dimension_perms):
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset[coord])
        expand_cluster(dataset[comb], epsilon, border_points[comb],
                       dimension_perms, links_list, comb, tmp_mat[comb],
                       file_id, *neigh_squares)
    print links_list
    print "Time elapsed: " + str(time.time()-initial_time)
    return 1


if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
