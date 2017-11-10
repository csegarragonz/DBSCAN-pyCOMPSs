# DBSCAN 4 PyCOMPSs
# carlos.segarra @ bsc.es

# Imports
from pycompss.api.task import task
from pycompss.api.parameter import INOUT
from pycompss.api.api import compss_wait_on
from collections import defaultdict
from classes.DS import DisjointSet
from classes.Data import Data
from ast import literal_eval
import itertools
# import time
import numpy as np
import sys
import os


@task(data_pos=INOUT)
def init_data(data_pos, tupla, file_id):
    # path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_count)
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


def neigh_squares_query(square, epsilon):
    # This could be modified to include only borders. 0.1 is the side size
    dim = len(square)
    neigh_squares = []
    border_squares = int(max(epsilon/0.1, 1))
    perm = []
    for i in range(dim):
        perm.append(range(-border_squares, border_squares + 1))
    for comb in itertools.product(*perm):
        current = []
        for i in range(dim):
            if square[i] + comb[i] in range(10):
                current.append(square[i]+comb[i])
        if len(current) == dim and current != square:
            neigh_squares.append(tuple(current))
    neigh_squares.append(tuple(square))
    return tuple(neigh_squares)


@task(data=INOUT, returns=defaultdict)
def partial_scan(data, epsilon, min_points, *args):
    tmp_unwrap = [i.value[0] for i in args]
    tmp_unwrap_2 = [i.value[1] for i in args]
    neigh_points = np.vstack(tmp_unwrap)
    neigh_points_clust = np.concatenate(tmp_unwrap_2)
    non_assigned = defaultdict(int)
    for num, point in enumerate(data.value[0]):
        poss_neigh = np.linalg.norm(neigh_points - point, axis=1) - epsilon < 0
        neigh_count = np.sum(poss_neigh)
        if neigh_count > min_points:
            data.value[1][num] = -1
        elif neigh_count > 0:
            tmp = []
            for pos, proxim in enumerate(poss_neigh):
                if proxim:
                    # Adding three since later to detect core points we will
                    # require value > -1 and -0 is > 1 and not to confuse it
                    # with a noise point
                    if neigh_points_clust[pos] == -1:
                        data.value[1][num] = -(pos+3)
                        break
                    else:
                        tmp.append(pos)
            else:
                non_assigned[num] = tmp
    return non_assigned
#    return data, non_assigned


@task(data=INOUT, returns=int)
def merge_cluster(data, epsilon):
    cluster_count = 0
    core_points = [[num, p] for num, p in enumerate(data.value[0]) if
                   data.value[1][num] == -1]
    for pos, point in core_points:
        tmp_vec = (np.linalg.norm(data.value[0]-point, axis=1) - epsilon) < 0
        for num, poss_neigh in enumerate(tmp_vec):
            if poss_neigh and data.value[1][num] > -1:
                data.value[1][pos] = data.value[1][num]
                break
        else:
            data.value[1][pos] = cluster_count
            cluster_count += 1
    return [cluster_count]
#    return data, [cluster_count]


def dict_compss_wait_on(dicc, dimension_perms):
    for comb in itertools.product(*dimension_perms):
        dicc[comb] = compss_wait_on(dicc[comb])
    return dicc


@task(adj_mat=INOUT)
def sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc, *args):
    nice_args = [[args[i], neigh_sq_loc[i]] for i in range(len(neigh_sq_loc))]
    for num1, point1 in enumerate(data.value[0]):
        current_clust_id = int(data.value[1][num1])
        if current_clust_id > -1:
            for tmp2, loc2 in nice_args:
                tmp_vec = (np.linalg.norm(tmp2.value[0]-point1, axis=1) -
                           epsilon) < 0
                for num2, poss_neigh in enumerate(tmp_vec):
                    if poss_neigh:
                        clust_ind = int(tmp2.value[1][num2])
                        adj_mat_elem = [loc2, clust_ind]
                        if (clust_ind > -1) and (adj_mat_elem not in
                                                 adj_mat[current_clust_id]):
                            adj_mat[current_clust_id].append(adj_mat_elem)
#    return adj_mat


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


@task(data=INOUT)
def expand_cluster(data, epsilon, border_points, dimension_perms, *args):
    tmp_unwrap_2 = [i.value[1] for i in args]
    neigh_points_clust = np.concatenate(tmp_unwrap_2)
    for elem in border_points:
        for p in border_points[elem]:
            if neigh_points_clust[p] > -1:
                data.value[1][elem] = neigh_points_clust[p]
                break
    for num, elem in enumerate(data.value[1]):
        if elem < -2:
            clust_ind = -1*elem - 3
            data.value[1][num] = neigh_points_clust[clust_ind]
#    return data


def DBSCAN(epsilon, min_points, file_id):
    #   TODO: code from scratch the Disjoint Set
    #   TODO: comment the code apropriately
    #   TODO: remove hardcoded 0.1 side length
    #   TODO: use numpy masks.

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
    neigh_sq_coord = defaultdict()
    for comb in itertools.product(*dimension_perms):
        dataset[comb] = Data()
        init_data(dataset[comb], comb, file_id)
        neigh_sq_coord[comb] = neigh_squares_query(comb, epsilon)

    # Partial Scan And Initial Cluster merging
    adj_mat = defaultdict()
    border_points = defaultdict()
    for comb in itertools.product(*dimension_perms):
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset[coord])
        # TODO:Use border points and adj_mat as INOUT instead of OUT
        border_points[comb] = partial_scan(dataset[comb], epsilon,
                                           min_points, *neigh_squares)
        adj_mat[comb] = merge_cluster(dataset[comb],  epsilon)

    # Cluster Synchronisation
    adj_mat = dict_compss_wait_on(adj_mat, dimension_perms)
    border_points = dict_compss_wait_on(border_points, dimension_perms)
    import copy
    tmp_mat = copy.deepcopy(adj_mat)
    for comb in itertools.product(*dimension_perms):
        adj_mat[comb] = [[] for _ in range(max(adj_mat[comb][0], 1))]
        neigh_squares = []
        neigh_squares_loc = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares_loc.append(coord)
            neigh_squares.append(dataset[coord])
        sync_clusters(dataset[comb], adj_mat[comb], epsilon, comb,
                      neigh_squares_loc,  *neigh_squares)

    # Cluster list update
    adj_mat = dict_compss_wait_on(adj_mat, dimension_perms)
    links_list = unwrap_adj_mat(tmp_mat, adj_mat, neigh_sq_coord,
                                dimension_perms)
    for comb in itertools.product(*dimension_perms):
        neigh_squares = []
        for coord in neigh_sq_coord[comb]:
            neigh_squares.append(dataset[coord])
        expand_cluster(dataset[comb], epsilon, border_points[comb],
                       dimension_perms, *neigh_squares)
    print links_list
    return links_list


if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
