# DBSCAN 4 PyCOMPSs
# carlos.segarra @ bsc.es

# Imports
from pycompss.api.task import task
from pycompss.api.parameter import INOUT
from pycompss.api.api import compss_wait_on
from collections import defaultdict
from classes.DS import DisjointSet
from classes.Data import Data
import itertools
# import time
import numpy as np
import sys
# import math


@task(data_pos=INOUT)
def init_data(data_pos, square, num_points_max, means, std):
    for i in range(len(means)):
        for j in range(num_points_max):
            tmp = np.random.multivariate_normal(means[i], std[i], size=(1))
            for k in range(len(square)):
                if (tmp[0][k] <= float(square[k])/10 or
                        tmp[0][k] > float(square[k] + 1)/10):
                    break
            else:
                data_pos.value.append(tmp)
    if len(data_pos.value) > 0:
        data_pos.value = np.vstack(data_pos.value)
        tmp_vec = -2*np.ones(np.shape(data_pos.value)[0])
        data_pos.value = [data_pos.value, tmp_vec]
    else:
        data_pos.value = np.array([(np.random.uniform(low=float(square[0])/10,
                                    high=float(square[0] + 1)/10),
                                    np.random.uniform(low=float(square[1]) /
                                    10, high=float(square[1]+1)/10))])
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
            neigh_squares.append(current)
    return neigh_squares


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


@task(adj_mat=INOUT)
def sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc, *args):
    nice_args = [[args[i], neigh_sq_loc[i]] for i in range(len(neigh_sq_loc))]
    nice_args.append([data, coord])
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


def unwrap_adj_mat(tmp_mat, adj_mat):
    links_list = []
    for i in range(len(tmp_mat)):
        for j in range(len(tmp_mat[i])):
            for k in range(tmp_mat[i][j][0]):
                links_list.append([[i, j], k])
    mf_set = DisjointSet(links_list)
    for i in adj_mat:
        for j in i:
            # For an implementation issue, elements in adj_mat are
            # wrapped with an extra pair of [], hence the j[0]
            for k in range(len(j[0])-1):
                mf_set.union(j[0][k], j[0][k+1])
    out = mf_set.get()
    return out


@task(data=INOUT)
def expand_cluster(data, epsilon, border_points, *args):
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


def DBSCAN(epsilon, min_points):
    #   TODO: code from scratch the Disjoint Set
    #   TODO: comment the code apropriately
    #   TODO: remove all the hardcoded parts add a dim input argument

    # Initial Definitions
    epsilon = float(epsilon)
    min_points = int(min_points)
    # This variables are currently hardcoded
    num_grid_rows = 10
    num_grid_cols = 10
    num_points_max = 100

    # Data inisialitation
    centers = [[0.2, 0.3], [0.6, 0.7]]
    std = [[[0.01, 0], [0, 0.01]], [[0.01, 0], [0, 0.01]]]
    dataset = [[Data() for _ in range(num_grid_cols)] for __ in
               range(num_grid_rows)]
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            # This ONLY WORKS IN 2D
            init_data(dataset[i][j], [i, j], num_points_max,
                      centers, std)

    # Partial Scan And Initial Cluster merging
    adj_mat = [[[] for _ in range(num_grid_cols)] for __ in
               range(num_grid_rows)]
    border_points = [[[] for _ in range(num_grid_cols)] for __ in
                     range(num_grid_rows)]
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            neigh_sq_coord = neigh_squares_query([i, j], epsilon)
            neigh_squares = []
            for coord in neigh_sq_coord:
                neigh_squares.append(dataset[coord[0]][coord[1]])
            # Use border points and adj_mat as INOUT instead of OUT
            border_points[i][j] = partial_scan(dataset[i][j], epsilon,
                                               min_points, *neigh_squares)
            adj_mat[i][j] = merge_cluster(dataset[i][j],  epsilon)

    # Cluster Synchronisation
    adj_mat = compss_wait_on(adj_mat)
    border_points = compss_wait_on(border_points)
    import copy
    tmp_mat = copy.deepcopy(adj_mat)
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            adj_mat[i][j] = [[] for _ in range(max(adj_mat[i][j][0], 1))]
            neigh_sq_coord = neigh_squares_query([i, j], epsilon)
            neigh_squares = []
            neigh_squares_loc = []
            for coord in neigh_sq_coord:
                neigh_squares_loc.append([coord[0], coord[1]])
                neigh_squares.append(dataset[coord[0]][coord[1]])
            sync_clusters(dataset[i][j], adj_mat[i][j], epsilon, [i, j],
                          neigh_squares_loc,  *neigh_squares)

    # Cluster list update
    adj_mat = compss_wait_on(adj_mat)
    links_list = unwrap_adj_mat(tmp_mat, adj_mat)
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            neigh_sq_coord = neigh_squares_query([i, j], epsilon)
            neigh_squares = []
            for coord in neigh_sq_coord:
                neigh_squares.append(dataset[coord[0]][coord[1]])
            expand_cluster(dataset[i][j], epsilon, border_points[i][j],
                           *neigh_squares)
    print links_list
    return links_list


if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]))
