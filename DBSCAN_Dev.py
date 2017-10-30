# Implement adaptative space partition
# carlos.segarra @ bsc.es

# Imports
from pycompss.api.task import task
# from pycompss.api.parameter import *
from pycompss.api.parameter import INOUT
from pycompss.api.api import compss_wait_on
# from pycompss.api.api import barrier
from collections import defaultdict
# from classes.cluster import Cluster
from classes.DS import DisjointSet
from classes.Data import Data
import itertools
import time
import numpy as np
import sys
import math


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


@task(data=INOUT)
def partial_scan(data, core_points_list, epsilon, min_points, *args):
    tmp_unwrap = [i.value[0] for i in args]
    neigh_points = np.vstack(tmp_unwrap)
    for num, point in enumerate(data.value[0]):
        if np.sum((np.linalg.norm(neigh_points - point, axis=1) -
                   epsilon) < 0) > min_points:
            core_points_list.append(point)
            data.value[1][num] = -1
#    return data


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
#    adj_mat = [[] for _ in range(max(adj_mat[0], 1))]
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
                            print adj_mat
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


def update(links_list, adj_mat):
    mf_set = DisjointSet(range(len(links_list)))
    for i in links_list:
        for j in range(len(links_list[i])-1):
            mf_set.union(links_list[i][j], links_list[i][j+1])
    final_mf_set = mf_set.get()
    return final_mf_set


def expand_cluster(clusters, fragData, epsilon, minPoints, fragSize,
                    numParts, rangeToEps):
    """
    Expands all clusters contained in a list of clusters. Expansion means
    adding non-core points to the
    already established clusters.
    """
    addToClust = [[[] for _ in range(numParts)] for x in
        range(len(clusters))]
    for numClust,clust in enumerate(clusters):
        for k in range(numParts):
            neigh_expansion(addToClust[numClust][k], clust, fragData,
                fragSize, rangeToEps, epsilon)
    addToClust = compss_wait_on(addToClust)
    for i,m in enumerate(addToClust):
        addToClust[i] = [j for k in m for j in k]
    pointsToClust = defaultdict(list)
    links = defaultdict(list)
    for i,clust in enumerate(addToClust):
        for p in addToClust[i]:
            if str(p) in pointsToClust:
                for c in pointsToClust[str(p)]:
                    if not i in links[c]: links[c].append(i)
                    if not c in links[i]: links[i].append(c)
                pointsToClust[str(p)].append(i)
            else:
                pointsToClust[str(p)].append(i)
                clusters[i].addPoint(p)
    return update(clusters, links, False)

##@task(clustPoints = INOUT)


def neigh_expansion(clustPoints, clust, fragData, fragSize, rangeToEps,
                    epsilon):
    """
    Given a cluster of core points, returns all the points lying in the
    squares reachable with epsilon distance
    from the clusters square.
    :param cluster: cluster whose neighbour points will be added.
    :param fragData: dict containing space partitioning.
    :param epsilon: maximum distance between two points to be considered
    neighbours.
    :param minPoints: minimum number of neighbours for a point to be
    considered core point.
    :param fragSize: length of the side of each hypercube inside fragData.
    """
    squaresNot = clust.square[:]
    pointSet = []
    for sq in clust.square:
        #Segur que no ho sumem dos cops aleshores?
        pointSet += fragData[sq][:]
        dim = len(fragData[sq][0])
        k = rangeToEps[sq]
        size = pow(10,len(str(fragSize+1)))
        perm = []
        for i in range(dim):
            perm.append(range(-k[i],k[i]+1))
        for comb in itertools.product(*perm):
            current = sq
            for i in range(dim):
                current = current+comb[i]*math.pow(size,i)
            if current in fragData and not(current in squaresNot):
                pointSet = pointSet+fragData[current][:]
                squaresNot.append(current)
    tmp = [b for b in pointSet if not(b in clust.points)]
    for i,point in enumerate(tmp):
        for p in clust.points:
            if np.linalg.norm(point-p) < epsilon:
                clustPoints.append(point)
                break


def DBSCAN(epsilon, min_points):
    #   TODO: code from scratch the Disjoint Set
    epsilon = float(epsilon)
    min_points = int(min_points)
#    start = time.time()
    # This variables are currently hardcoded
    num_grid_rows = 10
    num_grid_cols = 10
    num_points_max = 100
    centers = [[0.2, 0.3], [0.6, 0.7]]
    std = [[[0.01, 0], [0, 0.01]], [[0.01, 0], [0, 0.01]]]
    dataset = [[Data() for _ in range(num_grid_cols)] for __ in
               range(num_grid_rows)]
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            # This ONLY WORKS IN 2D
            init_data(dataset[i][j], [i, j], num_points_max,
                      centers, std)
    core_points = [[[] for _ in range(num_grid_cols)] for __ in
                   range(num_grid_rows)]
    adj_mat = [[[] for _ in range(num_grid_cols)] for __ in
               range(num_grid_rows)]
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            neigh_sq_coord = neigh_squares_query([i, j], epsilon)
            neigh_squares = []
            for coord in neigh_sq_coord:
                neigh_squares.append(dataset[coord[0]][coord[1]])
            partial_scan(dataset[i][j], core_points[i][j], epsilon, min_points,
                         *neigh_squares)
            adj_mat[i][j] = merge_cluster(dataset[i][j],  epsilon)
    adj_mat = compss_wait_on(adj_mat)
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
    adj_mat = compss_wait_on(adj_mat)
    links_list = unwrap_adj_mat(tmp_mat, adj_mat)
#    links_list = update(links_list, adj_mat)
    print links_list
    return links_list


if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]))
