#Implement adaptative space partition
#carlos.segarra @ bsc.es

#Imports
#from pycompss.api.task import task
#from pycompss.api.parameter import *
#from pycompss.api.api import compss_wait_on
from collections import defaultdict
from classes.cluster import Cluster
from classes.DS import DisjointSet
import itertools
import time
import numpy as np
import sys
import math
import random
 
##This might change to include sigmoid normalisation, rather than dividing by
##the maximum value what might result in numerical inestability.
##Complete migration to NumPy?
#def normalize_data(data_file):
#    """
#    Given a dataset, divide each dimension by its maximum obtaining data 
#    in the range [0,1]^n
#    :param dataFile: name of the original data file.
#    :return newName: new name of the data file containing the normalized 
#    data.
#    """
#    dataset = np.loadtxt(data_file)
#    norm_data = np.where(np.max(dataset, axis=0) == 0, dataset, 
#        dataset*1./np.max(dataset, axis=0))
#    new_name = data_file[data_file.find('.') + 1:data_file.rfind('.')]
#    new_name = '.' + new_name + 'n.txt'
#    np.savetxt(new_name,norm_data)
#    return new_name

#TRASH function
def consume_data(dataset, filename):
    data_to_load = np.loadtxt(filename)
    for point in data_to_load:
        i = int(point[0]*10)
        j = int(point[1]*10)
        dataset[i][j].append(point)
    for i, row in enumerate(dataset):
        for j, square in enumerate(row):
            if len(square) > 0:
                dataset[i][j] = np.vstack(square)
            else:
                dataset[i][j] = np.array([(np.random.uniform(low = float(i)/10
                    , high = float(i+1)/10), np.random.uniform(low = 
                    float(j)/10, high = float(j+1)/10))])
    return dataset

#@task(data_pos = INOUT)
def init_data(data_pos, square, num_points):
    """
    FOR THE MOMENT IT ONLY WORKS IN A  10x10 MATRIX
    Initializes random generated data? in the given square.
    :inout data_pos: pointer to the square in the grid where points will be 
    stored.
    :param square: indications to the given square to generate data belonging
    only to this square. Square is an int and each number from right to left 
    denotes the cell number along that dimension (dim_n-1 dim_n-2 ... dim 0) 
    """
    dim=len(square)
    data_pos = np.random.sample([num_points, dim])
    for i in range(dim):
        data_pos[:,i] = data_pos[:,i]/10 + float(square[i])/10 
    #DELETE
    return data_pos

#Only works in 2D and 10 squares per side
def neigh_squares_query(square, epsilon):
    #This could be modified to include only borders. 0.1 is the side size
    dim = len(square)
    neigh_squares = []
    border_squares = int(max(epsilon/0.1, 1))
    perm = []
    for i in range(dim):
        perm.append(range(-border_squares, border_squares + 1))
    for comb in itertools.product(*perm):
        current=[]
        for i in range(dim):
            if square[i] + comb[i] in range(10):
                current.append(square[i]+comb[i])
        if len(current) == dim and current != square:
            neigh_squares.append(current)
    return neigh_squares    

#@task(core_points_list = INOUT)
def partial_scan(data, core_points_list, epsilon, min_points, *args):
    neigh_points = np.vstack(args)
    for point in data:
        if np.sum((np.linalg.norm(neigh_points - point, axis = 1) - 
            epsilon) < 0) > min_points:
            core_points_list.append(point)
    #DELETE
    return core_points_list

#@task(clusters = INOUT) 
def merge_cluster(clusters, core_points, epsilon):
    for point in core_points:
        for clust in clusters:
            if np.sum((np.linalg.norm(clust - point, axis = 1) - 
                epsilon) < 0) > 0:
                clust.append(point)
                break
        else:
            clusters.append([point])
    for num, clust in enumerate(clusters):
        np.vstack(clust)
    #DELETE
    return clusters

#@task(adj_mat = INOUT)
def sync_clusters(clusters, adj_mat, epsilon, *args):
    adj_mat = [[] for _ in range(len(clusters))]
    for num1, clust in enumerate(clusters):
        for neigh_with_coord in args:
            neigh_clusters = neigh_with_coord[0]
            neigh_coord = neigh_with_coord[1]
            for num2, neigh_clust in enumerate(neigh_clusters):
                for point in clust:
                    if np.sum((np.linalg.norm(neigh_clust - point, axis = 1) -
                        epsilon) < 0) > 0:
                        adj_mat[num1].append([neigh_coord[0], neigh_coord[1], 
                            num2])
                        break
    #DELETE
    return adj_mat

def unwrap_adj_mat(adj_mat):
    cardinals = [[len(adj_mat[i][j]) for j in range(len(adj_mat[i]))] for i in
        range(len(adj_mat))]
    cum_sum = np.cumsum(cardinals, axis = 1)
    for i in range(1, len(cum_sum) - 1):
        cum_sum[i] += cum_sum[i-1][-1]
    count = 0
    links = defaultdict(list)
    for i in range(len(adj_mat)):
        for j in range(len(adj_mat[i])):
            for k in range(len(adj_mat[i][j])):
                links[count].append(count)
                for num, l in enumerate(adj_mat[i][j][k]):
                    clust_id = cum_sum[i][j] + num
                    links[count].append(clust_id)
                count += 1
    return links

def update(links_list):
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
            
#@task(clustPoints = INOUT)
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
    epsilon = float(epsilon)
    min_points = int(min_points)
    start = time.time()
    #This variables are currently hardcoded
    num_grid_rows = 10
    num_grid_cols = 10
    num_points = 5
    dim = 2
    dataset = [[[] for _ in range(num_grid_cols)] for __ in
        range(num_grid_rows)]
#    dataset = np.empty([num_grid_rows, num_grid_cols, num_points, dim])
#    for i in range(num_grid_rows):
#        for j in range(num_grid_cols):
            #This ONLY WORKS IN 2D
            #init_data(dataset[i][j], [i,j], num_points)  
            #DELETE
#            dataset[i][j] = init_data(dataset[i][j], [i,j], num_points)  
    dataset = consume_data(dataset, "./data/blobs_0to1.txt")
    clusters = [[[] for _ in range(num_grid_cols)] for __ in 
        range(num_grid_rows)]
    core_points = [[[] for _ in range(num_grid_cols)] for __ in 
        range(num_grid_rows)]
    adj_mat = [[[] for _ in range(num_grid_cols)] for __ in 
        range(num_grid_rows)]
    for i in range(num_grid_rows):
        for j in range(num_grid_cols):
            neigh_sq_coord = neigh_squares_query([i,j], epsilon)
            neigh_squares = []
            for coord in neigh_sq_coord:
                neigh_squares.append(dataset[coord[0]][coord[1]])
#            partial_scan(dataset[i][j], core_points[i][j], epsilon, min_points,
##                 *neigh_squares) 
##            merge_cluster(clusters[i][j], core_points[i][j], epsilon)
            #DELETE
            core_points[i][j] = partial_scan(dataset[i][j], core_points[i][j], epsilon, min_points,
                 *neigh_squares) 
            clusters[i][j] = merge_cluster(clusters[i][j], core_points[i][j], epsilon)
            neigh_clusters = []
            for coord in neigh_sq_coord:
                neigh_clusters.append([clusters[coord[0]][coord[1]], coord])
            adj_mat[i][j] = sync_clusters(clusters[i][j], adj_mat[i][j], epsilon, 
                *neigh_clusters) 
##    adj_mat = compss_wait_on(adj_mat)
    link_list = unwrap_adj_mat(adj_mat)
    link_list = update(link_list)

    #What to do with the output?
    return link_list
