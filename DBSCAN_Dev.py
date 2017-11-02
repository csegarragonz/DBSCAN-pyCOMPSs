# Implement adaptative space partition
# carlos.segarra @ bsc.es

# Imports
from pycompss.api.task import task
from pycompss.api.parameter import INOUT
from pycompss.api.api import compss_wait_on
# from pycompss.api.api import barrier
from collections import defaultdict
# from classes.cluster import Cluster
from classes.DS import DisjointSet
from classes.Data import Data
import itertools
# import time
import numpy as np
import sys
import math
import random

#This might change to include sigmoid normalisation, rather than dividing by
#the maximum value what might result in numerical inestability.
#Complete migration to NumPy?
def normalize_data(data_file):
    """
    Given a dataset, divide each dimension by its maximum obtaining data
    in the range [0,1]^n
    :param dataFile: name of the original data file.
    :return newName: new name of the data file containing the normalized
    data.
    """
    dataset = np.loadtxt(data_file)
    norm_data = np.where(np.max(dataset, axis=0) == 0, dataset,
        dataset*1./np.max(dataset, axis=0))
    new_name = data_file[data_file.find('.') + 1:data_file.rfind('.')]
    new_name = '.' + new_name + 'n.txt'
    np.savetxt(new_name,norm_data)
    return new_name

def partition_space(dataset, fragSize, epsilon):
    """
    Returns:
    :returns fragData: dict with (key,value)=(spatial id, points in the
    region).
    :returns rangeToEps: list with the number of squares to check to
    preserve eps distance to all points along each dim.
    fragVec begins at 0
    """
    fragData = defaultdict(list)
    rangeToEps = defaultdict(list)
    dim = len(dataset[0])
    fragVec = [[np.max(np.min(dataset, axis=0)[i] - epsilon,0),
        np.mean(dataset, axis = 0)[i], np.max(dataset, axis=0)[i]]
        for i in range(dim)]
    size = pow(10, len(str(fragSize + 1)))
    for i in range(dim):
        for j in range(fragSize):
            tmpPoint = defaultdict(int)
            for point in dataset:
                k = 0
                while point[i]>fragVec[i][k]: k += 1
                tmpPoint[k] += 1
            ind = max(tmpPoint.iterkeys(), key=(lambda key:
                tmpPoint[key]))
            val = float((fragVec[i][ind-1] + fragVec[i][ind])/2)
            fragVec[i].insert(ind, val)
    for point in dataset:
        key = 0
        for i in range(dim):
            k = 0
            while point[i]>fragVec[i][k]: k += 1
            key += (k-1)*pow(size,i)
        fragData[key].append(point)
    for square in fragData:
        tmp = []
        for i in range(dim):
            pos = square % size
            a = [[j,x-fragVec[i][pos]] for j,x in enumerate(fragVec[i])
                if abs(x - fragVec[i][pos]) < epsilon]
            b = [[j,x-fragVec[i][pos+1]] for j,x in enumerate(fragVec[i])
                if abs(x - fragVec[i][pos+1]) < epsilon]
            maxa = abs(max(a, key=lambda x: x[1])[0] - pos)
            maxb = abs(max(b, key=lambda x: x[1])[0] - pos)
            tmp.append(max(maxa, maxb, 1))
            pos = pos/size
        rangeToEps[square] = tmp
    return (fragData, fragVec, rangeToEps)

def partial_scan(corePoints, square, epsilon, minPoints, fragData,
                fragSize, numParts, rangeToEps):
    """
    Looks for all the core points (over minPoints neighbours) inside a
    certain square.
    :param square: space region where core points are looked for.
    :param epsilon: maximum distance between two points to be considered
    neighbours.
    :param minPoints: minimum number of neighbours for a point to be
    considered core point.
    :param fragData: dict containing space partition and its points.
    :param fragSize: length of the side of each hypercube inside fragData.
    :return:    list of all the core points found in square.
    """
    dim = len(fragData[square][0])
    pointSet = fragData[square]
    pointSetReal = pointSet[:]
    k=rangeToEps[square]
    size = pow(10,len(str(fragSize+1)))
=======
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
>>>>>>> labeling_clusters
    perm = []
    for i in range(dim):
        perm.append(range(-border_squares, border_squares + 1))
    for comb in itertools.product(*perm):
        current = []
        for i in range(dim):
<<<<<<< HEAD
            current = current + comb[i]*math.pow(size,i)
        if current in fragData and current != square:
            pointSet = pointSet + fragData[current]
    for i in range(numParts):
        tmpPS = [p for j,p in enumerate(pointSetReal) if j % numParts ==
            i]
        scan_task(corePoints[i], pointSet, tmpPS, epsilon, minPoints)

@task(corePoints = INOUT)
def scan_task(corePoints, pointSet, pointSetReal, epsilon, minPoints):
    for point in pointSetReal:
        neighbourPts = 0
        for p in pointSet:
            if np.linalg.norm(point - p) < epsilon:
                neighbourPts = neighbourPts+1
                if neighbourPts >= minPoints:
                    corePoints.append(point)
                    break

@task(clusters = INOUT)
def merge_cluster(clusters, corePoints, square, epsilon):
    """
    Given a list of core points, decides wether each core point forms a
    new cluster or belongs to an existing one.
    :inout clusters:    list of clusters found.
    :param corePoints:   list of all the core points found in square.
    :param square:     working square.
    :param epsilon:    input parameter.
    """
    for point in corePoints:
        possibleClusters = []
        for clust in clusters:
            for clustP in clust.points:
                if np.linalg.norm(point-clustP) < epsilon:
                    possibleClusters.append(clust)
                    if len(possibleClusters) > 1:
                        clusters.remove(clust)
                    break
        if len(possibleClusters) > 0:
            master = possibleClusters.pop(0)
            master.addPoint(point)
            for slave in possibleClusters:
                master.merge(slave)
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

def sync_clusters(clusters, epsilon, numParts):
    """
    Returns a matrix of booleans. Pos [i,j]=1 <=> clusters -i and -j
    should be merged.
    :inout clusters: list of all clusters and their points.
    :param numParts: list of clusters to be merged with the general list.
    :param epsilon: input parameter, maximum distance under which two
    points are considered neighbours.
    :return possibleClusters: result matrix.
    """
    possibleClusters = defaultdict(list)
    n = 10
    tmpDel = (n - len(clusters) % n) % n
    if not len(clusters) % n ==0:
        for i in range(n - len(clusters) % n): clusters.append(Cluster())
    enable = [[[] for z in range(len(clusters)/n)] for w in
        range(len(clusters)/n)]
    for i in range(len(clusters)/n):
        tmp1 = [clusters[n*i + k].points for k in range(n)]
        for j in range(len(clusters)/n):
            tmp2=[clusters[n*j + k].points for k in range(n)]
            sync_task(enable[i][j], tmp1, tmp2, epsilon)
    a = time.time()
    enable = compss_wait_on(enable)
    print "Long Wait on Time: %.6f" % float(time.time()-a)
    for i in range(len(clusters)/n):
        for j in range(len(clusters)/n):
            for k in range(n):
                for l in range(n):
                    if enable[i][j][n*k +l]:
                        possibleClusters[i*n+k].append(j*n+l)
    l = len(clusters)
    for i in range(tmpDel):
        possibleClusters.pop(l-1-i, None)
        clusters.pop()
    return possibleClusters

@task(enable = INOUT)
def sync_task(enable, hosts, visits, epsilon):
    """
    Given two lists, checks wether the distance between the two sets is
    less than epsilon.
    """
    for host in hosts:
        for visit in visits:
            for p in host:
                for point in visit:
                    if np.linalg.norm(point-p) < epsilon:
                        enable.append(1)
                        break
                else: continue
                break
            else: enable.append(0)

def update(clusters, possibleClusters, returnCluster):
    """
    Given a list of clusters, it reassigns each point to a cluster making
    sure a point does not appear in more than one different clusters. If
    so, those clusters will be merged.
    :param clusters: provisional cluster list result of the cluster
    expansion.
    :return defClusters: final cluster list.
    """
    MF_set = DisjointSet(range(len(clusters)))
    for i in possibleClusters:
        for j in range(len(possibleClusters[i])-1):
            MF_set.union(possibleClusters[i][j], possibleClusters[i][j+1])
    a = MF_set.get()
    if returnCluster:
        defCluster = [Cluster() for _ in range(len(a))]
        for i,lst in enumerate(a):
            for elem in lst:
                defCluster[i].merge(clusters[elem])
        return defCluster
    defCluster = [[] for _ in range(len(a))]
    for i,lst in enumerate(a):
        for elem in lst:
            for point in clusters[elem].points:
                defCluster[i].append(point)
    return defCluster

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

@task(clustPoints = INOUT)
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
    print len(links_list)
    return links_list


def DBSCAN(data_file, fragSize, epsilon, minPoints, numParts):
    """
    Main DBSCAN algorithm.
    :param dataFile: file where the data is stored.
    :param fragSize: size used for space partition.
    :param epsilon: maximum distance for two points to be considered
    neighbours.
    :param minPoints: minimum number of neighbours for a point to be
    considered core point.
    :param numParts: number of parts in which fragData is divided for
    processing.
    :return defClusters: result of the application, dict where values
    are points corresponding to the <key> cluster.
    """
    start = time.time()
    print "Density Based Scan started."
    norm_data = normalize_data(data_file)
    dataset = np.loadtxt(norm_data)
    [fragData, fragVec, rangeToEps] = partition_space(dataset, fragSize,
        epsilon)
    print "Starting partial scan..."
    clusters = [[[] for _ in range(len(fragData))] for __ in
        range(numParts)]
    corePMatrix = [[[] for _ in range(numParts)] for __ in
        range(len(fragData))]
    for i, (square,value) in enumerate(fragData.viewitems()):
        partial_scan(corePMatrix[i], square, epsilon, minPoints, fragData,
            fragSize, numParts, rangeToEps)
        for j in range(numParts):
            merge_cluster(clusters[j][i], corePMatrix[i][j], square,
                epsilon)
    clusters = compss_wait_on(clusters)
    print "Initial Proposal Finished"
    iniP = time.time()
    clusters = [clust for rows in clusters for squares in rows for clust
        in squares]
    print "Length of clusters found: "+str(len(clusters))
    possibleClusters = sync_clusters(clusters, epsilon, numParts)
    syncT = time.time()
    print "Syncing Finished"
    halfDefClusters = update(clusters, possibleClusters, True)
    updateTime = time.time()
    defClusters = expand_cluster(halfDefClusters, fragData, epsilon,
        minPoints, fragSize, numParts, rangeToEps)
    print "Expansion Finished"
    print "DBSCAN Algorithm finished succesfully."
    end = time.time()
    print "Exec Times:"
    print "----------------------------------"
    print "Initial Proposal Time: \t %.6f" % float(iniP-start)
    print "Syncing Time: \t \t  %.6f" % float(syncT-iniP)
    print "Update Time: \t \t  %.6f" % float(updateTime-syncT)
    print "Expand: \t \t  %.6f" % float(end-updateTime)
    print "----------------------------------"
    print "Time elapsed: \t \t  %.6f" % float(end-start)
    print "Number of clusters found: "+str(len(defClusters))
    sys.stdout.flush()
    sys.stderr.flush()
    return [defClusters, fragVec]

if __name__ == "__main__":
    DBSCAN(float(sys.argv[1]), int(sys.argv[2]))
