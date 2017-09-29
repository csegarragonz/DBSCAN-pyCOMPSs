#Implement adaptative space partition
#carlos.segarra @ bsc.es

#Imports
from pycompss.api.task import task
from pycompss.api.parameter import *
from pycompss.api.api import compss_wait_on
from collections import defaultdict
from classes.cluster import Cluster
from classes.DS import DisjointSet
import itertools
import time
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
    perm = []
    for i in range(dim):
        perm.append(range(-k[i],k[i] + 1))
    for comb in itertools.product(*perm):
        current = square
        for i in range(dim):
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
        else:
            #Class Cluster is defined at /classes/cluster.py
            tmp = ([point],[square])
            tmpc = Cluster()
            tmpc.add(*tmp)
            clusters.append(tmpc)

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
                break

def DBSCAN(dataFile, fragSize, epsilon, minPoints, numParts):
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
    normData = normalize_data(dataFile) 
    print "Normalize data: " + str(normData)
    dataset = np.loadtxt(normData)
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
