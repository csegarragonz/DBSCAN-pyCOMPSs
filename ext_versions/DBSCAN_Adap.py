#DBSCAN Tasques curtes. Carlos Segarra. THIS IS A PRELIMINARY VERSION
#Implement adaptative space partition
#carlos.segarra @ bsc.es

#Imports
from pycompss.api.task import task
from pycompss.api.parameter import *
from pycompss.api.api import compss_wait_on
from collections import defaultdict
from classes.cluster import Cluster
import itertools
import time
import numpy as np
import sys
import math
import random
 
def normalizeData(dataFile):
    """
    Given a dataset, divide each dimension by its maximum obtaining data in the range [0,1]^n
    :param dataFile: name of the original data file.
    :return newName: new name of the data file containing the normalized data.
    """
    dataset=np.loadtxt(dataFile)
    normData=np.where(np.max(dataset, axis=0)==0, dataset, dataset*1./np.max(dataset, axis=0))
    newName=dataFile[dataFile.find('.')+1:dataFile.rfind('.')]
    newName='.'+newName+'n.txt'
    np.savetxt(newName,normData)
    return newName

def chunksByNumP(list, numParts):
    """
    Given a list, returns the same list chunked in numParts parts.
    :param list: list to be chunked.
    :param numParts: number of parts for the list to be divided.
    :return: list of lists.
    """
    return [[k for i,k in enumerate(list) if i%numParts==j] for j in range(numParts)]

def partitionSpace(dataset, fragSize, epsilon):
    """
    Returns: 
    :returns fragData: dict with (key,value)=(spatial id, points in the region).
    :returns rangeToEps: list with the number of squares to check to preserve eps distance to all points along each dim.
    fragVec begins at 0
    """
    fragData=defaultdict(list)
    rangeToEps=defaultdict(list)
    dim = len(dataset[0])
    fragVec=[[np.max(np.min(dataset, axis=0)[i]-epsilon,0),np.mean(dataset, axis=0)[i],np.max(dataset, axis=0)[i]] for i in range(dim)]
    size = pow(10,len(str(fragSize+1)))
    print "Size: "+str(size)
    for i in range(dim):
        for j in range(fragSize):
            tmpPoint = defaultdict(int)
            for point in dataset:
                k = 0
                while point[i]>fragVec[i][k]: k += 1
                tmpPoint[k] += 1        
            ind = max(tmpPoint.iterkeys(), key=(lambda key: tmpPoint[key]))
            val=float((fragVec[i][ind-1]+fragVec[i][ind])/2)
            fragVec[i].insert(ind, val)
    for point in dataset:
        key=0
        for i in range(dim):
            k=0
            while point[i]>fragVec[i][k]: k += 1
            key += (k-1)*pow(size,i)        
        fragData[key].append(point)
    for square in fragData:
        tmp=[]
        for i in range(dim):
            pos=square % size
            a=[[j,x-fragVec[i][pos]] for j,x in enumerate(fragVec[i]) if abs(x - fragVec[i][pos]) < epsilon]
            b=[[j,x-fragVec[i][pos+1]] for j,x in enumerate(fragVec[i]) if abs(x - fragVec[i][pos+1]) < epsilon]
            maxa=abs(max(a, key=lambda x: x[1])[0]-pos)
            maxb=abs(max(b, key=lambda x: x[1])[0]-pos)
            tmp.append(max(maxa, maxb, 1)) 
            pos = pos/size
        rangeToEps[square]=tmp
    print "Frag vec: "+str(fragVec)
    return (fragData, fragVec, rangeToEps)

@task(clusters = INOUT) 
def mergeCluster(clusters, corePoints, square, epsilon):
    """
    Given a list of core points, decides wether each core point forms a new cluster or belongs to an existing one.
    :inout clusters:    list of clusters found.
    :param corePoints:   list of all the core points found in square.
    :param square:     working square.
    :param epsilon:    input parameter.
    """
    for point in corePoints:
        possibleClusters=[]
        for clust in clusters:
            for clustP in clust.points:
                if np.linalg.norm(point-clustP) < epsilon:
                    possibleClusters.append(clust)
                    if len(possibleClusters) > 1:
                        clusters.remove(clust)
                    break
        if len(possibleClusters) > 0:
            master=possibleClusters.pop(0)
            master.addPoint(point)
            for slave in possibleClusters:
                master.merge(slave)
        else:
            #Class Cluster is defined at /classes/cluster.py
            tmp=([point],[square])
            tmpc=Cluster()
            tmpc.add(*tmp)
            clusters.append(tmpc)

def DBSCAN(dataFile, fragSize, epsilon, minPoints):
    """
    Main DBSCAN algorithm.
    :param dataset: normalized dataset. Is here because we may supress the main().
    :param epsilon: maximum distance for two points to be considered neighbours.
    :param minPoints: minimum number of neighbours for a point to be considered core point.
    :param fragData: dict with the partitioned space and the points classified.
    :param fragSize: size used for space partition.
    :param numParts: number of parts in which fragData is divided for processing.
    :return defClusters: result of the application, dict where values are points corresponding to the <key> cluster.
    """
    start=time.time()
    print "Density Based Scan started."
    #numParts is used to avoid having too many tasks. numParts=numProcessesGiven*k
    numParts=30
    normData=normalizeData(dataFile) 
    dataset=np.loadtxt(normData)
    [fragData, fragVec, rangeToEps]=partitionSpace(dataset, fragSize, epsilon) 
    print "Starting partial scan..."
    clusters=[[] for x in range(numParts)]
    corePMatrix=[[[[] for x in range(numParts)] for y in range(len(fragData))] for z in range(numParts)]
    for j in range(numParts):
        fragDataSub={key: value for i, (key, value) in enumerate(fragData.viewitems()) if i % numParts == j}
        for k,(square,value) in enumerate(fragDataSub.viewitems()):
            partialScan(corePMatrix[j][k], square, epsilon, minPoints, fragData, fragSize, numParts, rangeToEps)
            for n in range(numParts):
                mergeCluster(clusters[j], corePMatrix[j][k][n], square, epsilon)
    clusters=compss_wait_on(clusters)
    print "Initial Proposal Time: %.6f" % float(time.time()-start)
    clusters=[j for i in clusters for j in i] 
    np.set_printoptions(threshold=sys.maxint)
    print "Length of clusters found: "+str(len(clusters))
    expandCluster(clusters, fragData, epsilon, minPoints, fragSize, numParts, rangeToEps)
    print "Expansion Time: %.6f" % float(time.time()-start)
    possibleClusters=syncClusters(clusters, epsilon, numParts)
    print "Syncing Time: %.6f" % float(time.time()-start)
    print possibleClusters
    defClusters=update(clusters, possibleClusters)
    print "DBSCAN Algorithm finished succesfully."
    end = time.time()
    print "Time elapsed: %.6f" % float(end-start)
    sys.stdout.flush()
    sys.stderr.flush()
    return [defClusters, fragVec] 

def syncClusters(clusters, epsilon, numParts):
    """
    Returns a matrix of booleans. Pos [i,j]=1 <=> clusters -i and -j should be merged.
    :inout clusters: list of all clusters and their points.
    :param numParts: list of clusters to be merged with the general list.
    :param epsilon: input parameter, maximum distance under which two points are considered neighbours.
    :return possibleClusters: result matrix.
    """
    possibleClusters=defaultdict(list)
    tmp=range(len(clusters))
    n=1
    enable=[[[[[] for z in range(n)] for w in range(n)] for x in tmp] for y in tmp]
    for i,clust in enumerate(clusters):
        for j,visit in enumerate(clusters):
#            print "i: "+str(i)+", j: "+str(j)
            for k, partialClust in enumerate(chunksByNumP(clust.points,n)):
                for l, partialClust2 in enumerate(chunksByNumP(visit.points, n)):
                    syncTask(enable[i][j][k][l], partialClust, partialClust2, epsilon)
    a=time.time()
    enable = compss_wait_on(enable)
    print "Long Wait on Time: %.6f" % float(time.time()-a)
    for i, clust in enumerate(clusters): 
        for j, visit in enumerate(clusters):
#            print "i: "+str(i)+", j: "+str(j)
            if myAny(enable[i][j]): 
                possibleClusters[i].append(j)
    return possibleClusters

def myAny(m):
    """
    Checks wether a 2-D matrix contains a True value.
    """
    for i in m:
        for j in i:
            if j[0]:
                return True
    return False

@task(enable=INOUT)
def syncTask(enable, host, visit, epsilon):
    """
    Given two lists, checks wether the distance between the two sets is less than epsilon.
    """
    for p in host:
        for point in visit:
            if np.linalg.norm(point-p) < epsilon:
                return enable.append(1)   
    return enable.append(0)

def partialScan(corePoints, square, epsilon, minPoints, fragData, fragSize, numParts, rangeToEps):
    """
    Looks for all the core points (over minPoints neighbours) inside a certain square.
    :param square: space region where core points are looked for.
    :param epsilon: maximum distance between two points to be considered neighbours.
    :param minPoints: minimum number of neighbours for a point to be considered core point.
    :param fragData: dict containing space partition and its points.
    :param fragSize: length of the side of each hypercube inside fragData.
    :return:    list of all the core points found in square.  
    """
    dim=len(fragData[square][0])
    pointSet=fragData[square]
    pointSetReal=pointSet[:]
    k=rangeToEps[square]
    size = pow(10,len(str(fragSize+1)))
    perm=[]
    #To check all the possible neighbours, all the possible permutations and dimensions are considered.
    for i in range(dim):
        perm.append(range(-k[i],k[i]+1))
    for comb in itertools.product(*perm):
        current=square
        for i in range(dim):
            current=current+comb[i]*math.pow(size,i)
        if current in fragData and current != square:
            pointSet=pointSet+fragData[current]
    #Only points in the current square should be considered as possible centroids.
    for i in range(numParts):
        tmpPS = [p for j,p in enumerate(pointSetReal) if j % numParts == i]
        scanTask(corePoints[i], pointSet, tmpPS, epsilon, minPoints)

@task(corePoints=INOUT) 
def scanTask(corePoints, pointSet, pointSetReal, epsilon, minPoints):
    for point in pointSetReal:
        neighbourPts=0
        for p in pointSet:
            if np.linalg.norm(point-p) < epsilon:
                neighbourPts=neighbourPts+1
                if neighbourPts >= minPoints:
                    corePoints.append(point)
                    break  

def expandCluster(clusters, fragData, epsilon, minPoints, fragSize, numParts, rangeToEps):
    """
    Expands all clusters contained in a list of clusters. Expansion means adding non-core points to the 
    already established clusters.
    """
    #This matrix contains all the possible neighbour points of all clusters
    #Worst case scenario it could have as many elements as points in the Dataset (which imo is not that much)
    addToClust=[[[] for x in range(numParts)] for y in range(len(clusters))]
    clustPoints = np.asarray([i for j in clusters for i in j.points])
    for numClust,clust in enumerate(clusters):
        possiblePoints = neighExpansion(clustPoints, clust, fragData, fragSize, rangeToEps)
        div=int(np.ceil(float(len(possiblePoints)/numParts)))
        if div:
            for part in range(numParts):
                tmpVec = [p for i,p in enumerate(possiblePoints) if i/div == part]
                #addToClust[clust][part]=[]
                expandTask(addToClust[numClust][part], clust, tmpVec, epsilon) 
    addToClust = compss_wait_on(addToClust)
    addPoints(clusters, addToClust)     
            
def neighExpansion(clustPoints, clust, fragData, fragSize, rangeToEps):
    """
    Given a cluster of core points, returns all the points lying in the squares reachable with epsilon distance
    from the clusters square.
    :param cluster: cluster whose neighbour points will be added.
    :param fragData: dict containing space partitioning.
    :param epsilon: maximum distance between two points to be considered neighbours.
    :param minPoints: minimum number of neighbours for a point to be considered core point.
    :param fragSize: length of the side of each hypercube inside fragData.
    """
    squaresNot=clust.square[:]
    pointSet=[]
    for sq in clust.square:
        #Segur que no ho sumem dos cops aleshores?
        pointSet+=fragData[sq][:]
        dim=len(fragData[sq][0])    
        k=rangeToEps[sq]
        size = pow(10,len(str(fragSize+1)))
        perm=[]
        for i in range(dim):
            perm.append(range(-k[i],k[i]+1))
        for comb in itertools.product(*perm):
            current=sq
            for i in range(dim):
                current=current+comb[i]*math.pow(size,i)
            if current in fragData and not(current in squaresNot):
                pointSet=pointSet+fragData[current][:]
                squaresNot.append(current)
    temp=[b for b in pointSet if not(b in clustPoints)]
    return temp

def addPoints(clusters,addToClust):
    """
    """
    start = time.time()
    for numClust, cluster in enumerate(clusters):
        for vec in addToClust[numClust]:
            for pos in vec:
                cluster.addPoint(pos)
    end = time.time()
    print "addPoints method lasts %.6f seconds" % float(end-start)

@task(res = INOUT)
def expandTask(res, cluster, points, epsilon):
    for i,point in enumerate(points):
        for p in cluster.points:
            if np.linalg.norm(point-p) < epsilon:
                res.append(point)

def update(clusters, possibleClusters):
    """
    Given a list of clusters, it reassigns each point to a cluster making sure a point does not appear in more than one 
    different clusters. If so, those clusters will be merged.
    :param clusters: provisional cluster list result of the cluster expansion.
    :return defClusters: final cluster list.
    """
    visited=defaultdict(int)
    defCluster=defaultdict(list)
    #Clusters numbered from 1,...,N
    for key,value in enumerate(clusters):
        #Do a task out of this.
        points=[p for a,b in enumerate(clusters) for p in b.points if a in possibleClusters[key]]
        for p in points:            
            if visited[str(p)] and visited[str(p)] != key:
                #This follows the paper's implementation
                for j in points:
                    if visited[str(j)] != visited[str(p)]:
                        visited[str(j)] = visited[str(p)]
                        defCluster[visited[str(p)]].append(j)
                defCluster.pop(key, None)
                break
            elif not visited[str(p)]:
                visited[str(p)]=key
                defCluster[key].append(p)
    return defCluster
        
if __name__ == "__main__":
    main(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))

