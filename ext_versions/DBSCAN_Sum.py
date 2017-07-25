#DBSCAN Tasques curtes. Carlos Segarra. THIS IS A PRELIMINARY VERSION
#This is Est with expand improvement.
#carlos.segarra @ bsc.es

#Imports
from pycompss.api.task import task
from pycompss.api.parameter import *
from pycompss.api.api import compss_wait_on
from pycompss.api.local import local
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
	:param	dataFile: name of the original data file.
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

def partitionData(dataset, fragSize):
	"""
	Divide the dataset in fragSize squares.
	:param 	dataset: data loaded from the normalized datafile. Raw.
	:param 	fragSize: size of the space partition.
	:return	fragData: dict where the key is a unique identifier for each hypercube in the grid and values are the points
				there contained.
	"""
	fragData=defaultdict(list)
	size=float(1/fragSize)
	#For the moment, size must be a multiple of 10. And all sides are of equal length.
	#Hence only negative powers of 10 are valid fragSizes.
	dim=len(dataset[0])
	#For each point we create its key. A key is a number that from right to left contains the -size first digits
	#of the point coordinate in each dimension. 
	for point in dataset:
		key = 0.0
		for i in range(dim):
			key += (math.pow(size,i))*(math.floor(point[i]/fragSize)) 
		fragData[key].append(point)
	return fragData

@task(clusters = INOUT) 
def mergeCluster(clusters, corePoints, square, epsilon):
	"""
	Given a list of core points, decides wether each core point forms a new cluster or belongs to an existing one.
	:inout 	clusters: 	list of clusters found.
	:param	corePoints: 	list of all the core points found in square.
	:param 	square:		working square.
	:param	epsilon: 	input parameter.
	"""
	for point in corePoints:
		#Check wether for each new core point, it already exists a cluster that lies in our point radius.
		possibleClusters=[]
		for clust in clusters:
			for clustP in clust.points:
				if np.linalg.norm(point-clustP) < epsilon:
					possibleClusters.append(clust)
					#Alternitavely, all clusters in possibleClusters could be removed and then the master
					#appended again. I would have to think wether copies would be safely transfered.
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

#def dumpPoints(matrix, clusters):
#	for i,k in enumerate(matrix):
#		clusters[i]=[m for n in matrix[i] for m in n]

def DBSCAN(dataFile, fragSize, epsilon, minPoints):
	"""
	Main DBSCAN algorithm.
	:param	dataset: normalized dataset. Is here because we may supress the main(). TODO
	:param 	epsilon: maximum distance for two points to be considered neighbours.
	:param 	minPoints: minimum number of neighbours for a point to be considered core point.
	:param	fragData: dict with the partitioned space and the points classified.
	:param	fragSize: size used for space partition.
	:param 	numParts: number of parts in which fragData is divided for processing.
	:return	defClusters: result of the application, dict where values are points corresponding to the <key> cluster.
	"""
	start=time.time()
	print "Density Based Scan started."
	#numParts is used to avoid having too many tasks. numParts=numProcessesGiven*k
	numParts=20
	normData=normalizeData(dataFile) 
        dataset=np.loadtxt(normData)
	fragData=partitionData(dataset, float(fragSize)) 
	print "Starting partial scan..."
	#numParts must be a power of 2
	#The neeighbour retreival is divided in -numParts different branches as to reduce data dependencies.
	clusters=[[] for x in range(numParts)]
#	pointsToClust=[[[[] for x in range(numParts)] for y in range(len(fragData)] for z in range(numParts)]
	for j in range(numParts):
		fragDataSub={key: value for i, (key, value) in enumerate(fragData.viewitems()) if i % numParts == j}
		corePMatrix=[[[] for x in range(numParts)] for y in range(len(fragDataSub))]
		for k,(square,value) in enumerate(fragDataSub.viewitems()):
			partialScan(corePMatrix[k], square, epsilon, minPoints, fragData, fragSize, numParts)
			for n in range(numParts):
				mergeCluster(clusters[j], corePMatrix[k][n], square, epsilon)
#		for square in fragDataSub:
#			coreP=partialScan(square, epsilon, minPoints, fragData, fragSize) 
#			mergeCluster(clusters[j], coreP, square, epsilon) 
#	pointsToClust=compss_wait_on(pointsToClust)
#	dumpPoints(clusters, pointsToClust)
	clusters=compss_wait_on(clusters)
	clusters=[j for i in clusters for j in i] 
	print len(clusters)
	#Two loops since nesting is not supported
	expandCluster(clusters, fragData, epsilon, minPoints, fragSize, numParts)
	possibleClusters=syncClusters(clusters, epsilon, numParts)
	print possibleClusters
	defClusters=update(clusters, possibleClusters)
	print "DBSCAN Algorithm finished succesfully."
	end = time.time()
	print "Time elapsed: %.6f" % float(end-start)
	sys.stdout.flush()
	sys.stderr.flush()
	return defClusters 

def syncClusters(clusters, epsilon, numParts):
	"""
	Returns a matrix of booleans. Pos [i,j]=1 <=> clusters -i and -j should be merged.
	:inout	clusters: list of all clusters and their points.
	:param 	numParts: list of clusters to be merged with the general list.
	:param	epsilon: input parameter, maximum distance under which two points are considered neighbours.
	:return possibleClusters: result matrix.
	"""
	possibleClusters=defaultdict(list)
	tmp=range(len(clusters))
	n=2
	enable=[[[[0 for z in range(n)] for w in range(n)] for x in tmp] for y in tmp]
	for i,clust in enumerate(clusters):
		for j,visit in enumerate(clusters):
			for k, partialClust in enumerate(chunksByNumP(clust.points,n)):
				for l, partialClust2 in enumerate(chunksByNumP(visit.points, n)):
					enable[i][j][k][l] = syncTask(partialClust, partialClust2, epsilon)
	enable = compss_wait_on(enable)
	for i, clust in enumerate(clusters): 
		for j, visit in enumerate(clusters):
				if myAny(enable[i][j]): 
					possibleClusters[i].append(j)
	return possibleClusters

def myAny(m):
	"""
	Checks wether a 2-D matrix contains a True value.
	"""
	for i in m:
		for j in i:
			if j:
				return True
	return False

@task(returns=bool)
def syncTask(host, visit, epsilon):
	"""
	Given two lists, checks wether the distance between the two sets is less than epsilon.
	"""
	for p in host:
		for point in visit:
			if np.linalg.norm(point-p) < epsilon:
				return True	
	return False
		 
def partialScan(corePoints, square, epsilon, minPoints, fragData, fragSize, numParts):
	"""
	Looks for all the core points (over minPoints neighbours) inside a certain square.
	:param	square: space region where core points are looked for.
	:param 	epsilon: maximum distance between two points to be considered neighbours.
	:param 	minPoints: minimum number of neighbours for a point to be considered core point.
	:param 	fragData: dict containing space partition and its points.
	:param	fragSize: length of the side of each hypercube inside fragData.
	:return: 	list of all the core points found in square.	
	"""
	dim=len(fragData[square][0])
	pointSet=fragData[square]
	pointSetReal=pointSet[:]
	#To check wether a point is a core point or not, we check all the square's neighbour (reachable in
	# epsilon distance) squares and their points.
	k=int(max(epsilon/fragSize, 1))
	perm=[]
	#To check all the possible neighbours, all the possible permutations and dimensions are considered.
	for i in range(dim):
		perm.append(range(-k,k+1))
	for comb in itertools.product(*perm):
		current=square
		for i in range(dim):
			current=current+comb[i]*math.pow(1/fragSize,i)
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
#	corePoints = np.asarray(corePoints)

def expandCluster(clusters, fragData, epsilon, minPoints, fragSize, numParts):
	"""
	Expands all clusters contained in a list of clusters. Expansion means adding non-core points to the 
	already established clusters.
	"""
	#This matrix contains all the possible neighbour points of all clusters
	#Worst case scenario it could have as many elements as points in the Dataset (which imo is not that much)
	addToClust=[[[] for x in range(numParts)] for y in range(len(clusters))]
	for numClust,clust in enumerate(clusters):
		possiblePoints = neighExpansion(clust, fragData, epsilon, minPoints, fragSize, numParts)
		div=int(np.ceil(float(len(possiblePoints)/numParts)))
		#addToClust=[[] for x in range(numParts)]
		for part in range(numParts):
			tmpVec = [p for i,p in enumerate(possiblePoints) if i/div == part]
			#addToClust[clust][part]=[]
			expandTask(addToClust[numClust][part], clust, tmpVec, epsilon) 
	addToClust = compss_wait_on(addToClust)
	#could this be a task aswell?
	addPoints(clusters, possiblePoints, addToClust, div) 	
			
def neighExpansion(cluster, fragData, epsilon, minPoints, fragSize, numParts):
	"""
	Given a cluster of core points, returns all the points lying in the squares reachable with epsilon distance
	from the clusters square.
	:param	cluster: cluster whose neighbour points will be added.
	:param 	fragData: dict containing space partitioning.
	:param 	epsilon: maximum distance between two points to be considered neighbours.
	:param 	minPoints: minimum number of neighbours for a point to be considered core point.
	:param	fragSize: length of the side of each hypercube inside fragData.
	"""
	squaresNot=cluster.square[:]
	for sq in cluster.square:
		pointSet=fragData[sq]
		dim=len(fragData[sq][0])	
		k=int(max(epsilon/fragSize, 1))
		perm=[]
		for i in range(dim):
			perm.append(range(-k,k+1))
		for comb in itertools.product(*perm):
			current=sq
			for i in range(dim):
				current=current+comb[i]*math.pow(1/fragSize,i)
		        if current in fragData and not(current in squaresNot):
                               	pointSet=pointSet+fragData[current]
				squaresNot.append(current)
        temp=[b for b in pointSet if not(b in cluster.points)]
	return temp

def addPoints(clusters,possiblePoints, addToClust, div):
	"""
	"""
	start = time.time()
	for numClust,cluster in enumerate(clusters):
		for i, vec in enumerate(addToClust[numClust]):
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

#@local 
def update(clusters, possibleClusters):
	"""
	Given a list of clusters, it reassigns each point to a cluster making sure a point does not appear in more than one 
	different clusters. If so, those clusters will be merged.
	:param	clusters: provisional cluster list result of the cluster expansion.
	:return defClusters: final cluster list.
	"""
	visited=defaultdict(int)
	defCluster=defaultdict(list)
#	random.shuffle(clusters)
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

