#DBSCAN v3. Carlos Segarra.
#TODO we do not re-use calculated distances, we just try to minimize the number of distances we have to compute.

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
			#Consider modifying ClusterOld init TODO
			#Class Cluster is defined at /classes/cluster.py
			tmp=([point],[square])
			tmpc=Cluster()
			tmpc.add(*tmp)
			clusters.append(tmpc)

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
	print "Density Based Scan started."
	normData=normalizeData(dataFile) 
        dataset=np.loadtxt(normData)
	fragData=partitionData(dataset, float(fragSize)) 
	clusters=[]	
	print "Starting partial scan..."
	#For the moment this is arbitrarily fixed to 10, user fixed/learnt parameter? TODO
	numParts=10
	for j in range(numParts):
		clustersSub=[]
		fragDataSub={key: value for i, (key, value) in enumerate(fragData.viewitems()) if i % numParts == j}
		for square in fragDataSub:
			coreP=partialScan(square, epsilon, minPoints, fragData, fragSize) 
			mergeCluster(clustersSub, coreP, square, epsilon) 
		syncClusters(clusters, clustersSub, epsilon)
	#Might be a barrier instead?
	clusters=compss_wait_on(clusters) 
	for clust in clusters: 
		expandCluster(clust, fragData, epsilon, minPoints, fragSize) 
	defClusters=update(clusters)
	print "DBSCAN Algorithm finished succesfully."
	return defClusters 

@task(clusters = INOUT) 
def syncClusters(clusters, clustersSub, epsilon):
	"""
	Synchronizes the general list of clusters (master) with the ones found by a speciphic job.
	:inout	clusters: list of all clusters and their points.
	:param 	clustersSub: list of clusters to be merged with the general list.
	:param	epsilon: input parameter, maximum distance under which two points are considered neighbours.
	"""
	for candidate in clustersSub:
		possibleClusters=[]
		for host in clusters: 
			for p in candidate.points:
				exitFlag=False
				for point in host.points:
					if np.linalg.norm(point-p) < epsilon:
						possibleClusters.append(host)
						if len(possibleClusters) > 1:
							clusters.remove(host)
						exitFlag=True
						break
				if exitFlag==True:
					break
		if len(possibleClusters) > 0:
			master=possibleClusters.pop(0)
			master.merge(candidate)
			for i in possibleClusters:
				master.merge(i)
		else:
			tmp=(candidate.points, candidate.square)
			tmpc=Cluster()
			tmpc.add(*tmp)
			clusters.append(tmpc)		
		 
@task(returns=list) 
def partialScan(square, epsilon, minPoints, fragData, fragSize):
	"""
	Looks for all the core points (over minPoints neighbours) inside a certain square.
	:param	square: space region where core points are looked for.
	:param 	epsilon: maximum distance between two points to be considered neighbours.
	:param 	minPoints: minimum number of neighbours for a point to be considered core point.
	:param 	fragData: dict containing space partition and its points.
	:param	fragSize: length of the side of each hypercube inside fragData.
	:return: 	list of all the core points found in square.	
	"""
	corePoints=[]
	dim=len(fragData[square][0])
	pointSet=fragData[square]
	pointSetReal=pointSet[:]
	#To check wether a point is a core point or not, we check all the square's neighbour squares and their points.
	k=int(max(epsilon/fragSize, 1))
	perm=[]
	for i in range(dim):
		perm.append(range(-k,k+1))
	for comb in itertools.product(*perm):
		current=square
		for i in range(dim):
			current=current+comb[i]*math.pow(1/fragSize,i)
	        if current in fragData and current != square:
                        pointSet=pointSet+fragData[current]
	#Only points in the current square should be considered as possible centroids.
	for point in pointSetReal:
		#Still re-calculating plenty of distances, TODO.
		#In order not to recalculate distances, let's do it the long way round.
		neighbourPts=0
		for p in pointSet:
			if np.linalg.norm(point-p) < epsilon:
				neighbourPts=neighbourPts+1
				if neighbourPts >= minPoints:
					corePoints.append(point)
					break	
	return np.asarray(corePoints)

@task(c=INOUT) 		
def expandCluster(cluster, fragData, epsilon, minPoints, fragSize):
	"""
	Given a cluster of core points, adds to the cluster all the core points' neighbours.
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
        for point in temp:
        	for p in cluster.points:
                	if np.linalg.norm(point-p) < epsilon:
                               cluster.addPoint(point)
                               break

@local 
def update(clusters):
	"""
	Given a list of clusters, it reassigns each point to a cluster making sure a point does not appear in more than one 
	different clusters. If so, those clusters will be merged.
	:param	clusters: provisional cluster list result of the cluster expansion.
	:return defClusters: final cluster list.
	"""
	visited=defaultdict(int)
	defCluster=defaultdict(list)
	random.shuffle(clusters)
	#Clusters numbered from 1,...,N
	i = 1
	for value in clusters:
		for p in value.points:			
			if visited[str(p)] and visited[str(p)] != i:
				#This follows the paper's implementation
				for j in value.points:
					defCluster[visited[str(p)]].append(j)
				defCluster.pop(i, None)
				break
			visited[str(p)]=i
			defCluster[i].append(p)
		i+=1		
	return defCluster

def main(dataFile, fragSize, epsilon, minPoints):
	plot = False
	defCluster=DBSCAN(dataFile, fragSize, epsilon, minPoints) 
       # if plot:
	#	import matplotlib.pyplot as plt
        # 	from matplotlib import colors
	#        colours = [hex for (name, hex) in colors.cnames.iteritems()]
	#        fig, ax = plt.subplots()
		#ax.set_xticks(np.arange(0, 2, fragSize))
		#ax.set_yticks(np.arange(0, 2, fragSize))
		#ax.tick_params(colors='lightgrey')
		#For debuggin it might be usefull to let ticks on, otherwise they might get really annoying.
	#	ax.set_xticklabels([])
        #        ax.set_yticklabels([])
	#        for c in defCluster:
        #     		ax.scatter([p[0] for p in defCluster[c]], [p[1] for p in defCluster[c]    ], color=colours[c], s=1)
	#	plt.grid()
        # 	plt.savefig('plot/clusters.png')
	#	plt.close()
	f= open('outDBSCAN.txt', 'w')
	for num in defCluster:
		f.write('Cluster '+str(num)+':\n')
		#defCluster[num]=sorted(defCluster[num])
		for point in defCluster[num]:
			f.write(str(point)+'\n')	
	f.close()

if __name__ == "__main__":
	main(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))

