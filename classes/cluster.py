import numpy
class Cluster(object):

	def __init__(self,  *args, **kwargs):
        	super(Cluster, self).__init__(*args, **kwargs)
        	self.square=[]
        	self.points=[]

	def add(self, points, square):
		self.points=numpy.asarray(points)
		self.square=square

	def addPoint(self,p):
		self.points=numpy.vstack((self.points, numpy.asarray(p)))

	def addSquare(self,s):
		self.square=list(set(self.square) | set(s))

	def merge(self, q):
                if len(self.points):
                    for p in q.points:
                            self.addPoint(p)
                    self.addSquare(q.square)
                else: 
                    tmp=(q.points, q.square)
                    self.add(*tmp)

	def notInCluster(self, point):
        	"""
         	Returns wether a point belongs to a given cluster or not.
         	"""
         	tmp=self.points[:]-point
         	tmp=tmp[:]==[0,0]
         	return not numpy.any(numpy.multiply(tmp[:,0],tmp[:,1]))


class ClusterOld(object):

	def __init__(self,  *args, **kwargs):
        	super(ClusterOld, self).__init__(*args, **kwargs)
        	self.square=0
        	self.points=[]

	def add(self, points, square):
		self.points=numpy.asarray(points)
		self.square=square

	def addPoint(self,p):
		self.points=numpy.vstack((self.points, numpy.asarray(p)))

	def addSquare(self,s):
		self.square=list(set(self.square) | set(s))

	def merge(self, q):
		for p in q.points:
			self.addPoint(p)
