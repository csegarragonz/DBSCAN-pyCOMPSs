# DBSCAN-PyCOMPSs
## 1. Introduction
DBSCAN for PyCOMPSs is a distributed approach to the well known clustering algorithm proposed for the first time in 1996 [here](https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf "Original Implementation"). The algorithm parallelisation is  performed by the [COMPSs framework](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/ "COMPSs Homepage"). 
### How it works?
  1. Initially, the dataset is chunked according to the point density, obtaining equally loaded chunks. 
  2. A core point retrieval is performed simultaneously at chunks of each chunk.
  3. Results are synced using an adjacency matrix.
  4. For each resulting cluster, a reachable point retrieval is performed and results lastly updated.
![](https://github.com/csegarragonz/DBSCAN-pyCOMPSs/blob/master/img/animation.gif "Main Workflow")
## 2. Files
## 3. Requirements
  1. Python 2.7.x (with NumPy) **_COMPSs won't work with Python 3.x_**
  2. [COMPSs 2.1.](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/downloads "Download COMPSs") 
  3. Matplotlib (in case plotting is desired)
