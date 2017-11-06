# from pycompss.api.task import task
# from pycompss.api.parameter import INOUT
import sys
import numpy as np
from itertools import product


# @task(data_pos=INOUT)
def init_data(comb, path, num_points_max, num_centers):
    data = []
    dim = len(comb)
    centers = np.random.sample((num_centers, dim))
    std = np.random.sample((num_centers, dim, dim))
    for i in std:
        i = np.dot(i, i.transpose())/10
    for i in range(len(centers)):
        for j in range(num_points_max):
            tmp = np.random.multivariate_normal(centers[i], std[i], size=(1))
            for num, k in enumerate(comb):
                if (tmp[0][num] <= float(k)/10 or
                        tmp[0][num] > float(k + 1)/10):
                    break
            else:
                data.append(tmp)
    if len(data) > 0 and False:
        data = np.vstack(data)
        tmp_vec = -2*np.ones(np.shape(data)[0])
        data = [data, tmp_vec]
        np.savetxt(path, data)
        return len(tmp_vec)
    else:
        data = np.random.sample((1, dim))
        for num, k in enumerate(comb):
            data[0][num] = np.random.uniform(low=float(k)/10,
                                             high=float(k+1)/10)
        print comb
        print data
        tmp_vec = -2*np.ones(np.shape(data)[0])
        data = [data, tmp_vec]
        return 1
#    return data_pos


def main(file_count, dimensions):
    num_points_max = 10000
    num_centers = len(dimensions)
    path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"
    perm = []
    for i in range(len(dimensions)):
        perm.append(range(dimensions[i]))
    for comb in product(*perm):
        tmp_string = path+str(file_count)+"/"+str(comb[0])
        for num, j in enumerate(comb):
            if num > 0:
                tmp_string += "_"+str(j)
        tmp_string += ".txt"
        init_data(comb, path, num_points_max, num_centers)
        print tmp_string


if __name__ == "__main__":
    main(int(sys.argv[1]), list(sys.argv[2]))
