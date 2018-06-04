#Imports
from pandas import read_csv # General Imports
from collections import defaultdict
import os
import itertools
import numpy as np
from pycompss.api.task import task # PyCOMPSs Imports
from classes.Data import Data # DBSCAN Imports

def count_lines(tupla, file_id, is_mn):
    if is_mn:
        path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_id)
    else:
        path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(tupla[0])
    for num, j in enumerate(tupla):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += ".txt"
    with open(tmp_string) as infile:
        for i, line in enumerate(infile):
            pass
        return i+1

def orquestrate_init_data(tupla, file_id, len_data, quocient,
                          res, fut_list, TH_1, is_mn,
                          count_tasks):
    if (len_data/quocient) > TH_1:
        [fut_list,
         count_tasks] = orquestrate_init_data(tupla, file_id, len_data,
                                              quocient*2, res*2 + 0, fut_list,
                                              TH_1, is_mn, count_tasks)
        [fut_list,
         count_tasks] = orquestrate_init_data(tupla, file_id, len_data,
                                              quocient*2, res*2 + 1, fut_list,
                                              TH_1, is_mn, count_tasks)
    else:
        count_tasks += 1
        tmp_f = init_data(tupla, file_id, quocient, res, is_mn)
        fut_list.append(tmp_f)
    return fut_list, count_tasks

@task(returns=1)
def init_data(tupla, file_id, quocient, res, is_mn):
    data_pos = Data()
    if is_mn:
        path = "/gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/"+str(file_id)
    else:
        path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(tupla[0])
    for num, j in enumerate(tupla):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += ".txt"
    df = read_csv(tmp_string, sep=' ', skiprows=lambda x: (x % quocient)
                  != res, header = None)
    data_pos.value = df.values.astype(float)
    tmp_vec = -2*np.ones(np.shape(data_pos.value)[0])
    data_pos.value = [data_pos.value, tmp_vec]
    return data_pos

@task(returns=1)
def merge_task_init_data(*args):
    tmp_data = np.vstack([i.value[0] for i in args])
    return tmp_data