import subprocess
import os
import time
import csv

path = "/home/csegarra/DBTRASH/mem/COMPSsWorker/"
orig_dir = os.path.realpath(__file__)
orig_dir = orig_dir[0:orig_dir.rfind("/")]

os.chdir(path)
values = []
mem = 1
while mem != "8.0K":
    try:
        output = subprocess.check_output(['du', '-h'])
    except subprocess.CalledProcessError as e:
        pass
    beg = output.rfind("\n", 0, output.rfind("\n"))
    if beg == -1:
        beg = 0
    else:
        beg += 1
    end = output.rfind("\n") - 2
    mem = output[beg:end]
    values.append(mem)
    time.sleep(2)
os.chdir(orig_dir)
tmp_vec = [2*i for i in range(len(values))]
values = [tmp_vec, values]
values = zip(*values)
with open('mem_consum.txt', 'w') as myfile:
    wr = csv.writer(myfile)
    for _list in values:
        wr.writerow(_list)
