def sync_clusters(data, adj_mat, epsilon, coord, neigh_sq_loc, quocient,
                  res, len_neighs, *args):
    # TODO: change *args
    adj_mat_copy = [deque() for _ in range(max(adj_mat[0], 1))]
    neigh_data = [np.vstack([i.value[0] for i in args]),
                  np.concatenate([i.value[1] for i in args])]
    neigh_data = [np.vstack([i for num, i in
                             enumerate(neigh_data[0])
                             if ((num % quocient) == res)]),
                  np.array([i for num, i in enumerate(neigh_data[1])
                            if ((num % quocient) == res)])]
    print "Jelou: "+str(len(data.value[0]))+" "+str(len(neigh_data[0]))
    tmp_unwrap = [neigh_sq_loc[i] for i in range(len(neigh_sq_loc))
                  for j in range(len(args[i].value[1]))]
    tmp_unwrap = [i for num, i in enumerate(tmp_unwrap) if
                  ((num % quocient) == res)]
    for num1, point1 in enumerate(data.value[0]):
        current_clust_id = int(data.value[1][num1])
        if current_clust_id > -1:
            tmp_vec = (np.linalg.norm(neigh_data[0]-point1, axis=1) -
                       epsilon) < 0
            for num2, poss_neigh in enumerate(tmp_vec):
                loc2 = tmp_unwrap[num2]
                if poss_neigh:
                    clust_ind = int(neigh_data[1][num2])
                    adj_mat_elem = [loc2, clust_ind]
                    if ((clust_ind > -1) and
                        (adj_mat_elem not in
                            adj_mat_copy[current_clust_id])):
                        adj_mat_copy[current_clust_id].append(adj_mat_elem)
    return adj_mat_copy

if __name__ == "__main__":
    # file = 3/data.txt
    file = sys.argv[0]
    epsilon=0.1
    quocient = 1
    res = 0
    len_neighs = 0
    coord = 0
    data_pos = Data()
    path = "~/DBSCAN/data/"+str(file_id)
    path = os.path.expanduser(path)
    tmp_string = path+"/"+str(tupla[0])
    for num, j in enumerate(tupla):
        if num > 0:
            tmp_string += "_"+str(j)
    tmp_string += ".txt"
    from pandas import read_csv
    df = read_csv(tmp_string, sep=' ', header = None)
    data_pos.value = df.values.astype(float)
    tmp_vec = -2*np.ones(np.shape(data_pos.value)[0])
    data_pos.value = [data_pos.value, tmp_vec]
