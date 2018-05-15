import os
from itertools import product
from collections import defaultdict
from task_src.init_data import count_lines, orquestrate_init_data
from task_src.init_data import merge_task_init_data


class Square(object):

    def __init__(self, coord, epsilon, dimensions):
        self.coord = coord
        self.epsilon = epsilon
        self.dimensions = dimensions
        self.offset = defaultdict()
        self.len = defaultdict()
        self.__neigh_squares_query()


    def __neigh_squares_query(self):
        # Only multiples of 10 are supported as possible number
        # of squares per side.
        dim = len(self.coord)
        neigh_squares = []
        border_squares = [int(min(max(self.epsilon*i, 1), i-1)) for i in 
                          self.dimensions]
        perm = []
        for i in range(dim):
            perm.append(range(-border_squares[i], border_squares[i] + 1))
        for comb in product(*perm):
            current = []
            for i in range(dim):
                if self.coord[i] + comb[i] in range(self.dimensions[i]):
                    current.append(self.coord[i]+comb[i])
            if len(current) == dim and current != list(self.coord):
                neigh_squares.append(tuple(current))
        neigh_squares.append(tuple(self.coord))
        self.neigh_sq_id = tuple(neigh_squares)


    def init_data(self, file_id, is_mn, TH_1, count_tasks):
        fut_list = []
        prev = 0
        count_tasks
        for comb in self.neigh_sq_id:
            self.offset[comb] = prev
            self.len[comb] = count_lines(comb, file_id, is_mn)
            prev += self.len[comb]
            fut_list, count_tasks = orquestrate_init_data(comb, file_id,
                                                          self.len[comb], 1,
                                                          0, fut_list, TH_1,
                                                          is_mn, count_tasks)
        count_tasks += 1
        self.points = merge_task_init_data(*fut_list)
        return count_tasks

    def get_points(self, comb):
        return self.points[self.offset[comb]: self.offset[comb] + self.len[comb]]

