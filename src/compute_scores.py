#!/usr/bin/python
# _*_ coding: utf-8 _*_

import pandas as pd
import numpy as np
import sys
import os
import glob
import itertools

def compute_success_rates(ifile, domain2indexes):
    idx2class = pd.read_csv(ifile, header=0, index_col=0, delimiter = '\t')
    idx2class = idx2class.as_matrix()
    domain2class = {}
    for domain, indexes in domain2indexes.items():
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        for index in indexes:
            sum_x += float(idx2class[index][0])
            sum_y += float(idx2class[index][1])
            sum_z += float(idx2class[index][2])
        center_x = 0 if sum_x/len(indexes) < 0.5 else 1
        center_y = 0 if sum_y/len(indexes) < 0.5 else 1
        center_z = 0 if sum_z/len(indexes) < 0.5 else 1
        domain2class[domain] = [center_x, center_y, center_z]
    total_score = 0.0
    combs = list(itertools.combinations(domain2class.keys(),2))
    for comb in combs:
        x1 = domain2class[comb[0]][0]
        y1 = domain2class[comb[0]][1]
        z1 = domain2class[comb[0]][2]
        x2 = domain2class[comb[1]][0]
        y2 = domain2class[comb[1]][1]
        z2 = domain2class[comb[1]][2]           
        if comb == ('D1', 'D4') or comb == ('D4', 'D1'):
            score = 1 if abs(x1 - x2) + abs(y1 - y2) == 2 else 0
        elif comb == ('D2', 'D3') or comb == ('D3', 'D2'):
            score = 1 if abs(z1 - z2) == 1 else 0   
        else:
            score = 1 if abs(x1 - x2) + abs(y1 - y2) == 1 else 0
        total_score += score
    score = float(total_score) / len(combs)
    if score >= 1:
        score = 1
    else:
        score = 0
    return(score)

def compute_var(ifile, domain2indexes):
    idx2class = pd.read_csv(ifile, header=0, index_col=0, delimiter = '\t')
    idx2class = idx2class.as_matrix()
    domain2class = {}
    total_var = 0.0
    for domain, indexes in domain2indexes.items():
        x = []
        y = []
        z = []
        for index in indexes:
            x.append(float(idx2class[index][0]))
            y.append(float(idx2class[index][1]))
            z.append(float(idx2class[index][2]))
        var_x = np.var(np.array(x), ddof=1)
        var_y = np.var(np.array(y), ddof=1)
        var_z = np.var(np.array(z), ddof=1)
        total_var += (var_x + var_y + var_z) / 3.0
    total_var = total_var / float(len(domain2indexes))
    return(total_var)


#====================================================
# Main function
#====================================================
if __name__ == '__main__':
    idir = sys.argv[1]
    ofile = sys.argv[2]
    
    domain2indexes = {}
    domain2indexes['D1'] = [0,2,9,13,17,21,25,29,33]
    domain2indexes['D2'] = [4,5,7,8,11,12,15,16,19,20,24]
    domain2indexes['D3'] = [23,27,28,31,32,35,36,37,39,40]
    domain2indexes['D4'] = [1,3,6,10,14,18,22,26,30,34,38]

    num_cell = len(domain2indexes['D1']) + len(domain2indexes['D2']) + len(domain2indexes['D3']) + len(domain2indexes['D4'])
    
    epsilon = 10**(-6)
    itr = 100
    
    fout = open(ofile, 'w')
    fout.write('sample\tnum_cell\tsuccess_rate\tvar\n')
    directories = os.listdir(idir)
    for dir in sorted(directories):
        print(dir)
        total_success_rate = 0.0
        total_var = 0.0
        total_score = 0.0
        for i in range(itr):
            ifile = idir + '/' + dir + '/idx2class/idx2class.' + str(i) + '.txt'
            sr = compute_success_rates(ifile, domain2indexes)
            cs = compute_var(ifile, domain2indexes)
            total_success_rate += sr
            total_var += cs
        success_rate = total_success_rate / itr
        var = total_var / itr
        fout.write(dir + '\t' + str(num_cell) + '\t' + str(success_rate) + '\t' + str(var) + '\n')
    fout.close()