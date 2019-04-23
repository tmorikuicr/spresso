#!/usr/bin/python
# _*_ coding: utf-8 _*_

import numpy as np
import pandas as pd
import sys
import random
import itertools
import os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import pairwise
from numpy import argmax, argsort
import math


# ====================================================
# SOM class
# ====================================================
class cubicSOM():
    
    def __init__(self, exprs, dim, radii, rlen, method, seed=None):

        self.exprs = np.array(exprs.dataframe2numpy())  # matrix of gene expressions
        self.genes = sorted(list(exprs.genenames))      # genes of input data
        self.cells = exprs.cellnames                    # cell names of input data
        self.num_cells = self.exprs.shape[0]            # the number of cells (= 41)
        self.dim = dim                                  # the number of units for each axis
        self.radii = radii                              # learning radius
        self.rlen = rlen                                # the number of learning steps
        self.method = method                            # clustering method ("normal" only)
        self.seed = seed                                # random seed for map initialization

        # ----- set random seed for Numpy -----
        if not seed is None:
            np.random.seed(seed)

        # -------------------------------------
        # generate a grid coordinates "c":
        # self.c = array([[0,0,0],
        #                 [0,0,1],
        #                 [1,0,0],
        #                 [1,0,1],
        #                 [0,1,0],
        #                 [0,1,1],
        #                 [1,1,0],
        #                 [1,1,1]])
        # --------------------------------------
        x, y, z = np.meshgrid(range(self.dim), range(self.dim), range(self.dim))
        self.c = np.hstack((x.flatten()[:, np.newaxis],
                            y.flatten()[:, np.newaxis], 
                            z.flatten()[:, np.newaxis]))

        self.map_vector = self._map_init()  # initialized map
        self.idx2class = {}                 # correspondence table (sample -> (x,y,z))
        self.score = None                   # reconstruction score ( 0 <= score <= 1)
        self.rlenIdx2class = {}             # rlenIdx2class[(rlen, idx)] = (x, y, z)

    def _map_init(self):
        # ------------------------------------------------------------------------------------------
        # Map initialization
        #   * "map_vector" is initialized by zero-matrix whose size is 8 * the number of genes.
        #   * Each unit of map is initialized by the gene expressions of a cell randomly selected.
        # ------------------------------------------------------------------------------------------
        map_vector = np.zeros((self.dim * self.dim * self.dim, self.exprs.shape[1]))
        for i in range(self.dim * self.dim * self.dim):
            r = np.random.choice(range(self.num_cells))
            map_vector[i] = self.exprs[r]
        return map_vector

    def _map_distance(self, bmu):
        # -------------------------------------------------------------------------------------------
        # Computation of distances between units and BMU (best matching unit).
        #   * d_xyz: Euclidean distances on 3D-space between units and BMU
        #   * d_xy: Euclidean distances on 2D-space between units and BMU
        #   * self.c[0:self.dim*self.dim*self.dim,0:2]: Extract x,y coordinates
        #   * bmu[0:2]: Extract x,y coordinates of BMU
        #   * d_xy_designed: If the distance between a unit and BMU on xy-plane is
        #                       greater than 1 (that is, the unit and BMU are located diagonally),
        #                       set infinity as distance, otherwise set zero.
        #   * d_xyz + d_xy_designed: Euclidean distance on 3D-space,
        #                       where the distance between units located diagonally on xy-plane is infinity.
        # -------------------------------------------------------------------------------------------
        d_xyz = np.linalg.norm(self.c - bmu, axis=1)
        d_xy = np.linalg.norm(self.c[0:self.dim*self.dim*self.dim,0:2] - bmu[0:2], axis=1)
        d_xy_designed = [inf if i > 1 else 0 for i in d_xy]
        return(d_xyz + d_xy_designed)

    def batch_learning(self):
        # ---------------------------------------------------------------------------------------------
        # Run Batch SOM
        # ---------------------------------------------------------------------------------------------
        for i in range(self.rlen): # iterate learning steps rlen times
            # print(str(i+1) + 'th learning')
            w = np.zeros((self.dim*self.dim*self.dim, 2, self.exprs.shape[1]))  # temporal weight vector used for update of map
            for j, expr in enumerate(self.exprs):
                bmu = self._best_matching_unit(expr) # search BMU (Eq.(4))
                self.idx2class[j] = bmu              # store BMU for cell j
                self.rlenIdx2class[(i,j)] = bmu      # store BMU for a pair (learning step i, cell j)
                d = self._map_distance(bmu)          # compute distances between units and BMU
                h = self._neighbourhood(i, d)        # update the effect the neighbourhoods are received from BMU (Eq.(3))
                w[:, 0] += h[:, np.newaxis] * expr   # numerator of update function of map_vector (Eq.(5), Summation of h_{c_j(t)i(t)} * x_j_)
                w[:, 1] += h.repeat(self.exprs.shape[1]).reshape((self.dim * self.dim * self.dim, self.exprs.shape[1])) # denominator of update function of map_vector (Eq.(5), Summation of h_{c_j(t)i(t)})
            self.map_vector = w[:, 0] / (w[:, 1] + 10**(-10)) # update of map_vector (Eq.(5))

        # store the xyz-coordinates of unit which cell j belongs to after "rlen" iteration step.
        for j, expr in enumerate(self.exprs):
            bmu = self._best_matching_unit(expr)
            self.idx2class[j] = bmu
            self.rlenIdx2class[(self.rlen, j)] = bmu
        return self.map_vector

    def _best_matching_unit(self, expr):
        # -----------------------------------------------------------------------------------------------------
        # Search BMU
        #   * Eq.(4) on the draft
        #   * norms: Euclidean distance (dissimilarity) between weight of units and gene expression of cell j
        #   * BMU: best matching unit (the most similar unit)
        #   * tuple(self.c[bmu]): xyz-coordinates of BMU
        # -----------------------------------------------------------------------------------------------------
        norms = np.linalg.norm(self.map_vector - expr, axis=1) 
        bmu = np.argmin(norms) 
        return tuple(self.c[bmu])

    def _learning_radius(self, t):
        # -----------------------------------------------------------------------------------------------------
        # Update the learning radius
        #   * the learning radius monotonically decrease as the learning step progresses.
        # -----------------------------------------------------------------------------------------------------
        s = self.radii * (1 - t / float(self.rlen - 1))  # radii decreases lineally
        # s = self.radii * math.exp(-t/float(self.rlen)) # radii decreases exponentially
        return s

    def _neighbourhood(self, t, d):
        # -----------------------------------------------------------------------------------------------------
        # Neighbourhood function
        #   * Eq.(3) on the draft
        # -----------------------------------------------------------------------------------------------------
        if t == (self.rlen - 1):
            h = np.array([1.0 if di == 0.0 else 0.0 for di in d])
        else:
            s = self._learning_radius(t)
            max = 1.0
            min = 0.5
            d = ((max - min) * np.random.rand() + min) * d    # stochastic SOM
            h = np.exp(-d**2/(2*s**2))
        return h
    
    def write_idx2class(self, odir):
        # -----------------------------------------------------------------------------------------------------
        # Output clustering results
        #   * index: index of each cell (0 ~ 40)
        #   * x, y, z: xyz-coordinates ({0,1})
        # -----------------------------------------------------------------------------------------------------
        fout = open(odir + '/idx2class.' + str(self.seed) + '.txt', 'w')
        fout.write('index\tx\ty\tz\n')
        for i, c in self.idx2class.items():
            fout.write(str(i) + '\t' + str(c[0]) + '\t' + str(c[1]) + '\t' + str(c[2]) + '\n')
        fout.close()

    def write_map_vector(self, odir): # TM2018Aug
        # -----------------------------------------------------------------------------------------------------
        # Output map-vector of SOM result
        #   * index vector indicates xyz-coordinates
        # -----------------------------------------------------------------------------------------------------
        df = pd.DataFrame(data=self.map_vector, columns=self.genes,
                          index=["000", "001", "100", "101", "010", "011", "110", "111"]).T
        df.to_csv(odir + '/map_vector.' + str(self.seed) + '.txt', sep='\t')

    def plot_difference_curve(self, pp):
        # -----------------------------------------------------------------------------------------------------
        # Plot difference curve of cell positions through learning steps
        # -----------------------------------------------------------------------------------------------------
        differences = []
        for i in range(self.rlen):
            d = 0
            for j in range(self.num_cells):
                if not self.rlenIdx2class[i, j] == self.rlenIdx2class[i+1, j]:
                    d += 1
            differences.append(d)
        x = np.array(range(1, self.rlen+1))
        y = np.array(differences)
        fig = plt.figure()
        plt.plot(x, y, linewidth=2)
        plt.xlim([-5, self.rlen + 5])
        plt.ylim([-3, 45])
        plt.title("Seed = " + str(self.seed))
        plt.xlabel("Learning steps")
        plt.ylabel("Difference (# cells)")
        plt.grid(True)
        pp.savefig(fig)

    def evaluation(self):
        # -----------------------------------------------------------------------------------------------------
        # Compute reconstruction score
        #   * If the reconstruction is success, self.score becomes 1.
        # -----------------------------------------------------------------------------------------------------
        domain2indexes = {}
        domain2indexes['D1'] = [0,2,9,13,17,21,25,29,33]
        domain2indexes['D2'] = [4,5,7,8,11,12,15,16,19,20,24]
        domain2indexes['D3'] = [23,27,28,31,32,35,36,37,39,40]
        domain2indexes['D4'] = [1,3,6,10,14,18,22,26,30,34,38]
        domain2class = {}
        for domain, indexes in domain2indexes.items():
            sum_x = 0.0
            sum_y = 0.0
            sum_z = 0.0
            for index in indexes:
                sum_x += float(self.idx2class[index][0])
                sum_y += float(self.idx2class[index][1])
                sum_z += float(self.idx2class[index][2])
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
            # If D1 and D4 are positioned diagonally in the xy-direction, the total_score += 1
            if comb == ('D1', 'D4') or comb == ('D4', 'D1'):
                score = 1 if abs(x1 - x2) + abs(y1 - y2) == 2 else 0
            # If D2 and D3 are adjacent in the z-direction, the total_score += 1
            elif comb == ('D2', 'D3') or comb == ('D3', 'D2'):
                score = 1 if abs(z1 - z2) == 1 else 0
            # other cases
            else:
                score = 1 if abs(x1 - x2) + abs(y1 - y2) == 1 else 0
            total_score += score
        self.score = float(total_score) / len(combs)
        print(self.score)

    def write_results(self, fout):
        # -----------------------------------------------------------------------------------------------------
        # Output parameters and reconstruction score
        # -----------------------------------------------------------------------------------------------------
        out = self.method + '\t' + str(self.seed) + '\t' + str(self.radii) + '\t' + str(self.rlen) + '\t' + str(self.score) + '\n'
        fout.write(out)


# ====================================================
# Expression class
# ====================================================
class Expression():
    
    def __init__(self, ifile, odir):
        self.ifile = ifile
        self.odir = odir
        self.dataframe = pd.read_csv(self.ifile, header=0, index_col=0, delimiter='\t').T
        self.cellnames = self.dataframe.index
        self.genenames = self.dataframe.columns
        
    def dataframe2numpy(self):
            return self.dataframe.as_matrix()
    
    def write_idx2cell(self):
        fout = open(odir + 'idx/idx2cell.txt', 'w')
        fout.write('index\tcell\n')
        for i, c in enumerate(self.cellnames):
            fout.write(str(i) + '\t' + str(c) + '\n')
        fout.close()
    
    def write_idx2section(self):
        fout = open(odir + 'idx/idx2section.txt', 'w')
        fout.write('index\tsection\n')
        for i, s in enumerate(self.cellnames):
            fout.write(str(i) + '\t' + str(s[-1]) + '\n')
        fout.close()

# ====================================================
# Main function
# ====================================================
if __name__ == '__main__':

    # ----- command line arguments -----
    ifile = sys.argv[1]
    odir = sys.argv[2]
    if len(sys.argv) - 1 != 2:
        sys.exit('ERROR: # arguments must be 2')

    # ----- parameter setting -----
    dim    = 2              # the number of units for x,y,z-directions
    radii  = 1.5            # initial learning radius 1.5
    rlen   = 100            # the number of learning steps
    inf    = float('inf')

    # ----- preparing the output directories -----
    try:
        os.makedirs(odir)
        os.makedirs(odir + '/idx')
        os.makedirs(odir + '/idx2class')
        os.makedirs(odir + '/map_vector')
    except OSError:
        pass
    odir = odir + '/'

    # ----- generate an instance of gene expression table -----
    exprs = Expression(ifile, odir)
    exprs.write_idx2cell()      # output cell information
    exprs.write_idx2section()   # output section information

    # ----- open scores.txt file -----
    ofile_scores = odir.rstrip('/') + '/scores.txt'
    fout_scores = open(ofile_scores, 'w')

    # ----- display parameters -----
    print('----- Cubic SOM -----')
    print('Dim: ' + str(dim) + '*' + str(dim) + '*' + str(dim))
    print('Initial learning radius: ' + str(radii))
    print('Learning steps: ' + str(rlen))
    print('---------------------')

    # ----- start SOM -----
    with PdfPages(odir + 'diff_curves.pdf') as pp:
        for i in range(100):
            som = cubicSOM(exprs, dim=dim, radii=radii, rlen=rlen, method='normal', seed=i)    # initialization
            som.batch_learning()    # batch learning
            som.write_idx2class(odir + '/idx2class')    # output SOM results
            som.write_map_vector(odir + '/map_vector')  # output map-vector of SOM # TM2018Aug
            som.evaluation()    # compute reconstruction score
            som.write_results(fout_scores)  # output evaluation results
            som.plot_difference_curve(pp)

    fout_scores.close()

