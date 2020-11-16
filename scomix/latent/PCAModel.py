#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 15, 2020
# License: MIT License

# PCA Plots

#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files
#   (the "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject
#   to the following conditions:

#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.

#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
#   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
#   TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
#   SHALL THE AUTHORS, COPYRIGHT HOLDERS OR ARIZONA BOARD OF REGENTS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
#   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#   OR OTHER DEALINGS IN THE SOFTWARE.

from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import gc
import logging
from datetime import datetime

from ..rna import RnaCounts

_LOGGER = logging.getLogger(__name__)

class PCAModel:
    """
    PCA Model of nxm data set

    Parameters
    -------------
    matrix: `pd.DataFrame`
        `matrix to perform principal component on

    n_components: `int`
        Number of components for PCA

    transform: `str`
        Transform method for transforming data matrix

    svd_solver: `str`
        Types of Singular Value Decomposition (SVD) Solver

    seed: `int`
        Seed for reproducibility

    Attributes
    -----------
    n_components: `int`
        Number of components for PCA

    transform: `str`
        Transform method for transforming data matrix

    svd_solver: `str`
        Types of Singular Value Decomposition (SVD) Solver

    seed: `int`
        Seed for reproducibility

    rnacounts: `RnaCounts`
        RnaCounts object

            
    """

    def __init__(self, matrix: RnaCounts, n_components = 25, transform = None, svd_solver = "randomized", seed = 0) -> None:
        """        
        """

        self.n_components = n_components
        self.rnacounts = matrix
        self.transform = transform
        self.svd_solver = svd_solver
        self.seed =seed


    def fit(self, transpose = False):
        """
        Train the model and return principal components

        Parameters
        ------------
        """

        if not np.issubdtype(self.rnacounts.cellxgene.values.dtype, np.floating):
            self.rnacounts.cellxgene = self.rnacounts.cellxgene.astype(np.float32, copy=False)
            _LOGGER.info("Converted matrix to float32 data type.")


        self.rnacounts.scale(self.rnacounts.median, inplace = True)


        transform = "freeman-tukey"
        if self.transform is not None:
            transform = self.transform

        tranformed_counts = self.rnacounts.transform(func =transform, inplace= True)
        
        gc.collect()

        pca_model = PCA(n_components  = self.n_components, svd_solver=self.svd_solver, random_state=self.seed)

        # make X contiguous
        X = np.array(tranformed_counts.values, order = 'C', copy = False)

        if transpose:
            _LOGGER.info("Calculating PCA of transposed matrix.")
            X = X.transpose()

        Y = pca_model.fit_transform(X)
        Y = np.array(Y, order='F', copy=False)

        total_variance = X.var(axis = 0, ddof = 1).sum()
        explained_variance = Y.var(axis =0, ddof =1)
        explained_variance_ratio = explained_variance/total_variance

        dim_labels = [
                'PC_%d (%.3f)' % (c+1, v)
                for c, v in zip(
                    range(self.n_components), explained_variance_ratio)]


        if transpose:
            pc_scores = pd.DataFrame(Y, index=self.rnacounts.cells, columns=dim_labels)
        else:
            pc_scores = pd.DataFrame(Y, index=self.rnacounts.genes, columns=dim_labels)

        self.rnacounts.PC = pc_scores

        self.rnacounts.records["PC"] ={}

        self.rnacounts.records["PC"]["sequence"] = self.rnacounts.sequence 
        self.rnacounts.sequence  =  self.rnacounts.sequence  + 1
        self.rnacounts.records["PC"]["timestamp"] =datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S.%f')
        self.rnacounts.records["PC"]["n_components"] = self.n_components
        self.rnacounts.records["PC"]["transform"] = self.transform
        self.rnacounts.records["PC"]["svd_solver"] = self.svd_solver
        self.rnacounts.records["PC"]["seed"] = self.seed
        self.rnacounts.records["PC"]["tranposed"] = transpose

        return pc_scores
        
        
        



