#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 15, 2020
# License: MIT License

# Neat is a method of denoising the cellxgene dataset

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

import logging

_LOGGER = logging.getLogger(__name__)

import numpy as np
import pandas as pd
from ..rna import RnaCounts, rna
from ..latent import PCAModel

from scipy.stats import poisson

class NEAT:
    """
    NEAT is used for removing technical noise from cellxgene matrix

    """
    def __init__(self, rnacounts: RnaCounts) -> None:
        
        self.rnacounts = rnacounts
        self.n_components = None
        pass

    def get_latent_dimension(self, max_component = 100, variance_multiplier = 2.0):
        
        self.rnacounts.scale(self.rnacounts.median, inplace = True)
        mean = self.rnacounts.cellxgene.mean(axis = 0)

        empty_matrix = np.empty(self.rnacounts.cellxgene.shape, dtype=np.uint16)
        for k in range(0, self.rnacounts.n_cells):
            
            # generate synthetic counts drawn from poisson distribution
            empty_matrix[:,k] = poisson.rvs(mean.values)

        synthetic_counts = RnaCounts(data = empty_matrix, genes = self.rnacounts.genes, cells = self.rnacounts.cells)
        synthetic_pca_model = PCAModel(rnacounts=synthetic_counts, n_components=max_component)

        synthetic_pca_model.fit(scaler = self.rnacounts.median)
        variance_first_pc = synthetic_pca_model.pca_model.explained_variance_[0]
        variance_threshold  =  variance_first_pc*variance_multiplier

        # Now, we want to know how many PC components are above variance threshold
        original_pca_model = PCAModel(rnacounts=self.rnacounts, n_components=max_component)
        original_pca_model.fit()
        required_n_components = np.sum(original_pca_model.explained_variance >= variance_threshold)
        _LOGGER.info("The number of PC components estimates is {}".format(required_n_components))

        self.n_components = required_n_components

        return required_n_components

    def denoise(self):
        pass


