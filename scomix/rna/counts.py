#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 14, 2020
# License: MIT License

# This code is inspired from Monet's Expression Matrix

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

__author__ = 'Rahul Bhadani'
__email__  = 'rahulbhadani@email.arizona.edu'

import modin.pandas as pd
import numpy as np
import logging
import scipy
import scanpy

_LOGGER= logging.getLogger(__name__)

class counts(pd.DataFrame):
    """

    Class to implement count matrix fo single-cell RNA transcriptomes

    Each column of the dataframe represents gene and each row represents a cell sample

    """
    def __init(self, *args, **kwargs):
        """
        Parameters
        ------------
        genes:
            Specify the name of genes as list, numpy array or pandas series

        cells:
            Specify the cell identified or bar code as list, numpy array or pandas series
            
        """

        pd.DataFrame.__init__(self, *args, **kwargs)
        
        genes = kwargs.get("genes", None)
        cells = kwargs.get("cells", None)
        if genes is not None:
            self.columns = genes
        if cells is not None:
            self.index = cells

        self.index.name = 'Cells'
        self.columns.names = 'Genes'
        
    
    @property
    def n_genes(self):
        """
        Number of genes
        """
        return self.shape[1]

    @property
    def n_cells(self):
        """
        Number of cells
        """
        return self.shape[0]

    @property
    def genes(self):
        """
        `list`
            Return a list of genes
        """
        return self.columns.values

    @property
    def X(self):
        """
        CellxGene matrix
        """
        return self.values

    
    def transform(self, func, inplace= False,*args, **kwargs):
        """
        Apply a transformation to the transcriptomics
        
        func:`str`
             Type of transformation. Current transformations include "freeman-tuket", "anscombe", "log" or any expression in string format, e.g. "np.sqrt(x) + np.sin(x)", x must be used as a variable in the expression.

        """

        new_counts = None
        if func.lower() == "anscombe":
            new_counts = 2* np.sqrt(self + 3.0/8.0)
        elif func.lower() == "freeman-tukey":
            new_counts = np.sqrt(self) + np.sqrt(self+1)

        elif func.lower() == "log":
            new_counts = np.log(self + 1)

        else:
            try:
                func_self = func.replace("x", "self")
                new_counts = eval(func_self)
            except SyntaxError:
                raise SyntaxError("Syntax is not valid for transformation function {} ".format(func))

        if inplace:
            self._update_inplace(new_counts)
            new_counts = self

        return new_counts

    @classmethod
    def from_anndata(cls, adata, dtype=None):
        """
        Import trascriptmomics data from scanpy's `AnnData` object

        Parameters
        -------------
        adata: `scanpy.AnnData`
            scanpys AnnData object

        dtype:  `str`
            Data type to cast, e.g. np.float32

        """
        transcriptomics= adata.X.copy(deep=True)

        if dtype is not None:
            transcriptomics = transcriptomics.astype(dtype)

        if scipy.sparse.issparse(transcriptomics):
                _LOGGER.info("Count matrix is sparse, covert it to proper matrix")
                transcriptomics  = scipy.sparse.csr_matrix.toarray(transcriptomics)

        genes = adata.var_names.copy(deep = True)
        genes.name = 'Genes'

        cells = adata.obs_names.copy(deep = True)
        cells.name = 'Cells'

        counts = cls(data=transcriptomics, genes = genes, cells = cells)

        return counts


