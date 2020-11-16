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

import pandas as pd
import numpy as np
import tarfile
import io
import gzip
import ntpath
import logging
from pandas.core.arrays import sparse
import scipy
import scipy.io
import csv
from datetime import datetime

import scanpy

_LOGGER= logging.getLogger(__name__)

class RnaCounts:
    """

    Class to implement count matrix fo single-cell RNA transcriptomes

    Each column of the dataframe represents gene and each row represents a cell sample

    Parameters
    --------------

    data:
        specify data to be read
    genes:
        Specify the name of genes as list, numpy array or pandas series

    cells:
        Specify the cell identified or bar code as list, numpy array or pandas series

    Attributes
    ------------
    cellxgene: `pd.DataFrame`
        A matrix/DataFrame with transcriptome counts

    X: `np.array`
        Transcritomic values in numpy array format without gene or cell labels

    genes:
        A list available genes in the transcritome counts

    cells:
        Cell identifier

    sequence: `int`
        Sequence number to identify which operation was performed in what order

    records:`dict`
        A dictionary to keep track of workflow

    media:`float`
        Median count value by aggregating counts for each cell
        We first sum counts for each cells and then take median

    PC: `pd.DataFrame`
        Prinicipal Component of cellxgene matrix. Call `PCAModel` for calculating principal components
    
    raw: `pd.DataFrame`
        original cellxgene datframe; None if cellxgene not modified.

    """
    def __init__(self,  data = None, genes = None, cells = None, *args, **kwargs):
        """
            
        """
        
        self.cellxgene = pd.DataFrame(data)

        if genes is not None:
            self.cellxgene.columns = genes
        
        if cells is not None:
            self.cellxgene.index = cells

        self.cellxgene.index.name = 'Cells'
        self.cellxgene.columns.name = 'Genes'

        self.sequence = 1
        self.records = {}

        self.PC = None
        self.raw = None
        
    
    @property
    def n_genes(self):
        """
        Number of genes
        """
        return self.cellxgene.shape[1]

    @property
    def n_cells(self):
        """
        Number of cells
        """
        return self.cellxgene.shape[0]

    @property
    def genes(self):
        """
        `list`
            Return a list of genes
        """
        return self.cellxgene.index


    @property
    def cells(self):
        """
        `list`
            Return a list of cells
        """
        return self.cellxgene.columns
        

    @property
    def X(self):
        """
        CellxGene matrix
        """
        return self.cellxgene.values

    
    def transform(self, func, inplace= False, *args, **kwargs):
        """
        Apply a transformation to the transcriptomics
        
        func:`str`
             Type of transformation. Current transformations include "freeman-tukey", "anscombe", "log" or any expression in string format, e.g. "np.sqrt(x) + np.sin(x)", x must be used as a variable in the expression.

        """

        new_counts = None
        if func.lower() == "anscombe":
            new_counts = 2* np.sqrt(self.cellxgene + 3.0/8.0)
        elif func.lower() == "freeman-tukey":
            new_counts = np.sqrt(self.cellxgene) + np.sqrt(self.cellxgene+1)

        elif func.lower() == "log":
            new_counts = np.log(self.cellxgene + 1)

        else:
            try:
                func_self = func.replace("x", "self.cellxgene")
                new_counts = eval(func_self)
            except SyntaxError:
                raise SyntaxError("Syntax is not valid for transformation function {} ".format(func))

        if inplace:
            if self.raw is None:
                self.raw = self.cellxgene

            self.cellxgene = new_counts

        self.records["transform"] = {}
        self.records["transform"]["sequence"] = self.sequence 
        self.sequence  =  self.sequence  + 1
        self.records["transform"]["timestamp"] =datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S.%f')
        self.records["transform"]["type"] = func
        self.records["transform"]["inplace"] = inplace

        return new_counts

    @property
    def median(self):
        """
        The median count value of  counts from cellxgene.
        Median count value is calculated by aggregating counts for each cell. We first sum counts for each cells and then take median
        """
        return float(self.cellxgene.sum(axis = 1).median())

    def scale(self, scale_value, inplace= False):
        """
        Scale all transcriptome counts to `upper_limit`

        Parameters
        -------------
        scale_value: `float`
            Scale value


        """

        if self.records and 'scale' in self.records.keys:
            _LOGGER.info("cellxgene count matrix has already been scaled. Not re-scaling. Returning existinf cellxgene")
            return self.cellxgene
        
        count_sum_by_cell = self.cellxgene.sum(axis = 1)

        scaled_transcriptome = (scale_value/count_sum_by_cell)*self.cellxgene.T

        if inplace:
            if self.raw is None:
                self.raw = self.cellxgene

            self.cellxgene = scaled_transcriptome

        self.records["scale"] = {}
        self.records["scale"]["sequence"] = self.sequence 
        self.sequence  =  self.sequence  + 1
        self.records["scale"]["timestamp"] =datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S.%f')
        self.records["scale"]["scale_value"] = scale_value
        self.records["scale"]["inplace"] = inplace
        
        return scaled_transcriptome

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

        Returns
        ----------
        `RnaCounts`
        Returns `RnaCounts` object

        """
        transcriptomics= adata.X.copy()

        if dtype is not None:
            transcriptomics = transcriptomics.astype(dtype)

        if scipy.sparse.issparse(transcriptomics):
                _LOGGER.info("Count matrix is sparse, covert it to proper matrix")
                transcriptomics  = scipy.sparse.csr_matrix.toarray(transcriptomics)

        genes = adata.var_names.copy()
        genes.name = 'Genes'

        cells = adata.obs_names.copy()
        cells.name = 'Cells'

        counts = cls(data=transcriptomics, genes = genes, cells = cells)

        return counts

    @classmethod
    def from_npz(cls, npz):
        """
        Import data from numpy compressed file format
        

        Parameters
        -------------

        npz: `str`
            File path of numpy compressed data format

        Returns
        ----------
        `RnaCounts`
            Returns `RnaCounts` object

        """

        data = np.load(npz, allow_pickle=True)
        genes = data['genes']
        cells = data['cells']
        data = np.array(data['matrix'], order='F', copy=False)
        matrix = cls(data=data, genes = genes, cells = cells  )

        _LOGGER.info("Loaded expression matrix from {}".format(npz))

        return matrix

    @classmethod
    def from_tsv(cls, tsv, sep='\t', encoding='utf-8', index_col =0, header = 0, **kwargs):
        '''
        Import data from tsv formatted file

        Parameters
        -------------
        tsv: `str`
            Full file path were tsv file is located

        sep = `str`
            Value separators

        encoding = `str`
            File encoding

        index_col = `int`
            Which column is index

        header = `int`
            Specify header

        '''

        data = pd.read_csv(tsv, encoding=encoding, index_col=index_col, header=header,  sep = sep, **kwargs)

        data = data.astype('uint16')


        # Check if column is gene or row is genes
        # if rows are genes then transpose it
        if data.index.name.lower() == "genes":
            data = data.T
        

        matrix = cls(data = data)

        _LOGGER.info("Loaded expression matrix from {}".format(tsv))

        return matrix


    @classmethod
    def from_csv(cls, csv, sep=',', encoding='utf-8', index_col =0, header = 0, **kwargs):
        '''
        Import data from tsv formatted file

        Parameters
        -------------
        tsv: `str`
            Full file path were tsv file is located

        sep = `str`
            Value separators

        encoding = `str`
            File encoding

        index_col = `int`
            Which column is index

        header = `int`
            Specify header

        '''

        data = pd.read_csv(csv, encoding=encoding, index_col=index_col, header=header, sep = sep, **kwargs)

        data = data.astype('uint16')

        # Check if column is gene or row is genes
        # if rows are genes then transpose it
        if data.index.name.lower() == "genes":
            data = data.T

        matrix = cls(data = data)

        _LOGGER.info("Loaded expression matrix from {}".format(csv))

        return matrix

    @classmethod
    def from_10x(cls, path, use_ensemble_ids = False, dtype=np.uint16, skip_gene_chars=0):
        """
        Import data from 10x CellRanger tarball

        Parameters
        -------------
        path: `str`
            Full file path of the tar ball. File extension should end .tar.gz or zip

        use_ensemble_ids: `bool`

        dtype:
            Data types from numpy

        """

        file_extension = path.split(sep = '.')[-1]
        file_extension2 = path.split(sep = '.')[-2]

        data = None
        genes = None
        cells = None

        if file_extension2 + '.'  + file_extension== 'tar.gz':
            with tarfile.open(path, mode = 'r:gz') as tar:
                for tarinfo in tar.getmembers():
                    if "matrix.mtx" in tarinfo.name:
                        if tarinfo.size and tarinfo.isreg()> 0:

                            ti = tar.getmember(tarinfo.name)
                            with gzip.open(tar.extractfile(ti)) as tarmember:
                                sparse_data = scipy.io.mmread(tarmember)
                                sparse_data = sparse_data.astype(dtype)
                                data = sparse_data.todense(order = 'F')
                    
                    if ("features" in tarinfo.name) or  ("genes" in tarinfo.name):
                        if tarinfo.size and tarinfo.isreg()> 0:

                            ti = tar.getmember(tarinfo.name)
                            with gzip.open(tar.extractfile(ti)) as tarmember:
                                wrapper = io.TextIOWrapper(tarmember, encoding = "ascii")
                                i = 1
                                if use_ensemble_ids:
                                    i = 0
                                genes = [row[i][skip_gene_chars:] for row in csv.reader(wrapper, delimiter = '\t')]

                    if ("cells" in tarinfo.name) or  ("barcodes" in tarinfo.name):
                        if tarinfo.size and tarinfo.isreg()> 0:
                            ti = tar.getmember(tarinfo.name)
                            with gzip.open(tar.extractfile(ti)) as tarmember:
                                wrapper = io.TextIOWrapper(tarmember, encoding = "ascii")
                                cells = [row[0] for row in csv.reader(wrapper, delimiter = '\t')]


        print(data.shape)
        print(len(genes))
        print(len(cells))
        if data.shape[0] != len(genes):
            raise ValueError('Number of genes does not match! Something is wrong with the data.') 
        if data.shape[1] != len(cells):
            raise ValueError('Number of cells does not match! Something is wrong with the data.')
        
        data = pd.DataFrame(data.T)
        matrix = cls(data=data, genes = genes, cells = cells  )
        _LOGGER.info("Loaded expression matrix from  10x data located at {}".format(path))

        return matrix