#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 15, 2020
# License: MIT License

# Clustering Algorithms

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

import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
import math

def dbscan(latentdf: pd.DataFrame, min_samples_frac: float = 0.01, max_dist_frac: float = 0.03,
    metric: str = 'euclidean', seed: int = 0, **kwargs):
    """

    Perform clustering on latent embedding of cellxgene matrix using DBSCAN algorithm.

    Parameters
    --------------

    latentdf: `pd.DataFrame`
        two dimension dataframe to be clustered

    min_samples_frac: `float`
        The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.

    max_dist_frac: `float`
        The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.

    metric: `str`
        Distance metric used in DBSCAN

    seed: `int`
        Seed for reproducibility

    algorithm: `str`
        Available algorithms from ["auto", "ball_tree", "kd_tree", "brute"]

    Returns
    ----------

    """
    value_range = np.ptp(latentdf.values, axis = 0)
    eps = max_dist_frac*math.sqrt(np.sum(np.power(value_range, 2.0)))
    min_samples =int(math.ceil(min_samples_frac * latentdf.shape[0]))

    algorithm = kwargs.get("algorithm", "brute")

    _LOGGER.info("DBSCAN clustering started with eps {} and min_samples {} ".format(eps, min_samples))

    dbscan = DBSCAN(algorithm=algorithm, min_samples =min_samples, eps = eps, metric = metric )
    cluster_labels = dbscan.fit_predict(latentdf.values)

    clusters = []
    for i, label in enumerate(cluster_labels):
        if label == -1:
            clusters.append('Outliers')
        else:
            clusters.append('Cluster {}'.format(label))
    
    cluster_df = pd.DataFrame(index = latentdf.index, data =clusters, columns = ['dbclusters'] )
    _LOGGER.info("Cluser labels saved in DataFrame with column name dbclusters")

    return cluster_df

def optics(latentdf: pd.DataFrame, min_samples_frac: float = 0.01, max_dist_frac: float = 0.03,
    metric: str = 'minkowski', seed: int = 0, **kwargs):
    """

    Perform clustering on latent embedding of cellxgene matrix using OPTICS algorithm.

    Parameters
    --------------

    latentdf: `pd.DataFrame`
        two dimension dataframe to be clustered

    min_samples_frac: `float`
        The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.

    max_dist_frac: `float`
        The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.

    metric: `str`
        Distance metric used in DBSCAN. Default is minkowski

         ["braycurtis", "canberra", "chebyshev", "correlation", "dice", "hamming", "jaccard", "kulsinski", "mahalanobis", "minkowski", "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean", "yule"]


    seed: `int`
        Seed for reproducibility

    algorithm: `str`
        Available algorithms from ["auto", "ball_tree", "kd_tree", "brute"]

    Returns
    ----------

    """
    value_range = np.ptp(latentdf.values, axis = 0)
    eps = max_dist_frac*math.sqrt(np.sum(np.power(value_range, 2.0)))
    min_samples =int(math.ceil(min_samples_frac * latentdf.shape[0]))

    algorithm = kwargs.get("algorithm", "brute")

    _LOGGER.info("OPTICS clustering started with eps {} and min_samples {} ".format(eps, min_samples))

    dbscan = OPTICS(algorithm=algorithm, min_samples =min_samples, eps = eps, metric = metric )
    cluster_labels = dbscan.fit_predict(latentdf.values)

    clusters = []
    for i, label in enumerate(cluster_labels):
        if label == -1:
            clusters.append('Outliers')
        else:
            clusters.append('Cluster {}'.format(label))
    
    cluster_df = pd.DataFrame(index = latentdf.index, data =clusters, columns = ['opticlusters'] )
    _LOGGER.info("Cluser labels saved in DataFrame with column name opticlusters")

    return cluster_df

    