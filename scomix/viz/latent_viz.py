#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 15, 2020
# License: MIT License

# Visualization Plots

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

from typing import KeysView, Union
from ..rna import RnaCounts
from ..latent import PCAModel
from .colors import COLOR_LONG
import logging

_LOGGER = logging.getLogger(__name__)

from sklearn.manifold import TSNE
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import plotly.offline as pyo
# Set notebook mode to work in offline
pyo.init_notebook_mode()

def tsne_plot(rnacounts, n_components = 30, perplexity = 30, seed = 0, latent_model = "PCA", exaggerated_tsne = False, random_order = False,  title= "TSNE Plot", cluster_labels = None, cluster_col_names = "Cell Types", **kwargs):
    """
    Make a TSNE plot

    Parameters
    -------------
    
    rnacounts: `RnaCounts`
        `RnaCounts` object

    n_components: `int`
        Number of components for TSNE plot

    perplexity: `float`
        Perplexity for TSNE plot

    latent_model : `str`
        What model to use for latent embedding for TSNE plot
        Options: "PCA" or a fully-qualified class name for the model

    seed: `int`
        Seed for reproducibility

    Attributes
    ------------

    """

    embedding = None
    init = "random"
    scores = None

    transpose = kwargs.get("transpose", False)
    if latent_model == "PCA":
        embedding = PCAModel(matrix = rnacounts, n_components = n_components)
        scores = embedding.fit(transpose = transpose)
    

    tsne_kwargs = {}
    if exaggerated_tsne:
        tsne_kwargs["early_exaggeration"] = 4
        tsne_kwargs["n_iter"] = 250
        init = "random"

    tsne_model = TSNE(perplexity = perplexity, random_state = seed, init = init, **tsne_kwargs)

    Z = tsne_model.fit_transform(scores.values)

    dim_labels = None
    if exaggerated_tsne:
        dim_labels = ['t-SNE* dimension %d' % (l+1) for l in range(Z.shape[1])]
    else:
        dim_labels = ['t-SNE dimension %d' % (l+1) for l in range(Z.shape[1])]

    tsne_scores = pd.DataFrame(data = Z, index =scores.index, columns = dim_labels)

    marker_size = kwargs.get("marker_size", 5)
    width = kwargs.get("width", 1200)
    height = kwargs.get("height", 800)
    font_size = kwargs.get("font_size", 20)

    fig = plot_cells(df = tsne_scores, cluster_labels = cluster_labels,title  = title, cluster_col_names = cluster_col_names, width = width, height = height, marker_size = marker_size, font_size = font_size)

    return fig, tsne_scores


def plot_cells(
    df: pd.DataFrame,
    title: str,
    seed: int = 0, 
    cluster_labels = None,
    **kwargs) -> go.Figure:
    """
    Two-dimensional plot of cells features as clusters

    Parameters
    -------------
    df: `pd.DataFrame`
        Dataframe of latent embedding to plot on 2-dimensional plot as a scatter plot
    
    title: `str`
        Title of the plot

    cluster_labels: `pd.DataFrame`
        Cluster labels as dataframes. Index of the dataframe must match the index of `df` above

    cluster_col_names: `str`
        Name of cluster columns

    
    """


    marker_size = kwargs.get("marker_size", 5)
    width = kwargs.get("width", 1200)
    height = kwargs.get("height", 800)
    font_size = kwargs.get("font_size", 20)
    
    show_plot = kwargs.get("show_plot", False)
    xaxis_label = df.columns[0]
    yaxis_label = df.columns[1]
    fig = go.Figure()

    if cluster_labels is None:
        # TSNE component for x-axis
        x = df[:,0].values
        # TSNE components for y-axis
        y = df[:,1].values

        cell_id = df.index.tolist()
        scatter_trace = go.Scatter(x = x, y = y, mode = 'markers',
                                text = cell_id, marker = dict(size = marker_size,
                                color = "#4d194d"
                                ))
        fig.add_trace(scatter_trace)
        fig.update_layout(title = title,
            xaxis_title = xaxis_label,
            yaxis_title= yaxis_label
            )

    else:
        index_df = df.index
        index_cluster =cluster_labels.index
        
        if not index_df.equals(index_cluster):
            _LOGGER.error("Cell id/Index of Cluster dataframe is not the  same as of the latent dataframe.")
            raise ValueError("Index mismatch between cluster label and latent dataframe")
        
        cluster_col_names = kwargs.get("cluster_col_names", "cluster_label")

        if cluster_col_names not in cluster_labels.columns:
            _LOGGER.error("Cluster label column not present in cluster dataframe")
            raise ValueError("Cluster label column must be one of the columns of cluster_labels dataframe")

        df_with_cluster = pd.concat([df, cluster_labels], axis= 1)

        labels = df_with_cluster[cluster_col_names].unique()

        for label_index, onelabel in enumerate(labels):
            df_onelabel = df_with_cluster[df_with_cluster[cluster_col_names] == onelabel]
            cell_id = df_onelabel.index.tolist()
            scatter_trace = go.Scatter(x = df_onelabel[df_with_cluster.columns[0]], 
                                                y = df_onelabel[df_with_cluster.columns[1]],
                                                mode = 'markers', text = cell_id,
                                                name = onelabel, 
                                                marker = dict(size = marker_size,
                                                color = COLOR_LONG[label_index],
                                                opacity = 0.7
                                    ))
                
            fig.add_trace(scatter_trace)
        
        fig.update_layout(title = title,
            height = height,
            width = width,
            xaxis_title = xaxis_label,
            yaxis_title= yaxis_label,
            legend_title = cluster_col_names,
            font_size = font_size
            )

    return fig
    


        




    return