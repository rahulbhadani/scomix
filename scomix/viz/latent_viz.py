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

from typing import Union
from ..rna import RnaCounts
from ..latent import PCAModel
from .colors import COLOR_LONG

from sklearn.manifold import TSNE
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import plotly.offline as pyo
# Set notebook mode to work in offline
pyo.init_notebook_mode()

def tsne_plot(rnacounts, n_components = 30, perplexity = 30, seed = 0, latent_model = "PCA", exaggerated_tsne = False, random_order = False, 
title= "TSNE Plot", **kwargs):
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
        dim_labels = ['t-SNE* dim. %d' % (l+1)
                      for l in range(Z.shape[1])]
    else:
        dim_labels = ['t-SNE dim. %d' % (l+1)
                      for l in range(Z.shape[1])]

    tsne_scores = pd.DataFrame(data = Z, index =scores.index, columns = dim_labels)

    fig = plot_cells(df = tsne_scores, title  = title)

    return fig, tsne_scores


def plot_cells(
    df: pd.DataFrame,
    title: str,
    cell_labels: pd.Series = None,
    seed: int = 0, 
    cluster_labels = None,
    **kwargs) -> go.Figure:
    
    width = 800
    height = 600
    marker_color = 'navy'
    marker_size  = 4
    show_scale = True
    colorscale = 'RdBu'
    reversescale = True
    opacity  = 0.75
    cluster_colors = {}
    n_cells = df.shape[0]
    
    cell_labels = pd.Series(index = df.index, data =["Cells"]*n_cells)

    # cluster order 
    vc = cell_labels.value_counts()
    cluster_order = vc.index.tolist()

    for i, cluster in enumerate(cluster_order):
        if cluster in cluster_colors:
            continue
        try:
            cluster_colors[cluster] = COLOR_LONG[i]
        except IndexError:
            cluster_colors[cluster] = None

    if cluster_labels is None:
        cluster_labels = dict([cluster, cluster] for cluster in cluster_order)

    trace_list =[]
    for cluster in cluster_order:
        label = cluster_labels[cluster]
        indices = (cell_labels== cluster)
        color = cluster_colors[cluster]

        select_df = df.loc[indices]
        x = select_df.iloc[:,0].values
        y = select_df.iloc[:,1].values
        
        text = select_df.index.tolist()
        
        colorbar = dict(
           len = 0.5,
           outlinecolor = 'red',
           outlinewidth = 1,
           ticklen = 6,
           ticks ='outside' ,
           titleside = 'top',
           separatethousands = True)
           
        trace = go.Scatter(x = x, y = y, text = text, name = label,  mode = 'markers', marker = dict(
           size = marker_size,
           color = color,
           showscale=show_scale,
           colorbar=colorbar,
           colorscale=colorscale,
           opacity = opacity))
        trace_list.append(trace)

    xaxis_label = df.columns[0]
    yaxis_label = df.columns[1]

    layout = go.Layout(
        width = width, height = height,
        xaxis = dict(linecolor = 'black', title = xaxis_label),
        yaxis = dict(linecolor = 'black', title = yaxis_label),
        title = title
    )

    fig = go.Figure(data = trace_list, layout = layout)

    return fig