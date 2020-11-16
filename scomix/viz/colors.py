#!/usr/bin/env python
# coding: utf-8

# single-cell omics scomix

# Author : Rahul Bhadani
# Initial Date: Nov 15, 2020
# License: MIT License

# Color Definition

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

COLOR_LONG = ["#46f0f0", '#3a86ff', '#463F1A','#60492C', '#ff0022', '#208AAE', '#0D2149',\
              '#9ff9c1','#cb9173','#fdcff3',\
              '#2c80ff','#ffab2b','#fffb00', '#06b4ff','#23d160',\
              '#264653','#E9C46A',"#7900D7", "#A77500", "#6367A9", \
              "#A05837", "#772600", "#9B9700",'#E76F51',\
              '#B07BAC',"#7A87A1",'#5F7367','#472836',\
              "#00489C","#CB7E98", "#6F0062", "#EEC3FF",\
          '#97D8B2','#17BEBB','#C4A69D','#a3a271','#000000','#FFFFC7',\
          '#DB5461','#016FB9','#E56399','#D17B0F','#FFFC31','#772F1A',\
          '#808080','#000000', "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", \
          "#008941", "#006FA6", "#A30059","#FFDBE5", "#7A4900", "#0000A6", \
          "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87","#5A0007", \
          "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF","#4A3B53", \
          "#FF2F80","#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", \
          "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",\
          "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",\
          "#456D75", "#B77B68", "#788D66",\
          "#885578", "#FAD09F", "#1F819A", "#D157A0", "#BEC459", "#000080", "#0086ED", "#886F4C",\
          "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",\
          "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",\
          "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",\
          "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",\
          "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#324E72"]