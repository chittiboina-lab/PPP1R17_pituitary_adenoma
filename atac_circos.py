#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pyBigWig
import pycircos

mpl.use("agg")

Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

tumor_1mb = pyBigWig.open('tumor_corts_1mb_SF.bw')
normal_1mb = pyBigWig.open('normal_corts_1mb_SF.bw')

colors = ['#00c9b8', '#cef5f2', '#0ace5c', '#d2f6e1', '#1fcf0f', '#d6f6d3', '#75c800', '#e2f4c8', '#c5b500', '#eae398', '#e4a679', '#f5dccb', '#ec9ea5', '#f4c8cc', '#eb9ac9', '#f4c6e1', '#e498ea', '#f0c5f3', '#c4a6ed', '#ddccf5', '#a9aeee', '#cdd1f5', '#7dbae5', "#adb4b9"]
print(len(colors))
#Set chromosomes
lengths = {}
circle = Gcircle() 
with open("chr_start_end_hg38.txt") as f:
    color_count =-1
    for line in f:
        color_count += 1
        line   = line.rstrip().split("\t") 
        name   = line[0]
        length = int(line[-1])
        lengths[name] = length
        arc    = Garc(arc_id=name, label=name[3:], size=length, interspace=3, facecolor = colors[color_count], raxis_range=(510,560), label_visible=True)
        circle.add_garc(arc) 
        
circle.set_garcs()

vmin = 0
vmax = 20
for k, v in lengths.items():
    positions = np.arange(1, v)
    T_values = tumor_1mb.values(k, 1, v, numpy=True)
    N_values = normal_1mb.values(k, 1, v,numpy=True)
    circle.fillplot(k, data=N_values, positions=positions, 
               base_value=0,
               raxis_range=[800,1000], facecolor="royalblue", linewidth=0.01,
               rlim=[vmin, vmax]
               )
    circle.fillplot(k, data=T_values, positions=positions, 
                   base_value=0,
                   raxis_range=[600,800], 
                   facecolor="crimson",
                   linewidth=0.01,
                   rlim=[vmin, vmax]
                   )
    print('Plotted:', k)

circle.figure.savefig('circos_7-23/circos_0-20.png')
