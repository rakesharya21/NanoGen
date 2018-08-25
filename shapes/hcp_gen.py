#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 03:29:13 2018
hcp
@author: rakesh
"""

import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

tag = 'hcp'

coords = []

coords.extend([np.array((0,0,0)), np.array((2,0,0))])

theta = np.radians(60)
R = np.array(((np.cos(theta), np.sin(theta), 0), (-np.sin(theta), np.cos(theta), 0), (0,0,1)))

#making hexagon
print(R, coords[-1])
for i in range(5):
    coords.append(np.matmul(coords[-1], R))

theta = np.radians(120)
R = np.array(((np.cos(theta), np.sin(theta), 0), (-np.sin(theta), np.cos(theta), 0), (0,0,1)))

#adding atoms on top (alternate)
coords.append((coords[0]+coords[1]+coords[2])/3)
coords[-1] = coords[-1] + np.array((0,0,np.sqrt(2/3)))

for i in range(2):
    coords.append(np.matmul(coords[-1], R))

#adding atoms at bottom (inverse of top)
coords.append((coords[0]+coords[2]+coords[3])/3)
coords[-1] = coords[-1] - np.array((0,0,np.sqrt(2/3)))

for i in range(2):
    coords.append(np.matmul(coords[-1], R))

#writing to file
with open(os.path.join(dir_path, 'hcp'), 'w') as f:
    f.write(tag + '\n')
    for i in coords:
        f.write(' '.join([str(j) for j in list(i)]) + '\n')