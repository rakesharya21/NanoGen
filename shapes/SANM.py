#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 03:33:45 2018
nanowire
@author: rakesh
"""

tag = 'sanm'

from numpy import array as ar
from numpy import sqrt, cos, sin

coords = []

coords.append(ar([0,0,0]), ar([1,1,1]))

l = int(input('Enter length of wire: '))

line = []
ring = []

#generating a single atom wire
for i in range(1, l):
    line.append(coords[1] + ar([l,l,l]))

#generating a single atom ring
ring.append(-0.5*coords[0])
ring.append(0.5*coords[1])

theta


with open('sanm_line','w') as f:
    s.write(tag+'\n')
    s.write('\t0\nNone\n')
    s.write('1\n')



print('nanowire generated')