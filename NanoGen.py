#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 03:19:23 2018
NanoGen
Generates desired NanoClusters
Adds adsorbates
Powerful Object Oriented Program
@author: rakesh
"""

from cast import read_poscar, read_instructions
from ads import *

"""
Input types:
    POSCAR: a valid poscar file to be read by pymatgen Structure module
    instructions.txt:
        order of file is:
            shape <shape_name-1> <atom>'
            layer1 <atom>
            layer2 <atom>
            ...
            ads <atom> (optional)
            shape <shape_name-2> <atom>'
            layer <atom>
            layer <atom>
            ...
            ads <atom> (optional)
"""

print('Welcome to NanoGen :)')
cmd = input('Read from POSCAR <p> or read from instructions <i>? ')

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

#initializing global layers list and atom dict
layers = []
atom_dict = {}

while True:
    if cmd == 'p':
        read_poscar(layers, atom_dict)
        cmd2 = input("Which molecule to adsorb? ")
        break
    elif cmd == 'i':
        read_instructions(layers, atom_dict)
        break
    else:
        print('command not accepted :/ try again!')
        continue