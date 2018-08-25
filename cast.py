#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 01:11:32 2018
Shape related operations
@author: rakesh
"""

def read_poscar(layers, atom_dict):
    from pymatgen import Structure, Element
    from objects import atom

    name = (input('Enter Name of POSCAR file in the same directory: '))
    poscar = Structure.from_file(name)
    for i in poscar:
        mg_atom = Element(str(i.specie))
        atom_dict[len(atom_dict)] = atom(mg_atom, list(i.coords), num = len(atom_dict), c_tag = 'b')
    layers.append([atom_dict[i] for i in atom_dict])
    print('Read POSCAR as a single layer')
    cmd = input('add layer <l> or add adsorbate <a>? ')
    while cmd == 'l':
        sym = input("Which atom to deposit? ")
        add_layer(Element(sym))
        cmd = input('add layer <l> or add adsorbate <a>? ')

    if cmd == 'a':
        add_ads(layers[-1], Element(sym))

def read_instructions(layers, atom_dict):
    import pymatgen as mg
    f = open('instructions.txt', 'r')
    ctr = 0
    for line in f:
        #Each command is sent a pymatgen object
        #commands are sent last layer or all layers as needed
        if line.startswith('shape'):
            shape = line.split()[1]
            shape_atom = mg.Elemet(line.split()[2])
            read_shape(shape, shape_atom, layers, atom_dict)
        elif line.startswith('layer'):
            add_layer(layers[-1], mg.Element(line.split()[1]))
        elif line.startswith('ads'):
            add_ads(layers, mg.Element(line.split()[1]))
        elif line.startswith('mads'):
            add_mads(layers, line.split[1])
        elif line.startswith('done'):
            create_poscar(layers, str(ctr))
        else: print('Could not read instructions file')

def read_shape(shape, shape_atom, layers, atom_dict):
    import numpy as np
    import pymatgen as mg
    import os
    from objects import atom

    dir_path = os.path.dirname(os.path.realpath(__file__))

    #opening shapes file
    with open(os.path.join(dir_path, 'shapes/' + shape)) as f:
        tag = f.readline()
        try: c = int(tag.split()[1])
        except IndexError: c = 0
        mg_atom = mg.Element(shape_atom)

        #initializing layers
        for j in range(c):
            coords = np.array(map(float, f.readline().split()))*mg_atom.atomic_radius
            atom_dict[len(atom_dict)] = atom(mg_atom, coords, num = len(atom_dict))

        layers.append([atom_dict[i] for i in range(c)])
        print('Core layer created, atoms:', len(layers[0]))
        for line in f:
            coords = np.array(map(float, f.readline().split()))*mg_atom.atomic_radius
            atom_dict[len(atom_dict)] = atom(mg_atom, coords, num = len(atom_dict))

        layers.append([atom_dict[i] for i in range(c, len(atom_dict))])
        print('Shell layer created, atoms:', len(layers[-1]))