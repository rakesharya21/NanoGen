#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 01:33:29 2018
Finding adsorbation sites
@author: rakesh
"""


from pyhull.delaunay import DelaunayTri as DT
from objects import ttd, atom
from numpy import linalg, array, dot


def unit(N):
    """
    Input: a numpy array object N
    Output: unit vector for N
    """
    return N/linalg.norm(N)


def gen_ads_sites(layer, ads):
    """
    Input:
        layer: an outer layer (atoms shouldn't be missing, extra is fine)
        ads: pymatgen atom object
    Output:
        list of all atom objects WITH ROOTS
    """

    # T package returns indices of atoms to be considered
    ttd_indices = DT([i.ic for i in layer.atom_list])

    tol_h = 0.01 # float(input('Enter minimum height of a tetrahedral formed on surface (0.01)\n'))
    tol_d = 3*layer.atom_list[0].radius # float(input('Enter maximum distance between two surface atoms (3*atomic_radius testing)\n'))

    # dict for all_ttd
    all_ttd = []

    # rc (return coord) returns the coords of concerned index
    rc = lambda x : layer.atom_list[x].ic

    # making all ttds in ttd objects
    for i in ttd_indices:

        ttd_coord = [rc(j) for j in i]
        print(ttd_coord)

        # adding ttd objects to all_ttd
        all_ttd.append(ttd(i, ttd_coord, tol_d, tol_h))

    # faces_indices will have list of lists of relative face indices
    faces_indices = []

    # checking if ttd holds for all requirements of making surface
    # ttd object KNOWS if it is legit and on surface
    for i in all_ttd:
        if i.is_surface(all_ttd):
            faces_indices.append[i.face]

    # deleting all_ttd and rc, as not needed now
    del all_ttd, rc

    ar = ads.atomic_radius
    lr = layer.atom_list[0].atomic_radius

    #  a list of sites [atom objects WITH ROOTS]
    sites = []

    # keeps track of previously considered atoms to avoid repitions
    roots = []
    
    def return_atom(x):
        # return atom indexed x from the atom_list of layer
        return layer.atom_list[x]
    
    # add ontop adsorbate sites
    # simply add on all atoms on surface (no repititions)
    for i in faces_indices:
        for j in i:
            root = j
            if root not in roots:
                roots.append(root)
                c = return_atom(j).ic
                cn = unit(c)
                d = (ar+lr)*cn
                new = c+d
                sites.append(atom(ads, new, root = [return_atom(j)]))

    # add bridge elements
    # initializing roots (not necessary)
    roots = []
    comb = [[0,1],[0,2],[1,2]]

    for F in faces_indices:
        for i in comb:

            # root is a set so order doesn't matter
            root = {F[j] for j in i}

            if root not in roots:

                # adding root to roots so that it is not repeated
                roots.append(root)

                # getting atoms to add a new bridge atom
                A = layer.atom_list[faces_indices[i[0]]].ic
                B = layer.atom_list[faces_indices[i[1]]].ic

                # midpoint
                C = (A + B)/2
                c = unit(C)
                d = ar + lr
                r = linalg.norm(A-B)

                # only if atom is too big to sit at midpoint
                if d > r/2:
                    new = C + c*pow((pow(d, 2)-pow(r/2, 2)), 0.5)

                # else same coordinates are used irrespective of size
                else:
                    new = C

                sites.append(atom(ads, new, root = [return_atom(x) for x in root]))

    # initializing roots (not necessary)
    roots = []

    # func checks if a face is square like (then no hcc/hcp site is added)
    def is_squarelike(F):
        C = [return_atom(x).ec for x in F]
        comb = [[0,1,2],[1,0,2],[2,1,0]]
        for i in comb:
            a = unit(C[i[1]] - C[i[0]])
            b = unit(C[i[2]] - C[i[0]])

            # if angle is more than 75
            if abs(dot(a,b)) < 0.258:
                return 1

            # else check other comb
            else:
                continue

        # if all angles less than 75, must be equilateral tri
        else: return 0

    # add fcc/hcp elements
    for F in faces_indices:
        root = set(F)
        if root not in roots:
            if is_squarelike(F):
                pass
            else:
                A = return_atom(F[0]).ic
                B = return_atom(F[1]).ic
                C = return_atom(F[2]).ic
                D = (A+B+C)/3
                d = unit(D)
                e = linalg.norm(A-B)
                new = D+d*(pow(2/3,0.5)*e)
                sites.append(atom(ads, new, c_tag = 'i', root = [return_atom(x) for x in root]))
    return sites