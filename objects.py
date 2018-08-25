#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 00:54:28 2018
Objects
@author: rakesh
"""

class layer():
    def __init__ (self, atom_list):
        """
        Input:
            atom_list: a list of atom objects in the layer
        """
        self.atom_list = atom_list
        self.num = len(atom_list)

class ttd():
    def __init__ (self, indices, coords, tol_d, tol_h):
        """
        Input:
            self
            indices: RELATIVE indices of atoms in the layer as a list
            coords: coords of atoms
            tol_d: tolerance for length of edges
            tol_h: tolerance for height of ttd
        """
        self.index_dict = {}
        for i in indices:
            self.index_dict[i] = coords[i]
        self.tol_h = tol_h
        self.tol_d = tol_d

    def is_surface(self, all_ttd):
        """
        checks of a ttd is on surface
        checks if a ttd is legit first

        creates self.face (list of relative indices) if ttd is_surface
        """

        #is_legit takes care of distance tolerance of ttd
        if self.is_legit():

            #this code only compares indices to find surface
            #index list is a list of current ttd considered
            index_list = [i for i in self.index_dict]
            comb = [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]

            #y returns the xth index from index list
            y = lambda x: index_list[x]

            for i in comb:
                face = {y(j) for j in i}
                for j in all_ttd:

                    #if face shares all its indices with some other ttd
                    if len(i - face) == 1 and set(index_list) != j:
                        #then it is not a surface face
                        return 0
                    else:
                        #it may be a surface
                        continue

                else:
                    #making face as a list(probably set is better)
                    self.face = sorted(face)

                    #tells the calling function that ttd has a face
                    return 1
            else: print('this ttd is neither face nor not face :/')
        else: print('illegit ttd detected, not adding to list')

    #distance Plane and point
    def dPp(base, pivot):
        from numpy import array as ar
        from numpy import dot, cross

        """
        Input:
            base: a list of three lists of coordinates
            pivot: list coordinates of a point
        Output:
            absolute distance between base and pivot
        """
        def unit(N):
            from numpy import linalg
            """
            Input: a numpy array object N
            Output: unit vector for N
            """
            return N/linalg.norm(N)

        a,b,c = (ar(i) for i in base)
        p = ar(pivot)

        #relative vectors
        r1 = b-a
        r2 = c-a
        r3 = p-a

        #calulating distance (dot product with surface base unit vector)
        n = unit(cross(r1, r2))
        d = abs(dot(r3, n))

        return d

    #distance point and point
    def dpp(a, b):
        from numpy import linalg
        from numpy import array as ar
        a, b = ar(a), ar(b)
        return linalg.norm(a-b)

    def is_legit(self):
        """
        Calculates all conditions for a legit ttd
        tol_h for height
        tol_d for length of edges
        """
        td = self.tol_d
        th = self.tol_h
        I = {0,1,2,3}
        comb = [{1,2,3}, {0,1,2}, {0,1,3}]
        a = True
        for i in comb:
            p = list(I-i)[0]
            cl = list(i)
            index_list = [i for i in self.index_dict]

            #noting the coordinates of base and pivot selected to check distance
            base = [self.index_dict[index_list[i]] for i in cl]
            pivot = self.index_dict[index_list[p]]

            #checking for distance
            if self.dPp(base, pivot) > th: a = False
            else: continue

        if a == False: return 0
        else:
            comb = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
            for i in comb:

                #noting coordinates of two points to check distance
                a = self.index_dict[index_list[i[0]]]
                b = self.index_dict[index_list[i[1]]]

                #checking distance between points
                if self.dpp(a, b) > td: return 0
                else: continue
            else: return 1



class atom():
    def __init__ (self, mg_atom, coords, c_tag = 'i', num = -1, root = None):
        """
        Input:
            mg_atom: a pymatgen Element object
            coords: coordinates in cartesian (numpy array)
            c_tag: weather coords are ideal or exact
            num: global number of atom
            root: which atoms does atom correspond to in adsorption structure, list of atom objects
        """

        self.root = root
        self.mg_atom = mg_atom
        self.symbol = mg_atom.symbol
        self.radius = mg_atom.atomic_radius
        if c_tag == 'b':
            ideal_coords = coords
            exact_coords = coords
        elif c_tag == 'i':
            ideal_coords = coords
            exact_coords = self.calculate_ec()
        else:
            print('c_tag not recognised for atom:', self.num)

        self.ec = exact_coords
        self.ic = ideal_coords
        self.num = num


    def calculate_ec(self):
            root = self.root

            ec = []

            rl = len(root)
            if rl == 1:
                #ontop element
                pass
            elif rl == 2:
                #bridge
                pass
            elif rl == 3:
                #equilateral triangle
                pass

            return None