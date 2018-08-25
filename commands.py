# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 11:54:28 2018
Commands
@author: rakesh
"""

def get_unique_atoms(mol):
    """
    Input:
        mol: pymatgen class Molecule(species, coords)
        use tolerance as system set
    Output:
        eq_sets: dict with keys as unequivalent sites
    """
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer as pga
    #defining a symmetry object
    symm = pga(mol)
    u_sites = symm.get_equivalent_atoms().eq_sets
    print(u_sites)
    

def add_ads(layer, ads):
    """
    Input:
        An outer layer of atom objects as a list
    Processing:
        Gets sites from gen_ads_sites
        Gets symmetrically unequivalent sites form pymatgen.symmetry.analyzer
    Outputs:
        Each symmetrically unequivalent adsorbed config as a POSCAR
    """
    #add mads to be done later
    from ads import gen_ads_sites
    
    #getting all sites (atoms with roots) from gen_ads_sites
    sites = gen_ads_sites(layer, ads)
    species = [i.symbol for i in layer.atom_list]
    
    #Using ideal_coords to get unique sites 
    coords = [i.ic for i in layer.atom_list]
    
    mol = (species, coords)
    
    get_unique_atoms(mol)