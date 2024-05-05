#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris Brasnett
"""


from vermouth.gmx import read_itp, write_molecule_itp
from vermouth.forcefield import ForceField
import os
import argparse
from argparse import RawTextHelpFormatter

if __name__ == '__main__':

    msg = ('Write out simplified molecule topologies to use cg_bonds-v5.tcl\n\n'
           'output .itp files will have the same name as the input ones, but\n'
           'only consist of [ moleculetype ], [ atoms ], and [ bonds ] directives.\n'
           'We also write a vis.top with correct paths and names for the new files.\n'
           'In VMD, first source cg_bonds.tcl, then run:\n'
           '\tcg_bonds -top vis.top\n'
           'and now you should be able to see your system with the CG bonds.'
           )

    parser = argparse.ArgumentParser(description=msg, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-p", type = str, help  = 'input .top file used', default = 'topol.top')
    args = parser.parse_args()

    #get files in the present directory
    files = os.listdir(os.getcwd())

    #get the topology file
    top = os.getcwd() + f'/{args.p}'

    #read the topology file, find the header lines which are not core martini itps
    with open(top) as f:
        topol_lines = f.readlines()

    topol_core_itps = [a.split('"')[1] for a in topol_lines if "martini" in a]
    topol_rest = [a for a in topol_lines if "#include" not in a]
    topol_editable_itps = [a.split('"')[1] for a in topol_lines if "#include" in a and "martini" not in a]

    #for each molecule in the system, read in the itp
    d = {}
    for i,j in enumerate(topol_editable_itps):
        with open(j) as f:
            d[i] = f.readlines()

    #read the molecules into the forcefield
    ff = ForceField('martini3001')

    for i in d.keys():
        read_itp(d[i], ff)

    #iterate over the molecules to make visualisation topologies
    keep = ['bonds', 'constraints', 'pairs', 'virtual_sitesn',
            'virtual_sites2', 'virtual_sites3']
    for molname, block in ff.blocks.items():
        #delete the interactions which are not bonds
        for interaction_type in list(block.interactions):
            if interaction_type not in keep:
                    del block.interactions[interaction_type]
        # remove meta (ie. the #IFDEF FLEXIBLE) from the bonds
        for bond in block.interactions['bonds']:
            bond.meta.clear()

        # this should then keep any constraints which don't have IFDEF statements
        # eg. alpha helices are described by constraints without these.
        # however, the remove_interactions function doesn't work atm.
        # for bond in block.interactions['constraints']:
        #     if bond.meta:
        #         block.remove_interaction('constraints', bond.atoms)

        #rewrite pairs as bonds for visualisation
        for bond in block.interactions['pairs']:
            block.add_interaction('bonds', bond.atoms, bond.parameters[:2] + ['10000'])

        del block.interactions['pairs']

        #make bonds between virtual sites and each of the constructing atoms
        for vs_type in ['virtual_sitesn', 'virtual_sites2', 'virtual_sites3']:
            for vs in block.interactions[vs_type]:
                site = vs.atoms[0]
                constructors = vs.atoms[1:]
                for constructor in constructors:
                    # completely arbitrary parameters, the bond just needs to exist
                    block.add_interaction('bonds', [site, constructor],
                                          ['1', '1', '10000'])
            del block.interactions[vs_type]

        #write out the molecule with an amended name
        mol_out = block.to_molecule()
        mol_out.meta['moltype'] = molname+'_vis'

        header = [f'Visualisation topology for {molname}', 'NOT FOR SIMULATIONS']

        with open(molname+'_vis.itp', 'w') as outfile:
            write_molecule_itp(mol_out, outfile=outfile, header = header)

    #make a new header for the .top file to write the absolute paths
    new_topol_head = []
    for i in topol_core_itps:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    #get the new vis files to write for the topology header
    vis_files = [i for i in os.listdir(os.getcwd()) if '_vis' in i]
    for i in vis_files:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    #correct the [ molecules ] directive of the top file to correct for the new molecule names.
    topol_rest_vis = []
    original_mols = [i.split('_vis')[0] for i in vis_files]
    for i in topol_rest:
        if any(mol in i for mol in original_mols):
            topol_rest_vis.append(i.split()[0]+'_vis\t' + i.split()[1]+'\n')
        else:
            topol_rest_vis.append(i)

    #combine the sections of the vis.top and write it out.
    vis_topol = new_topol_head + topol_rest_vis

    with open('vis.top', 'w') as f:
        f.writelines(vis_topol)



