from rdkit import Chem
from math import log2
from bottchscore.symmetry_tools import *

import logging
logger = logging.getLogger()

def calc_Di(atom):
    """ Calculate chemically equivalent neighbor term  of an atom """
    topo_set = set(x.GetIntProp('topo_id') for x in atom.GetNeighbors() if x.GetAtomicNum() != 1)
    di = len(topo_set)
    atom.SetIntProp('Di',di)
    return di

def calc_Vi(atom):
    """ Calculate the valence term of an atom """
    #TODO I think this conflicts with the definitions in the paper, which describes Vi as 'of a neutral element at position i'
    valence = atom.GetTotalValence()
    charge = atom.GetFormalCharge()
    vi = 8 - valence + charge
    atom.SetIntProp('Vi', vi)
    return vi

def calc_atomic_Si(atom):
    """ Calculate the chirality term of an atom """
    #TODO potential chiral centers are not accounted for here, and I'm not sure how we could (science-wise)
    #TODO reintegrate support for E/Z isomerism
    if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        si = 1
    else:
        si=2
    atom.SetIntProp('Si', si)
    return si

def calc_Bi(atom):
    """ Calculate the bond term of an atom """
    bi = 0
    for bond in atom.GetBonds():
        other_atom = bond.GetOtherAtom(atom)
        if other_atom.GetAtomicNum() == 1:
            continue
        bi += bond.GetBondTypeAsDouble()
    atom.SetProp('Bi',str(bi))
    return bi

def calc_Ei(atom):
    """ Calculate the elemental composition term of an atom """
    #TODO does not currently account for isotopic information as described in the paper
    element_list = [x.GetAtomicNum() for x in atom.GetNeighbors() if x.GetAtomicNum() != 1]
    element_list.append(atom.GetAtomicNum())
    element_set = set(element_list)
    ei = len(element_set)
    atom.SetIntProp('Ei', ei)
    return ei

def calc_local_terms(atom):
    di = calc_Di(atom)
    vi = calc_Vi(atom)
    bi = calc_Bi(atom)
    ei = calc_Ei(atom)
    logging.debug(f'Local atomic terms for {atom.GetSymbol()}{atom.GetIdx()}: {di=}, {ei=}, {vi=}, {bi=}')
    if bi == 0 or vi == 0: # Avoid log of 0 for reasons
        logging.debug(f'Atomic contribution for {atom.GetSymbol()}{atom.GetIdx()} is 0 by definition')
        partial_acm = 0
    else:
        partial_acm =di*ei*log2(bi*vi)

    if atom.GetIntProp('num_symmetry') > 0:
        partial_acm = 0.5*partial_acm
    
    logging.debug(f'Local atomic contribution for {atom.GetSymbol()}{atom.GetIdx()}: {partial_acm}')
    atom.SetDoubleProp('partial_acm', partial_acm)

def increment_atomic_Si(bond):
    """ Determines which atom has lower local Cm to increase Si"""
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    if atom1.GetDoubleProp('partial_acm') < atom2.GetDoubleProp('partial_acm'):
        atom1.SetIntProp('Si',atom1.GetIntProp('Si')+1)
    else:
        atom2.SetIntProp('Si',atom2.GetIntProp('Si')+1)

def calc_bond_Si(bond):
    """ Find bonds with labeled stereochemistry and increment lower atomic Si"""
    if bond.HasProp('Atropisomeric'):
        if bond.GetBoolProp('Atropisomeric'):
            increment_atomic_Si(bond)
    if bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE:
        #TODO this catches STEREOANY and I'm not sure if it should
        increment_atomic_Si(bond)

def distribute_chirality(mol):
    """ Assign si values to local (atomic) and nonlocal (bond) chirality """
    #TODO Curently only handles axial chirality if it embeds info in a bond (atropisomerism) or an atom (allene)
    #TODO Other axial chirality (e.g. helical chirality in Ru(bpy)3 or helicene) makes me nauseous
    #TODO Allene chiral flags exist in rdkit but are not parsed from smiles
    #TODO Also generally unsuited for metal centers
    #TODO Does not currently handle cases where a chiral center could exist but is unspecified
    for atom in mol.GetAtoms():
        calc_atomic_Si(atom)
    
    for bond in mol.GetBonds():
        calc_bond_Si(bond)   


def calc_cm(mol):
    """ Calculate the atomic contribution to the Bottcher Score """
    cm = 0
    for atom in mol.GetAtoms():
        cm += atom.GetDoubleProp('partial_acm')*atom.GetIntProp('Si')

    return cm

def score_mol(mol):
    """ Calculate per atom terms and score the mol """

    # Check if equivalent positions are labeled and if not do so
    if not mol.GetAtomWithIdx(0).HasProp('topo_id'):
        calc_atom_equivalence(mol)

    for atom in mol.GetAtoms():
        calc_local_terms(atom)
    
    distribute_chirality(mol)

    cm = calc_cm(mol)

    return cm

def score_mols(mols):
    """ Calculate bottcher scores across multiple mols """
    scores = []
    for mol in mols:
        scores.append(score_mol(mol))
    return scores
        