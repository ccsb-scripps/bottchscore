from rdkit import Chem
from math import log2

import logging
logger = logging.getLogger()

def calc_atom_equivalence(mol):
    """ Uses canonical ranking to describe equivalent atomic environments"""
    #TODO symmetry which arises from tautomers and protomers is not represented here
    topo_ids = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    for atom in mol.GetAtoms():
        topo_id = topo_ids[atom.GetIdx()]
        symm_count = topo_ids.count(topo_id) - 1
        atom.SetIntProp('topo_id',topo_id)
        atom.SetIntProp('num_symmetry',symm_count)
        logging.debug(f'Symmetry detection for {atom.GetSymbol()}{atom.GetIdx()}: ID={topo_id}, Number of symmetric sites={symm_count}')


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

def calc_Si(atom):
    """ Calculate the chirality term of an atom """
    #TODO potential chiral centers are not accounted for here, and I'm not sure how we could (science-wise)
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

def calc_contribution(atom):
    """ Calculate the atomic contribution to the Bottcher Score """
    di = calc_Di(atom)
    vi = calc_Vi(atom)
    si = calc_Si(atom)
    bi = calc_Bi(atom)
    ei = calc_Ei(atom)
    logging.debug(f'Atomic terms for {atom.GetSymbol()}{atom.GetIdx()}: {di=}, {ei=}, {si=}, {vi=}, {bi=}')
    if bi == 0 or vi == 0:
        logging.debug(f'Atomic contribution for {atom.GetSymbol()}{atom.GetIdx()} is 0 by definition')
        return 0
    acm =di*si*ei*log2(bi*vi)
    if atom.GetIntProp('num_symmetry') > 0:
        acm = 0.5*acm
    
    logging.debug(f'Atomic contribution for {atom.GetSymbol()}{atom.GetIdx()}: {acm}')
    return acm

def score_mol(mol):
    """ Calculate per atom terms and score the mol """
    calc_atom_equivalence(mol)

    cm = 0
    for atom in mol.GetAtoms():
        cm += calc_contribution(atom)
    
    return cm

def score_mols(mols):
    """ Calculate bottcher scores across multiple mols """
    scores = []
    for mol in mols:
        scores.append(score_mol(mol))
        return scores
        