from rdkit import Chem

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

def enumerate_unassigned_stereo(mol):
    """ Add tags to unspecified potential stereo centers to include in calculation """

    Chem.FindPotentialStereoBonds(mol) # In places adds STEREOANY

    cc = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for center in cc:
        atom = mol.GetAtomWithIdx(center[0])
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            #TODO Currently assumes TETRAHEDRAL tag will be appropriate, unimportant for this codebase but may cause issues elsewhere
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL)

def assign_atropisomerism(mol, atrop_smarts, atrop_indices):
    """ Adds an atropisomerism flag at designated bonds """
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(atrop_smarts))
    if len(matches) == 0:
        logging.warning(f'No match found for atropisomer smarts {atrop_smarts}')
        return None
    
    for m in matches:
        bond = mol.GetBondBetweenAtoms(m[atrop_indices[0]],m[atrop_indices[1]])
        if bond == None:
            logging.warning(f'No bond found for atropisomer smarts {atrop_smarts} and indices {atrop_indices}')
            continue
        bond.SetBoolProp('Atropisomeric', True)



def calulcate_mesomery(mol):
    """ Matches the molecular graph with itself in absence of bond orders to calculate equivalent atoms """

    return None

