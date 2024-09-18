import sys 
import argparse
import pathlib
from rdkit import Chem
import bottchscore.bottchscore as bs


import logging
logger = logging.getLogger()

def cmd_lineparser():
    parser = argparse.ArgumentParser(description='Calculates Bottcher score for a set of molecules',
                                    epilog="""
        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2023 Forli Lab, Center for Computational Structural Biology,
                             Scripps Research.""")
    
    parser.add_argument('-i', dest='input', required=True, action='store',
                        help='input structure, in .sdf format, or smiles, in string or .smi format')
    
    parser.add_argument('-o', dest='output', action='store',
                        help='Output filename')
    
    parser.add_argument('-d', dest='debug', action='store_true')
    
    return parser.parse_args()

def main():
    args = cmd_lineparser()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    # Parse input and create supplier
    extension = pathlib.Path(args.input).suffix
    if extension == ".sdf":
        supplier = Chem.SDMolSupplier(args.input)
    elif extension == ".mol":
        supplier = [Chem.MolFromMolFile(args.input)]
    elif extension == ".smi":
        supplier = Chem.SmilesMolSupplier(args.input)
    elif extension == ".cxsmiles":
        supplier = Chem.SmilesMolSupplier(args.input, delimeter='\t')
    else:
        mol = Chem.MolFromSmiles(args.input)
        if mol is None:
            print("Input parsed as SMILES string, but conversion to RDKit mol failed.")
            print("The SMILES might be incorrect.")
            print("If you want to pass a filename, its extension must be .sdf/.mol/.smi.")
            sys.exit()
        supplier = [mol]
    
    scores = bs.score_mols(supplier)

    if args.output == None:
        print(scores)
    else:
        with open(args.output, 'w') as fout:
            for score in scores:
                fout.write(f'{score}\n')
    return

if __name__ == "__main__":
    sys.exit(main())

    
