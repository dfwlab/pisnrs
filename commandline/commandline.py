import pisnrs
import argparse
import os
import pandas as pd
from rdkit import Chem

def collect_and_check_protein(protein_filename, pisnrs_model):
    result_list = []
    NR_protein_list = sorted(pisnrs_model.seq_descriptor.keys())
    if not (os.path.isfile(protein_filename)):
        print("Invalid protein filename, use all nuclear receptors")
        result_list = NR_protein_list
        return result_list
    f = open(protein_filename, 'r')
    NR_protein_list = sorted(pisnrs_model.seq_descriptor.keys())
    for line in f.readlines():
        if (line.strip() not in NR_protein_list):
            print("error when loading proteins \"{}\". skip this protein".format(line.strip()))
            continue
        result_list.append(line.strip())
    if not result_list:
        print("No proteins input in the file, use all nuclear receptors")
        result_list = NR_protein_list
    return result_list
    
def predict_by_smiles(input_mol_file, protein_list, output_filename, pisnrs_model):
    if not (os.path.isfile(input_mol_file)):
        raise Exception("Invalid moleclue filename")
    f = open(input_mol_file, 'r')
    result = {}
    result['smiles'] = []
    mol_index = 0
    for protein in protein_list:
        result[protein] = []
    for line in f.readlines():
        try:
            temp = Chem.MolFromSmiles(line.strip())
        except:
            print("error when loading smiles \"{}\". skip this smiles".format(line))
            continue
        result['smiles'].append(line.strip())
        for protein in protein_list:
            des = pisnrs_model.calPCMDecriptorFromMolText(line.strip(), protein, moltype='smiles')
            result[protein].append(round(pisnrs_model.preProba(des)[1], 3)) 
        mol_index += 1
    dataframe = pd.DataFrame(result, index = [i for i in range(mol_index)])
    dataframe.to_csv(output_filename, sep = '\t')

def predict_by_files(input_dir, protein_list, output_filename, pisnrs_model, file_type):
    if not (os.path.isdir(input_dir)):
        raise Exception("Invalid directory")
    file_list = []
    for filename in os.listdir(input_dir):
        if filename.endswith("." + file_type):
            file_list.append(filename)
    result = {}
    result[file_type] = []
    for protein in protein_list:
        result[protein] = []
    mol_index = 0
    for filename in file_list:
        try:
            temp = Chem.MolFromMolFile(os.path.join(input_dir, filename))
        except:
            print("error when loading file \"{}\". skip this file".format(os.path.join(input_dir, filename)))
            continue
        result[file_type].append(filename)
        for protein in protein_list:
            des = pisnrs_model.calPCMDecriptorFromMolText(os.path.join(input_dir, filename), protein, moltype=file_type)
            result[protein].append(round(pisnrs_model.preProba(des)[1], 3)) 
        mol_index += 1
    dataframe = pd.DataFrame(result, index = [i for i in range(mol_index)])
    dataframe.to_csv(output_filename, sep = '\t')
    
def main(args):
    protein_list = []
    pisnrs_model = pisnrs.pisnrs()
    protein_list = collect_and_check_protein(args.protein_input, pisnrs_model)
    print(protein_list)
    if args.mols_type == "none":
        print("the type of molecule file is not clear")
    elif args.mols_type == "smiles":
        predict_by_smiles(args.mols_input, protein_list, args.output_filename, pisnrs_model)
    else:
        predict_by_files(args.mols_input, protein_list, args.output_filename, pisnrs_model, args.mols_type)

parser = argparse.ArgumentParser(description="""
example: python pipeline.py --mols-type smiles --mols-input molecule.txt --protein-input protein.txt --output-filename output.txt\n
\t python pipeline.py --mols-type mol --mols-input ./molecule/ --protein-input protein.txt --output-filename output.txt\n

""", formatter_class=argparse.RawTextHelpFormatter)
model = pisnrs.pisnrs()
parser.add_argument("--mols-type", type=str, default="none", choices=["mol", "smiles", "sdf"],help="""
The type of input file.
Three types are allowed: molfile, smiles, sdf.
Smiles of molecules are stored in a txt file.
Mol files and sdf files are stored in a directory.

""")
parser.add_argument("--mols-input", type=str,default=None,help="""
Input filename or directory about molecules. 
For smiles of molecules, please input one filename and each line in the file is a complete SMILES string.
For Mol files and sdf files, please input one directory in which is filled with mol or sdf files.

""")
parser.add_argument("--protein-input", type=str, default="none",help="""
Input filename about protein. 
please input one filename and each line in the file is a protein id string.

""")
parser.add_argument("--output-filename", default="output.txt", type=str,help="""
The content of output has several columns. 
The first column represents the names of molecules. 
For SMILES, the names of molecules are smiles.
For sdf or mol, the names are filenames.
Other columns represents the positive possiblities to target at the corresponding proteins.
""")

args = parser.parse_args()
if __name__ == "__main__":
    main(args)







