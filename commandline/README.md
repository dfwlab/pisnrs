# PisNRs Commandline

The commandline.py file could process a batch of compounds on the commandline. 


## Prerequisite 

Currently, commandline.py requires the following dependencies:

- Python (>= 3.6)
- PisNRs 

After the installation of pisnrs, just download [commandline.py](https://raw.githubusercontent.com/ddhmed/pisnrs/master/commandline/commandline.py) and put it into a specific directory.

## Options of commandline

### 1. -h, --help
Show this help message and exit.
#### Example
~~~~~~~~~~~~~~~
    python commandline.py --help
~~~~~~~~~~~~~~~

### 2. --mols-type 
The type of input file.Three types are allowed: molfile, smiles, sdf.  

Smiles of molecules need to be stored in a txt file.     

Mol files and sdf files need to be stored in a directory.

### 3. --mols-input
Input filename or directory about molecules.

For smiles of molecules, please input one filename and each line in the file is a complete SMILES string.

For Mol files and sdf files, please input one directory in which is filled with mol or sdf files.
### 4. --protein-input
Input filename about protein. please input one filename and each line in the file is a protein id string.
### 5. --output-filename
The content of output has several columns. The first column represents the names of molecules. 
For SMILES, the names of molecules are smiles.
For sdf or mol, the names are filenames.
Other columns represents the positive possiblities to target at the corresponding proteins.

## Example 
~~~~~~~~~~~~~~~
    python commandline.py --mols-type smiles --mols-input molecule.txt --protein-input protein.txt --output-filename output.txt
~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~
    python commandline.py --mols-type mol --mols-input ./molecule/ --protein-input protein.txt --output-filename output.txt
~~~~~~~~~~~~~~~