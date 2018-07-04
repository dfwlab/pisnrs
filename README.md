# pisnrs

pisnrs is a Python module for NR inhibitors prediction built on top of
RDKit, scikit-learn and distributed under the MIT license.


## Installation

### Dependencies

pisnrs requires:

- Python (>= 2.6)
- rdkit (>= 2018.03.2.0)
- scikit-learn (>= 0.19.1)
- pandas (>= 0.23.1)

You can configure the environment of the pisnrs through Ancaconda.

### 1. Introduction to anaconda

Conda is an open-source, cross-platform, software package manager. It supports the packaging and distribution of software components, and manages their installation inside isolated execution environments. It has several analogies with pip and virtualenv, but it is designed to be more "python-agnostic" and more suitable for the distribution of binary packages and their dependencies.

### 2. How to get conda

The easiest way to get Conda is having it installed as part of the [Anaconda Python distribution](https://conda.io/docs/user-guide/install/index.html). A possible (but a bit more complex to use) alternative is provided with the smaller and more self-contained [Miniconda](https://conda.io/miniconda.html). The conda source code repository is available on [github](https://github.com/conda) and additional documentation is provided by the project [website](https://conda.io/docs/).

### 3. How to install RDKit with Conda

Creating a new conda environment with the RDKit installed requires one single command similar to the following:

~~~~~~~~~~~~~~~
  $ conda create -c rdkit rdkit
~~~~~~~~~~~~~~~

### 4. Install scikit-learn

Install scikit-learn using ``pip`` :

~~~~~~~~~~~~~~~
    pip install -U scikit-learn
~~~~~~~~~~~~~~~

or ``conda`` :

~~~~~~~~~~~~~~~
    conda install scikit-learn
~~~~~~~~~~~~~~~

### 5. Install pisnrs

If you already have a working installation of rdkit and scikit-learn, the easiest way to install pisnrs is using ``pip`` :

~~~~~~~~~~~~~~~
    pip install pisnrs
~~~~~~~~~~~~~~~

## Example

### 1. import pisnrs and load model

~~~~~~~~~~~~~~~
    import pisnrs
    model = pisnrs.pisnrs()
    print(model.getNRs()) #print NR proteins in model
    print(model.getLigandDescriptors()) # print Ligand descriptors in model
~~~~~~~~~~~~~~~

### 2. predict the activity and scaffold of one ligand smiles

~~~~~~~~~~~~~~~
    smiles = 'CCCCC'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(smiles, protein) # create descriptors
    print(model.preProba(des)) # predict
~~~~~~~~~~~~~~~

### 3. Batch mode

~~~~~~~~~~~~~~~
    smiles_list = ['CCCC', 'CCC']
    protein_list = ['NR1C1', 'NR1C2']
    print(model.preBatch(smiles_list, protein_list=protein_list) # predict activity of every ligands and proteins in list 
    print(model.preBatch(smiles_list) # predict activity of every ligands in list and all proteins in model
~~~~~~~~~~~~~~~

### 4. Create molecule images

~~~~~~~~~~~~~~~
    smile = 'CCCC'
    model.image_from_smile(smile, name='my_ligand.png', dir='temp/') # output image of smiles to 'temp/' folder
~~~~~~~~~~~~~~~


## Important links

- Official source code repo: https://github.com/ddhmed/pisnrs
- Download releases: https://pypi.org/project/pisnrs/

## Source code

You can check the latest sources with the command:

    git clone https://github.com/ddhmed/pisnrs.git

