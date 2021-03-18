# PisNRs

PisNRs is a Python module for potential inhibitor and scaffold prediction 
of nuclear receptors, which constructed on top of RDKit, scikit-learn 
under the MIT license.


## Installation

Currently, PisNRs requires the following dependencies:

- Python (>= 3.6)
- rdkit (>= 2018.03.2.0)
- scikit-learn (>= 0.19.1)
- pandas (>= 0.23.1)

Ancaconda is recommended for package management and environment configures.

### 1. Anaconda introduction

The user can download and install Anaconda at [Anaconda Python distribution](https://conda.io/docs/user-guide/install/index.html).
Also, Miniconda is acceptable in [here](https://conda.io/miniconda.html). The conda source code repository is avaliable at
 [github](https://github.com/conda) and project [website](https://conda.io/docs/).

### 2. Install RDKit with Anaconda

Creating a new conda environment with the RDKit installation with the following command:

~~~~~~~~~~~~~~~
  $ conda create -c rdkit rdkit
~~~~~~~~~~~~~~~

### 3. Install scikit-learn

User can install scikit-learn by using ``pip`` :

~~~~~~~~~~~~~~~
    pip install -U scikit-learn
~~~~~~~~~~~~~~~

or ``conda`` :

~~~~~~~~~~~~~~~
    conda install scikit-learn
~~~~~~~~~~~~~~~

### 4. Install PisNRs

After installation of RDKit and scikit-learn, PisNRs can be installed by using ``pip`` :

~~~~~~~~~~~~~~~
    pip install --upgrade pisnrs
~~~~~~~~~~~~~~~

## Example

### 1. import PisNRs and load model

~~~~~~~~~~~~~~~
    import pisnrs
    model = pisnrs.pisnrs()
    print(model.getNRs()) #print NR proteins in model
    print(model.getLigandDescriptors()) # print Ligand descriptors in model
~~~~~~~~~~~~~~~

### 2. Predict the activity and scaffold of query ligand

~~~~~~~~~~~~~~~
    ### moltype includes : smiles, mol, block, sdf
    ### You can find the example .mol and .sdf file in the 'example/' folder
    ### Example: https://github.com/ddhmed/pisnrs/tree/master/example
    # 1. smiles input
    smiles = 'CC1OC(C2=CC=CC=C2)=NC=1CN(CC1=CC(=C(C(=C1)C)OC(C(O)=O)(C)C)C)CC1OC=CC=1'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(smiles, protein, moltype='smiles') # create descriptors
    print(model.preProba(des)) # predict
    
    # 2. .mol file input
    molfile = 'example/2.2_MolFile.mol'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(molfile, protein, moltype='mol')
    print(model.preProba(des))
    
    # 3. mol block input
    block = '''
     RDKit          2D

 36 39  0  0  0  0  0  0  0  0999 V2000
    4.9515   -5.4554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4879   -6.8820    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3696   -8.0956    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4879   -9.3091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9515  -10.7357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4187  -11.0475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8822  -12.4741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8785  -13.5888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4113  -13.2770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9478  -11.8504    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0613   -8.8456    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.0613   -7.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8478   -6.4639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4775   -7.0740    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3207   -8.5658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0496   -9.1759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2631   -8.2942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6335   -8.9043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7902  -10.3961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5767  -11.2778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2064  -10.6676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7335  -12.7695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1606  -11.0062    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3174  -12.4980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4741  -13.9897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8445  -14.5999    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2606  -14.8714    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8091  -12.3412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8256  -12.6548    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8470   -8.0226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7360   -6.1923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5792   -4.7005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6939   -3.6968    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0838   -2.3265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4079   -2.4833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7198   -3.9505    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  4 11  2  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 16 17  2  0
 17 18  1  0
 18 19  2  0
 19 20  1  0
 20 21  2  0
 20 22  1  0
 19 23  1  0
 23 24  1  0
 24 25  1  0
 25 26  1  0
 25 27  2  0
 24 28  1  0
 24 29  1  0
 18 30  1  0
 14 31  1  0
 31 32  1  0
 32 33  1  0
 33 34  1  0
 34 35  2  0
 35 36  1  0
 12  2  2  0
 21 16  1  0
 36 32  2  0
 10  5  1  0
M  END
    '''
    des = model.calPCMDecriptorFromMolText(block, protein, moltype='block')
    print(model.preProba(des))
    
    # 4. .sdf file input
    sdffile = 'example/2.4_SDFFile.sdf'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(sdffile, protein, moltype='sdf')
    print(model.preProba(des))
~~~~~~~~~~~~~~~

### 3. Derive molecule scaffold of query ligand

~~~~~~~~~~~~~~~
    smiles = 'CC1OC(C2=CC=CC=C2)=NC=1CN(CC1=CC(=C(C(=C1)C)OC(C(O)=O)(C)C)C)CC1OC=CC=1'
    scaffold = model.calScaffoldFromSmiles(smiles)
    print(scaffold)
~~~~~~~~~~~~~~~

### 4. Create molucule image of query ligand

~~~~~~~~~~~~~~~
    smiles = 'CC1OC(C2=CC=CC=C2)=NC=1CN(CC1=CC(=C(C(=C1)C)OC(C(O)=O)(C)C)C)CC1OC=CC=1'
    model.image_from_smiles(smiles, name='4_OutputImage.png', dir='example/') # output image of smiles to 'example/' folder
~~~~~~~~~~~~~~~

### 5. Batch mode uploading

~~~~~~~~~~~~~~~
    smiles_list = ['CC1OC(C2=CC=CC=C2)=NC=1CN(CC1=CC(=C(C(=C1)C)OC(C(O)=O)(C)C)C)CC1OC=CC=1', 'C1=CC=CC=C1']
    protein_list = ['NR1C1', 'NR1C2']
    print(model.preBatch(smiles_list, protein_list)) # predict activity of every ligands and proteins in list 
    print(model.preBatch(smiles_list)) # predict activity of every ligands in list and all proteins in model
    ### load smiles list from file
    smiles_file = 'example/5_SmilesList.smiles'
    smiles_list = [i.strip() for i in open(smiles_file, 'r').readlines()]
    print(model.preBatch(smiles_list))
~~~~~~~~~~~~~~~

## Related links

- Official source code repo: https://github.com/ddhmed/pisnrs
- Download releases: https://pypi.org/project/pisnrs/
- web server : http://pisnrs.badd-cao.net/predict

## Source code

The latest sources can be checked by using the following command:

    git clone https://github.com/ddhmed/pisnrs.git
