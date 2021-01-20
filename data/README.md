# PisNRs dataset

This page describe the information of database used in the model construction and evaluation.


## ONRLDB
### Requirement
- Python
- pickle
### Description

There are two files storing the data of ONRLDB dataset.
The first file, [ONRLDB-chemical-smiles.pkl](https://github.com/ddhmed/pisnrs/blob/master/data/ONRLDB-chemical-smiles.pkl), needs to be opened by using the pickle package and contains the ID, the name and the smiles of chemicals.
The second file, [ONRLDB-ActivityData.txt](https://github.com/ddhmed/pisnrs/blob/master/data/ONRLDB-ActivityData.txt), contains the information about EC50 values targted at specific nuclear receptors. The first column represents the ID in the ONRLDB database. The second and third represent the name of the target protein and the chain. The fourth column represents the CID of the compound. The final column represents the EC50 value. 

## TCM database@Taiwan
The file, [TCM-Database@Taiwan.csv](https://github.com/ddhmed/pisnrs/blob/master/data/TCM-Database%40Taiwan.csv) contain the data of TCM database@Taiwan. 

## Reference
### ONRLDB
[ONRLDB.nbib](https://github.com/ddhmed/pisnrs/blob/master/data/ONRLDB.nbib)
### TCM database@Taiwan
[TCM-Database@Taiwan.nbib](https://github.com/ddhmed/pisnrs/blob/master/data/TCM-Database%40Taiwan.nbib)

