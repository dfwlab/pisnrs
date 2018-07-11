# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 21:29:54 2018

@author: Dingfeng Wu
"""
import os
import pickle
import numpy as np
import pandas as pd
from sklearn import ensemble
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

class pisnrs:
    def __init__(self):
        self.BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace('\\', '/')+'/pisnrs/'
        self.model_path = 'Model_seq_cut1_RFC.pkl'
        self.pdes_path = 'seq.pkl'
        self.loadNRModel(self.BASE_DIR+self.model_path)
        self.loadProteinDescriptor(self.BASE_DIR+self.pdes_path)
    
    ### Load Model and descriptors
    def loadNRModel(self, path):
        ### Load model : sklearn.ensemble.RandomForestClassifier, EC50 cutoff=1, Sequence similarity descriptor
        self.sklmodel, importance, self.features = pickle.load(open(path, 'rb'), encoding='bytes')
        return self.sklmodel, self.features
    
    def loadProteinDescriptor(self, path):
        ### Load protein descriptors : Sequence similarity descriptor 
        self.seq_descriptor = pickle.load(open(path, 'rb'), encoding='bytes')
        return self.seq_descriptor
    
    def getNRs(self):
        return sorted(self.seq_descriptor.keys())
    
    def getLigandDescriptors(self):
        ligand_feature = [i for i in self.features if 'NR'!=i[:2]]
        return ligand_feature
    
    ### read ligand file and calculate Descriptor
    def calLigandDescriptor(self, mol, method):
        ### calculate Ligand Descriptors by RDKit
        return eval('Descriptors.'+method.strip())(mol)
    
    def calLigandDescriptors(self, mol):
        ligand_feature = [i for i in self.features if 'NR'!=i[:2]]
        ligand_des = {}
        for method in ligand_feature:
            try:
                value = self.calLigandDescriptor(mol, method)
            except:
                value = 0
            ligand_des[method] = value
        return ligand_des
    
    def calPCMDescriptors(self, mol, protein):
        ### calculate protein and ligand descriptors
        ligand_des = self.calLigandDescriptors(mol)
        des = []
        for f in self.features:
            if 'NR'!=f[:2]:
                des.append(ligand_des[f])
            else:
                des.append(self.seq_descriptor[protein][f])
        return des
    
    def getMolFromSmiles(self, smile):
        return Chem.MolFromSmiles(smile)

    def checkLigandType(self, text):
        ### checking input ligand file
        if '.mol' in text:
            return 'mol'
        if '.sdf' in text:
            return 'sdf'
        if 'V2000' in text or 'V3000' in text:
            return 'block'
        return 'smiles'
    
    def getMolFromText(self, text, moltype=None):
        try:
            if moltype == None:
                moltype = self.checkLigandType(text)
            if moltype == 'smiles':
                return Chem.MolFromSmiles(text)
            if moltype == 'mol':
                return Chem.MolFromMolFile(text)
            if moltype == 'block':
                return Chem.MolFromMolBlock(text) 
            if moltype == 'sdf':
                return Chem.SDMolSupplier(text)[0]
        except:
            return False
        return False
    
    def calPCMDecriptorFromMolText(self, MolText, protein, moltype=None):
        ### moltype = 'smiles', 'mol', 'block', 'sdf'
        mol = self.getMolFromText(MolText, moltype)
        if mol == False or mol == None:
            return False
        des = self.calPCMDescriptors(mol, protein)
        return des
    
    ### Predict activities
    def preLabel(self, des):
        ### predict the class of ligand activity [1 : Positive, 0 : Negtive]
        try:
            value = self.sklmodel.predict([des])
            return value[0]
        except:
            return False
    
    def preProba(self, des):
        ### predict the probability of ligand activity [Negtive probability, Postive probability]
        try:
            value = self.sklmodel.predict_proba([des])
            return (value[0][0], value[0][1])
        except:
            return False
    
    def calScaffoldFromSmile(self, smile):
        ### calculate scaffold of ligand
        return MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smile)

    def image_from_smile(self, smile, name, dir='temp/'):
        ### create molecule structure image
        path = dir+name
        try:
            mol = Chem.MolFromSmiles(smile)
            AllChem.Compute2DCoords(mol)
            Draw.MolToFile(mol, path)
        except:
            return False
        return True
    
    def preBatch(self, smiles_list, protein_list = None):
        ### batch prediction
        if protein_list == None:
            protein_list = sorted(self.seq_descriptor.keys())
        Result = []
        index = 0
        for smiles in smiles_list:
            result_row = [index, smiles]
            try:
                scaffold = self.calScaffoldFromSmile(smiles)
                result_row.append(scaffold)
            except:
                result_row.append('')
            for protein in protein_list:
                try:
                    des = self.calPCMDecriptorFromMolText(smiles, protein, moltype='smiles')
                    if des:
                        predict_proba = self.preProba(des)
                        result_row.append(round(predict_proba[1], 3))
                    else:
                        result_row.append('NaN')
                except:
                    result_row.append('NaN')
            Result.append(result_row)
            index += 1
        columns = ['index', 'smiles', 'scaffold']
        columns.extend(protein_list)
        Result = pd.DataFrame(Result, columns=columns)
        return Result


if __name__ == '__main__':
    model = pisnrs()
    print(model.BASE_DIR)
    mol = model.getMolFromSmiles('CCCC')
    print(model.calPCMDescriptors(mol=mol, protein='NR1C1'))
    print(model.calPCMDecriptorFromMolText('CCCC', protein='NR1C1'))
    des = model.calPCMDecriptorFromMolText('CCCC', protein='NR1C1')
    print(model.preProba(des))
    
    print(model.preBatch(['CCCC', 'CCC'], protein_list=['NR1C1', 'NR1C2']))
    
    print(model.getLigandDescriptors())
    print(model.getNRs())
    
    #model.image_from_smile('CCCC', 'test.png', dir='Temp/')
    
    smiles = 'CC1OC(C2=CC=CC=C2)=NC=1CN(CC1=CC(=C(C(=C1)C)OC(C(O)=O)(C)C)C)CC1OC=CC=1'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(smiles, protein, moltype='smiles') # create descriptors
    print(model.preProba(des)) # predict
    
    
    molfile = '../example/example.mol'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(molfile, protein, moltype='mol')
    print(model.preProba(des))
    
    molfile = '../example/example.sdf'
    protein = 'NR1C1'
    des = model.calPCMDecriptorFromMolText(molfile, protein, moltype='sdf')
    print(model.preProba(des))
