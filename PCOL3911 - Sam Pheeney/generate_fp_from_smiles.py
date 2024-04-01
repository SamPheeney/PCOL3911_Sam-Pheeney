# these lines import the necessary libraries
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from rdkit import Chem, DataStructs
from rdkit.Chem import PandasTools, AllChem, rdFingerprintGenerator

# these lines import the csv, and then preview the first 2 rows
data = pd.read_csv("gsar_a_1049665_sm4518.csv")
print(data.head(2))

# these lines add the molecule itself as a new column to the smiles data, and then previews the first row
PandasTools.AddMoleculeColumnToFrame(data,'SMILES','Molecule')
# print(data[["SMILES","Molecule"]].head(1))

if data.Molecule.isna().sum() != 0:
    print("An issue occurred: found empty data in the SMILES csv")
    sys.exit()

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
mfp_list = []
for mol in data['Molecule']: # for every molecule in our data...
  mfp = mfpgen.GetFingerprintAsNumPy(mol)
  mfp_list.append(mfp) # ...we want to generate a fingerprint, and add it to the list of fingerprints

data['MFP'] = mfp_list # now we add the list of fingerprints back to our csv as a new column
data_fp = data["MFP"].apply(pd.Series) # converts these fingerprints to a 'Series' datatype
data_fp.insert(1024, "Expr", data["Expr"]) # add back a dropped column
# note that it was called LogA in Slade's csv, but for our new csv the column has a different name!
print(data_fp.head(2)) # preview again to see how it's changed, should have 1025 columns now

data_fp.to_csv('organophosphate_fp.csv',index=None) # save to a new csv with no indices (row names)
