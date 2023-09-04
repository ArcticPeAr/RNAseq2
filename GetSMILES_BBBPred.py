import sys
import requests
from urllib.parse import urlencode
import pandas as pd
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import chemspipy
import numpy as np
from molvs import standardize_smiles
#for timing and connection is up
import time
import socket
import urllib.error
import openpyxl
#############################################################################
#Functions
#############################################################################
def predictBBB(smiles):
    '''
    Function to predict BBB penetration from smiles by LightBBB
    '''
    serverURL = 'http://165.194.18.43:7020/bbb_permeability?'
    param = urlencode({'smiles': smiles})
    response = requests.get(serverURL + param)
    result = response.text.strip("\n").strip('\"')
    return result

def getGCTFile(filePath):
    '''
    Function to open GCT file from Clue.io and get values from it
    '''
    #Import gct file with p-values
    importedGCT = pd.read_csv(filePath, sep='\t', skiprows=2)
    #keep only columns 1 and 2
    GCT = importedGCT[['pert_id', 'pert_iname']]
    #remove 1st row
    GCT = GCT.drop(GCT.index[0])
    #remove duplicate rows
    GCT = GCT.drop_duplicates(subset=['pert_id'])
    GCT = GCT.drop_duplicates(subset=['pert_iname'])
    GCT2 = pd.DataFrame()
    GCT2['Pert'] = [val for pair in zip(GCT['pert_id'], GCT['pert_iname']) for val in pair]
    GCT2 = GCT2.drop_duplicates(subset=['Pert'])
    return GCT2


def check_network_connection():
    '''
    Function to check if network is connected
    '''
    while True:
        try:
            # Attempt to connect to a reliable host
            socket.create_connection(("google.com", 80))
            print("Network is connected.")
            break
        except OSError:
            # Network is unreachable
            print("Network is unreachable. Retrying in 5 seconds...")
            time.sleep(5)

import requests
from rdkit import Chem

def name_to_smiles(compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/IsomericSMILES/txt"
    response = requests.get(url)
    
    if response.status_code == 200:
        smiles = response.text.strip()
        return smiles
    else:
        return None

# Example usage
compound_name = "aspirin"
smiles = name_to_smiles(compound_name)

if smiles:
    # Load the molecule in RDKit for further processing
    molecule = Chem.MolFromSmiles(smiles)
    print(f"The SMILES string for {compound_name} is {smiles}")
else:
    print(f"Could not find SMILES string for {compound_name}")


def pubchemSmiles(pert):
    '''
    Function to get smiles from pubchempy and add to dataframe in new column
    '''
    try:
        compound = pcp.get_compounds(pert, 'name')
        smiles = compound[0].canonical_smiles
        return smiles
    except:
        return np.nan



def add_smiles_to_df(df):
    # Create a new column for storing the SMILES strings, initialized with "NA"
    df['SMILES'] = "NA"
    for index, row in df.iterrows():
        compound_name = row[0]  # Assuming names are in the first column
        smiles = name_to_smiles(compound_name)
        if smiles:
            df.at[index, 'SMILES'] = smiles
        else:
            df.at[index, 'SMILES'] = "NA"
    return df



# Define the name_to_smiles function (or you can import it if defined elsewhere)
def name_to_smiles(compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/IsomericSMILES/txt"
    response = requests.get(url)
    if response.status_code == 200:
        smiles = response.text.strip()
        return smiles
    else:
        return None
    
    
#############################################################################
#Load SMILES files
#############################################################################
KnownSmilesDF1 = pd.read_csv('/home/petear/MEGA/TormodGroup/InputData/compoundinfo_beta.txt', sep='\t')
#If SMILES for substances not found in KnownSmilesDF1, then try KnownSmilesDF2 and KnownSmilesDF3
KnownSmilesDF2 = pd.read_excel('/home/petear/MEGA/TormodGroup/InputData/SmilesPerts2.xlsx')
KnownSmilesDF3 = pd.read_csv('/home/petear/MEGA/TormodGroup/InputData/SmilesPerts.txt', sep='\t')
  

#############################################################################
#Load data
#############################################################################
#Clue got FDR in the 9th column
clue = pd.read_csv('/home/petear/MEGA/TormodGroup/InputData/AllAlgoSi09Aug.csv', header=None)


#############################################################################
#Find MoA perturbagens and merge with clue dataframe 
#############################################################################
wanted_MoA = "MTOR inhibitor"
unwanted_MoA = "PI3K inhibitor"

#File with perturbagens per MoA
pertMoA = pd.read_excel('/home/petear/MEGA/TormodGroup/InputData/MoACounts.xlsx', sheet_name='PertsPrMoAAll')

#Remove all columns that dont have MoA in them
cols_to_drop1 = [col for col in pertMoA if MoA not in col]
cols_to_drop2 = [col for col in pertMoA if unwanted_MoA in col]

pertMoA.drop(columns=cols_to_drop, inplace=True)
pertMoA.drop(columns=cols_to_drop2, inplace=True)

#take all values from pertMoA dataframe and put them in a list
pertMoA = pertMoA.values.tolist()

#flatten list
pertMoA = [item for sublist in pertMoA for item in sublist]

#remove duplicates
pertMoA = list(dict.fromkeys(pertMoA))

#remove NaN
pertMoA = [x for x in pertMoA if str(x) != 'nan']

#remove all rows from clue dataframe that are not in pertMoA list
clue = clue[clue[0].isin(pertMoA)]

#remove duplicates from clue dataframe and keep the highest value in the last column
clue = clue.drop_duplicates(subset=[0], keep='first')

#############################################################################
#First iteration of Smiles from pubchempy
#############################################################################

clue = add_smiles_to_df(clue)


#############################################################################
#Second iteration of Smiles from other files
#############################################################################
#add new column to clue dataframe and fill with NaN
clue['SMILESFromFiles1'] = np.nan
clue['SMILESFromFiles2'] = np.nan
clue['SMILESFromFiles3'] = np.nan

#Use KnownSmilesDF1["cmap_name"] to fill in SMILESFromFiles column for if these values are in clue[0] with for loop
for i in range(len(KnownSmilesDF1["cmap_name"])):
    if KnownSmilesDF1["cmap_name"][i] in clue[0].values:
        clue.loc[clue[0] == KnownSmilesDF1["cmap_name"][i], 'SMILESFromFiles1'] = KnownSmilesDF1["canonical_smiles"][i]
    else:
        pass
    
#For items in KnownSmilesDF2["pert_id"] that are in clue[0], add the corresponding KnownSmilesDF2["CanonicalSMILES"] to SMILESFromFiles2 column with for loop
for i in range(len(KnownSmilesDF2["pert_id"])):
    if KnownSmilesDF2["pert_id"][i] in clue[0].values:
        clue.loc[clue[0] == KnownSmilesDF2["pert_id"][i], 'SMILESFromFiles2'] = KnownSmilesDF2["CanonicalSMILES"][i]
    else:
        pass

#For items in KnownSmilesDF3["SM_Center_Canonical_ID"] that are in clue[0], add the corresponding KnownSmilesDF3["SM_Center_Canonical_SMILES"] to SMILESFromFiles3 column with for loop
for i in range(len(KnownSmilesDF3["SM_Center_Canonical_ID"])):
    if KnownSmilesDF3["SM_Center_Canonical_ID"][i] in clue[0].values:
        clue.loc[clue[0] == KnownSmilesDF3["SM_Center_Canonical_ID"][i], 'SMILESFromFiles3'] = KnownSmilesDF3["SM_Center_Canonical_SMILES"][i]
    else:
        pass



#############################################################################
#BBBPred
#############################################################################
#remove columns if all values are NaN
clue = clue.dropna(axis=1, how='all')

#For each smile in SMILESFromFiles1, predict BBB penetration and add to new column
clue['BBBPredFromFile'] = clue['SMILES'].apply(predictBBB)

clue['BBBPredSmileFromPubchem'] = clue['SMILES'].apply(predictBBB)


#############################################################################
#Save to excel
#############################################################################

clue.to_excel('/home/petear/MEGA/TormodGroup/InputData/BBBPred.xlsx', index=False)