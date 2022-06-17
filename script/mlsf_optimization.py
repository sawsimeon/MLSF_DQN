#### selecting optimal features

from msilib import OpenDatabase
import os, sys, time, fnmatch, tarfile, oddt
from glob import glob
from tqdm import tqdm
from joblib import Parallel, delayed
from oddt.pandas as opd
from oddt.pandas import ChemDataFrame
from oddt.fingerprints import PLEC, InteractionFingerprint, SimpleInteractionFingerprint, SPLIF
from oddt.toolkits import ob
from oddt.scoring.descriptors import binana


protein = '/home/simeons/Desktop/MLSF_DQN/script/data/protein/1M17_insert_addH_modelledHis_Pro.mol2'
sdf = '/home/simeons/Desktop/MLSF_DQN/script/data/compound/egfr_ligand_library.sdf'

class featurizers:
  def __init__(self, sdf, protein):
    self.sdf = sdf
    self.protein = protein
    self.mols = opd.read_sdf(sdf)['mol']
    self.ID = opd.read_sdf(sdf)['mol_name']
    '''
    Calculate Protein Ligand interaction Fingerprint
    '''
    def plec():
      receptor = next(ob.readfile('mol2', self.protein))
      receptor.protein = True
      features_list = mols.map
      features_list =  [PLEC(x, protein = receptor, size = 4092, depth_protein = 5, depth_ligand = 1, distance_cutoff = 4.5, sparse = False) for x in mols]
      features = pd.DataFrame([mol in self.mols])

      rec
    receptor = next(ob.readfile('mol2', self.protein))
    receptor.protein = True
    sdf = opd.read_sdf(se)