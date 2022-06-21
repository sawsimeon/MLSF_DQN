#### selecting optimal features

from msilib import OpenDatabase
import os, sys, time, fnmatch, tarfile, oddt
from pyexpat import features
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
    self.sdf = opd.read_sdf(sdf)
    self.protein = next(ob.readfile('mol2', self.protein)
    self.mols = opd.read_sdf(sdf)['mol']
    self.ID = opd.read_sdf(sdf)['mol_name']
    '''
    Calculate Protein Ligand interaction Fingerprint
    '''
    def PLEC(self):
      features =  pd.DataFrame([PLEC(x, protein = self.protein.protein = True, size = 4092, depth_protein = 5, depth_ligand = 1, distance_cutoff = 4.5, sparse = False) for x in mols])
      features.insert(1, "ChEMBLID", list(self.ID))
      return features

    '''
    Calculate Interaction Fingerprint
    '''
    def IF(self):
      features = pd.DataFrame([InteractionFingerprint(x, protein = self.protein.protein = True)])
      features.insert(1, 'ChEMBL', self.ID)
      return features
    
    '''
    Calculate Simple interaction Fingerprint 
    http://dx.doi.org/10.1016/j.csbj.2014.05.004
    '''
    def SFP(self):
      feautres = pd.DataFrame([SimpleinteractionFingerprint(x, protein = self.protein.protein = True)])
      features.insert(1,, 'ChEMBLID', list(self.ID))
      return features
    '''
    Calculate Binina Features 
    '''
    def BININA(self):
      binina_engine = binana.binana_descriptor(self.protein.protein = True)
      features = pd.DataFrame({name: value for name, value in zip(binana_engine.titles, binana_engine.build([x])[0])}, )
