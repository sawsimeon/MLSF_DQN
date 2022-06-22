#### selecting optimal features
import os, sys, time, fnmatch, tarfile, oddt
from pyexpat import features
from glob import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import oddt.pandas as opd
from oddt.pandas import ChemDataFrame
from oddt.fingerprints import PLEC, InteractionFingerprint, SimpleInteractionFingerprint, SPLIF
from oddt.toolkits import ob
from oddt.scoring.descriptors import binana, close_contacts_descriptor, oddt_vina_descriptor
import pandas as pd
from oddt.scoring import ensemble_descriptor
import numpy as np


protein = '/home/simeons/Desktop/MLSF_DQN/script/data/protein/1M17_insert_addH_modelledHis_Pro.mol2'
sdf = '/home/simeons/Desktop/MLSF_DQN/script/data/compound/egfr_ligand_library.sdf'

class featurizers:
  def __init__(self, sdf, protein):
    self.sdf = opd.read_sdf(sdf)
    self.protein = next(ob.readfile('mol2', protein))
    self.mols = opd.read_sdf(sdf)['mol']
    self.ID = opd.read_sdf(sdf)['mol_name']
  
  def PLEC(self):
    features =  pd.DataFrame([PLEC(x, protein = self.protein, size = 4092, depth_protein = 5, depth_ligand = 1, distance_cutoff = 4.5, sparse = False) for x in self.mols])
    features.insert(0, "ChEMBLID", list(self.ID))
    return features

  '''
  Calculate Interaction Fingerprint
  '''
  def IF(self):
    features = pd.DataFrame([InteractionFingerprint(x, protein = self.protein) for x in self.mols])
    features.insert(0, 'ChEMBLID', list(self.ID))
    return features
    
  '''
  Calculate Simple interaction Fingerprint 
  http://dx.doi.org/10.1016/j.csbj.2014.05.004
  '''
  def SFP(self):
    features = pd.DataFrame([SimpleInteractionFingerprint(x, protein = self.protein) for x in self.mols])
    features.insert(0, 'ChEMBLID', list(self.ID))
    return features
  '''
  Calculate Binina Features 
  '''
  def BINANA(self):
    binana_engine = binana.binana_descriptor(self.protein)
    features = pd.DataFrame([{name: value for name, value in zip(binana_engine.titles, binana_engine.build([x])[0])}  for x in self.mols])
    features.insert(0, 'ChEMBLID', list(self.ID))
    return features
  '''
  Calculate RF Score Features
  '''
  def RFSCORE(self, version):
    ligand_atomic_nums = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    protein_atomic_nums = [6, 7, 8, 16]
    cutoff = 12

    if version == 1:
      cutoff = 12
      rfscore_engine = close_contacts_descriptor(protein = self.protein, cutoff = cutoff, protein_types = protein_atomic_nums, ligand_types = ligand_atomic_nums)
      features = pd.DataFrame([{name: value for name, value in zip(rfscore_engine.titles, rfscore_engine.build([x])[0])}  for x in self.mols])
      features.insert(0, 'ChEMBLID', list(self.ID))
      return features
    elif version == 2:
      cutoff = np.array([0, 2, 4, 6, 8, 10, 12])
      rfscore_engine = close_contacts_descriptor(protein = self.protein, cutoff = cutoff, protein_types = protein_atomic_nums, ligand_types = ligand_atomic_nums)
      features = pd.DataFrame([{name: value for name, value in zip(rfscore_engine.titles, rfscore_engine.build([x])[0])}  for x in self.mols])
      features.insert(0, 'ChEMBLID', list(self.ID))
      return features
    elif version == 3:
      cutoff = 12
      cc = close_contacts_descriptor(protein = self.protein, cutoff = cutoff, protein_types = protein_atomic_nums, ligand_types = ligand_atomic_nums)
      vina_scores = ['vina_gauss1',
                           'vina_gauss2',
                           'vina_repulsion',
                           'vina_hydrophobic',
                           'vina_hydrogen',
                           'vina_num_rotors']
      vina = oddt_vina_descriptor(protein = self.protein, vina_scores = vina_scores)
      rfscore_engine = ensemble_descriptor((vina, cc))
      features = pd.DataFrame([{name: value for name, value in zip(rfscore_engine.titles, rfscore_engine.build([x])[0])} for x in self.mols])
      features.insert(0, 'ChEMBLID', list(self.ID))
      return features


        