import ROOT
import os
import numpy as np
import pandas as pd
import math
import glob
import itertools
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class Object():
    def __init__(self):
        pass

    @classmethod
    def fromDict(cls, d):
        obj = cls()
        for key in d.keys():
            setattr(obj, key, d[key])
        return obj

class PhysicsObjects(Module):

    def __init__(self, year=2018, isData=False):
        self.year = year
        self.isData = isData
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        ## Define a first minimum set of objects needed for the analysis
        #FIXME objects should have cross-reference to full collection


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def deltaPhi(self, phi1, phi2):
        dphi = phi2-phi1
        if  dphi > math.pi:
            dphi -= 2.0*math.pi
        if dphi <= -math.pi:
            dphi += 2.0*math.pi
        return abs(dphi)

    def deltaR2(self, l1, l2):
        return self.deltaPhi(l1.phi, l2.phi)**2 + (l1.eta - l2.eta)**2

    def deltaR(self, l1, l2):
        return math.sqrt(self.deltaR2(l1,l2))

    def invMass(self, o1, o2):
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(o1['pt'], o1['eta'], o1['phi'], o1['mass'])
        v2.SetPtEtaPhiM(o2['pt'], o2['eta'], o2['phi'], o2['mass'])
        return (v1+v2).M(), (v1+v2).Pt(), (v1+v2).Eta(), (v1+v2).Phi()

    def MCT2(self, b1, b2):
      return 2*b1['pt']*b2['pt']*(1+math.cos(self.deltaPhi(b1['phi'],b2['phi'])))

    def MCT(self, b1, b2):
      return math.sqrt(self.MCT2(b1,b2))

    def MT2(self, b1, metpt, metphi):
      return 2*b1['pt']*metpt*(1+math.cos(self.deltaPhi(b1['phi'],metphi)))

    def MT(self, b1, metpt, metphi):
      return math.sqrt(self.MT2(b1,metpt, metphi))

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        taus        = Collection(event, "Tau")
        jets        = Collection(event, "Jet")
        isotracks   = Collection(event, "IsoTrack")
        fatjets     = Collection(event, "FatJet")
        if not self.isData:
            genjets     = Collection(event, "GenJet")
        #    genW        = Collection(event, "W")
        
        # MET
        met_pt  = event.MET_pt
        met_phi = event.MET_phi


        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selector = lambda year, isData : PhysicsObjects( year=year, isData=isData )
