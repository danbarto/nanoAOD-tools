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

        # New collection of Muons
        self.out.branch("Lepton_pt", "F", lenVar="nLepton")
        self.out.branch("Lepton_eta", "F", lenVar="nLepton")
        self.out.branch("Lepton_phi", "F", lenVar="nLepton")
        self.out.branch("Lepton_mass", "F", lenVar="nLepton")
        self.out.branch("Lepton_pdgId", "I", lenVar="nLepton")
        self.out.branch("Lepton_miniIso", "F", lenVar="nLepton")
        self.out.branch("Lepton_muIndex", "I", lenVar="nLepton")
        self.out.branch("Lepton_elIndex", "I", lenVar="nLepton")

        # Counters
        self.out.branch("nLepton",      "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodJet(self, jet):
        return (jet.pt > 25 and abs(jet.eta)<2.4 and jet.jetId>1)

    def isFwdJet(self, jet):
        return (((jet.pt > 25 and abs(jet.eta)<2.7) or (jet.pt>60 and abs(jet.eta)>=2.7 and abs(jet.eta)<3.0) or (jet.pt>25 and abs(jet.eta)>=3.0 and abs(jet.eta)<5.0))  and jet.jetId>1)

    def isGoodBJet(self, jet):
        if self.year == 2018:
            threshold = 0.4184
        return (self.isGoodJet(jet) and jet.btagDeepB > threshold)

    def isVetoMuon(self, muon):
        return (muon.looseId and muon.pt>5 and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.2 and abs(muon.dxy)<0.1 and abs(muon.dz)<0.5)

    def isVetoElectron(self, electron):
        return (electron.cutBased>0 and electron.miniPFRelIso_all < 0.2)

    def isTightMuon(self, muon):
        return (muon.pt > 25 and muon.mediumId and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.1)

    def isTightElectron(self, electron):
        return (electron.pt > 30 and electron.cutBased >= 3 and abs(electron.eta) < 2.4 and electron.miniPFRelIso_all < 0.1)# and electron.sip3d < 4.0 and abs(electron.dxy) < 0.05 and abs(electron.dz) < 0.1)

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

    def getW(self, jet1, jet2):
        mass, pt, eta, phi = self.invMass(jet1,jet2)
        qgl_sum = jet1['qgl'] + jet2['qgl']
        qgl_prod = jet1['qgl'] * jet2['qgl']
        return {'chi2': abs(mass-80.)/20., 'mass':mass, 'pt':pt, 'eta':eta, 'phi':phi, 'qgl_sum':qgl_sum, 'qgl_prod': qgl_prod}

    def getWcandidates(self, jets):
        uniqueCombs = [(0,1,2,3), (1,2,3,0), (0,2,1,3)]
        minChi2 = 9999
        jetIndices = range(len(jets))
        sets = [ comb for comb in itertools.combinations(jetIndices, 4) ]
        W_cand = []
        #print "Number of jets:", len(jets)
        #print "All combinations of the 4 jets", sets
        for s in sets:
            # selection of 4 jets
            for combs in uniqueCombs:
                # combination of the 4 jets
                indices = [ s[x] for x in combs ]
                #print "One of the combinations", indices
                W = [self.getW(jets[indices[0]], jets[indices[1]]), self.getW(jets[indices[2]], jets[indices[3]])]
                chi2 = W[0]['chi2'] + W[1]['chi2']
                W_cand.append({'chi2': chi2, 'indices':indices, 'W':W})

        W_cand = sorted(W_cand, key=lambda x: x['chi2'])

        return W_cand
            
    def getRealWs(self, jets, genWs):
        combs = [ comb for comb in itertools.combinations(jets, 2) ]
        Wcands = []
        for j1, j2 in combs:
            Wcand = self.getW(j1, j2)
            if j1['WIdx']==j2['WIdx'] and j1['WIdx']>-1:
                Wcand['genPt'] = genWs[j1['WIdx']].pt
                Wcand['genEta'] = genWs[j1['WIdx']].eta
                Wcand['genPhi'] = genWs[j1['WIdx']].phi
                Wcand['genMass'] = genWs[j1['WIdx']].mass
            else:
                Wcand['genPt']      = -9999
                Wcand['genEta']     = -9999
                Wcand['genPhi']     = -9999
                Wcand['genMass']    = -9999
            Wcands.append(Wcand)

        return Wcands
                


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        jets        = Collection(event, "Jet")
        
        # tight lepton collection, will be sorted by pt
        leptons     = []
        for i,mu in enumerate(muons):
            leptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1, 'mass':mu.mass})

        for i,el in enumerate(electrons):
            leptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i, 'mass':el.mass})

        leptons = sorted(leptons, key = lambda i: i['pt'], reverse=True)
        
        # make pandas dataframe out of list
        leptons_pd = pd.DataFrame(leptons)

        self.out.fillBranch("nLepton",          len(leptons_pd) )
        if len(leptons_pd)>0:
            self.out.fillBranch("Lepton_pt",        leptons_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("Lepton_eta",       leptons_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("Lepton_phi",       leptons_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("Lepton_mass",      leptons_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("Lepton_pdgId",     leptons_pd.sort_values(by='pt', ascending=False)['pdgId'].tolist() )
            self.out.fillBranch("Lepton_miniIso",   leptons_pd.sort_values(by='pt', ascending=False)['miniIso'].tolist() )
            self.out.fillBranch("Lepton_muIndex",   leptons_pd.sort_values(by='pt', ascending=False)['muIndex'].tolist() )
            self.out.fillBranch("Lepton_elIndex",   leptons_pd.sort_values(by='pt', ascending=False)['elIndex'].tolist() )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selector = lambda year, isData : PhysicsObjects( year=year, isData=isData )
