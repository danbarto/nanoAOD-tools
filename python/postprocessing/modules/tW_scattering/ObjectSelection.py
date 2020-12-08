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

        ##
        self.out.branch("Muon_isVeto",      "I", lenVar="nMuon")
        self.out.branch("Muon_isTight",     "I", lenVar="nMuon")
        self.out.branch("Electron_isVeto",  "I", lenVar="nElectron")
        self.out.branch("Electron_isTight", "I", lenVar="nElectron")
        
        self.out.branch("Jet_isGoodJet",    "I", lenVar="nJet")
        self.out.branch("Jet_isGoodBJet",   "I", lenVar="nJet")
        self.out.branch("Jet_leptonClean",  "I", lenVar="nJet")
        self.out.branch("Jet_WIdx",  "I", lenVar="nJet")

        self.out.branch("FatJet_isGoodJet", "I", lenVar="nFatJet")
        self.out.branch("FatJet_fromW", "I", lenVar="nFatJet")
        self.out.branch("FatJet_fromH", "I", lenVar="nFatJet")
        self.out.branch("FatJet_leptonClean",  "I", lenVar="nFatJet")


        #MIRANDA'S ADD-ONS
        self.out.branch("Tau_isVeto",       "F", lenVar="nTau")
        self.out.branch("nVetoTau",         "I")
        self.out.branch("IsoTrack_isVeto",  "F", lenVar="nIsoTrack")
        self.out.branch("nVetoIsoTrack",    "I")

        # Counter for good b-tags
        self.out.branch("nLepton",      "I")
        self.out.branch("nVetoLepton",  "I")
        self.out.branch("nGoodJet",     "I")
        self.out.branch("nGoodBTag",    "I")

        self.out.branch("GoodJet_pt"    ,"F",  lenVar="nGoodJet")
        self.out.branch("GoodJet_eta"   ,"F",  lenVar="nGoodJet")
        self.out.branch("GoodJet_phi"   ,"F",  lenVar="nGoodJet")

        self.out.branch("BJet_pt"    ,"F",  lenVar="nBJet")
        self.out.branch("BJet_eta"   ,"F",  lenVar="nBJet")
        self.out.branch("BJet_phi"   ,"F",  lenVar="nBJet")

        self.out.branch("BFatJet_pt"    ,"F",  lenVar="nBFatJet")
        self.out.branch("BFatJet_eta"   ,"F",  lenVar="nBFatJet")
        self.out.branch("BFatJet_phi"   ,"F",  lenVar="nBFatJet")

        self.out.branch("GoodFatJet_pt"    ,"F",  lenVar="nGoodFatJet")
        self.out.branch("GoodFatJet_eta"   ,"F",  lenVar="nGoodFatJet")
        self.out.branch("GoodFatJet_phi"   ,"F",  lenVar="nGoodFatJet")

        self.out.branch("W_pt"              ,"F",  lenVar="nW")
        self.out.branch("W_eta"             ,"F",  lenVar="nW")
        self.out.branch("W_phi"             ,"F",  lenVar="nW")

        self.out.branch("H_pt"              ,"F",  lenVar="nH")
        self.out.branch("H_eta"             ,"F",  lenVar="nH")
        self.out.branch("H_phi"             ,"F",  lenVar="nH")
        
        self.out.branch("bb_pt"               ,"F",  lenVar="nbb")
        self.out.branch("bb_eta"               ,"F",  lenVar="nbb")
        self.out.branch("bb_phi"               ,"F",  lenVar="nbb")
        self.out.branch("bb_mass"               ,"F",  lenVar="nbb")
        self.out.branch("bb_mct"               ,"F",  lenVar="nbb")

        self.out.branch("jj_pt"               ,"F", lenVar="njj")
        self.out.branch("jj_eta"               ,"F", lenVar="njj")
        self.out.branch("jj_phi"               ,"F", lenVar="njj")
        self.out.branch("jj_mass"               ,"F", lenVar="njj")
        self.out.branch("jj_mct"               ,"F", lenVar="njj")

        self.out.branch("WHptMET"             ,"F",  lenVar="nWH")

        self.out.branch("ak4matchedH_pt"       ,"F",  lenVar="njmH")
        self.out.branch("ak4matchedH_eta"       ,"F",  lenVar="njmH")
        self.out.branch("ak4matchedH_phi"       ,"F",  lenVar="njmH")

        self.out.branch("bjetmatchedH_pt"       ,"F",  lenVar="nbmH")
        self.out.branch("bjetmatchedH_eta"       ,"F",  lenVar="nbmH")
        self.out.branch("bjetmatchedH_phi"       ,"F",  lenVar="nbmH")

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

    def isGoodFatJet(self, fatjet):
        return (fatjet.pt > 200 and fatjet.jetId>0)

    def isGoodBFatJet(self, fatjet):
        if self.year == 2016:
            threshold = 0.6321
        if self.year == 2017:
            threshold = 0.4941
        if self.year == 2018:
            threshold = 0.4184
        return(self.isGoodFatJet(fatjet) and fatjet.btagDeepB > threshold)

    def isVetoMuon(self, muon):
        return (muon.looseId and muon.pt>5 and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.2 and abs(muon.dxy)<0.1 and abs(muon.dz)<0.5)

    def isVetoElectron(self, electron):
        return (electron.cutBased>0 and electron.miniPFRelIso_all < 0.2)

    def isTightMuon(self, muon):
        return (muon.pt > 25 and muon.mediumId and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.1)

    def isTightElectron(self, electron):
        return (electron.pt > 30 and electron.cutBased >= 3 and abs(electron.eta) < 2.4 and electron.miniPFRelIso_all < 0.1)# and electron.sip3d < 4.0 and abs(electron.dxy) < 0.05 and abs(electron.dz) < 0.1)

    def isVetoTau(self, tau):
        return(tau.pt > 20 and abs(tau.eta) < 2.3 and tau.idDecayMode and tau.idMVAoldDM2017v1 == 8)

    def isVetoIsoTrack(self, isotrack):
        return(isotrack.pt > 10 and abs(isotrack.eta) < 2.4 and (isotrack.miniPFRelIso_all < (0.1*isotrack.pt) or isotrack.miniPFRelIso_all < 6))

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

    #Don't think we actually need any of the following, but keeping for now         
    def isGoodNonBJet(self,jet):
        if self.year == 2016:
            threshold = 0.6321
        if self.year == 2017:
            threshold = 0.4941
        if self.year == 2018:
            threshold = 0.4184
        return (self.isGoodJet(jet) and jet.btagDeepB <= threshold)

    def isFatJetfromW(self, fatjet):
        return(self.isGoodFatJet(fatjet) and fatjet.deepTagMD_WvsQCD > 0.9)

    def isFatJetfromH(self, fatjet):
        if self.year == 2016:
            threshold = 0.8945
        if self.year == 2017:
            threshold = 0.8695
        if self.year == 2018:
            threshold = 0.8365
        return(self.isGoodFatJet(fatjet) and fatjet.deepTagMD_HbbvsQCD > threshold)

    def MCT2(self, b1, b2):
      return 2*b1['pt']*b2['pt']*(1+math.cos(self.deltaPhi(b1['phi'],b2['phi'])))

    def MCT(self, b1, b2):
      return math.sqrt(self.MCT2(b1,b2))

    def MT2(self, b1, metpt, metphi):
      return 2*b1['pt']*metpt*(1+math.cos(self.deltaPhi(b1['phi'],metphi)))

    def MT(self, b1, metpt, metphi):
      return math.sqrt(self.MT2(b1,metpt, metphi))

    def Mjj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1['pt'], b1['eta'], b1['phi'], 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2['pt'], b2['eta'], b2['phi'], 0)
      return (bjet1 + bjet2).M()

    def Ptjj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1['pt'], b1['eta'], b1['phi'], 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2['pt'], b2['eta'], b2['phi'], 0)
      return (bjet1 + bjet2).Pt()

    def Etajj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1['pt'], b1['eta'], b1['phi'], 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2['pt'], b2['eta'], b2['phi'], 0)
      return (bjet1 + bjet2).Eta()

    def Phijj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1['pt'], b1['eta'], b1['phi'], 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2['pt'], b2['eta'], b2['phi'], 0)
      return (bjet1 + bjet2).Phi()


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

        # tight lepton collection, will be sorted by pt
        leptons     = []
        vleptons    = []

        isTightMuon = []
        isVetoMuon = []
        for i,mu in enumerate(muons):
            mu.isTight  = self.isTightMuon(mu)
            mu.isVeto   = self.isVetoMuon(mu)
            isTightMuon.append(mu.isTight)
            isVetoMuon.append(mu.isVeto)
            if self.isTightMuon(mu):
                leptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1, 'mass':mu.mass})
            if self.isVetoMuon(mu):
                vleptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1, 'mass':mu.mass})


        isTightElectron = []
        isVetoElectron = []
        for i,el in enumerate(electrons):
            el.isTight  = self.isTightElectron(el)
            el.isVeto   = self.isVetoElectron(el)
            isTightElectron.append(el.isTight)
            isVetoElectron.append(el.isVeto)
            if self.isTightElectron(el):
                leptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i, 'mass':el.mass})
            if self.isVetoElectron(el):
                vleptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i, 'mass':el.mass})


        isVetoTau   = []
        for t in taus:
            t.cleanmask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if self.deltaR(t, p) < 0.4:
                        t.cleanmask = 0
            isVetoTau.append(1 if (self.isVetoTau(t) and t.cleanmask) else 0)

            
        isVetoIsoTrack  = []
        for t in isotracks:
            t.cleanmask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if self.deltaR(t, p) > 0.4:
                        t.cleanmask = 0
            isVetoIsoTrack.append(1 if (self.isVetoIsoTrack(t)and t.cleanmask) else 0)


        cleanMaskV  = []
        isGoodJet   = []
        isGoodBJet  = []
        isGoodJetAll = []

        fwdjets     = []
        jets_out    = []
        bjets       = []
        nonbjets    = []
        alljets = []
        WIdx = []

        for j in jets:

            j.cleanMask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if p.isVeto:
                        if self.deltaR(j, p) < 0.4:
                            j.cleanMask = 0

            isGoodJet.append(1 if (self.isGoodJet(j) and j.cleanMask) else 0)
            isGoodJetAll.append(1 if (self.isFwdJet(j) and j.cleanMask) else 0)
            isGoodBJet.append(1 if (self.isGoodBJet(j) and j.cleanMask) else 0)
            
            cleanMaskV.append(j.cleanMask)


            # now get some info from the GenJets
            j.WIdx = -1
            if not self.isData:
                if j.genJetIdx>-1 and j.genJetIdx<len(genjets): # apparently lowest pt genjets are not stored...
                    j.WIdx = genjets[j.genJetIdx].WIdx
            WIdx.append(j.WIdx)

            alljets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx})

            # Fill the other collections
            if j.cleanMask:

                if self.isGoodJet(j):
                    jets_out.append({'pt':j.pt_nom, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx, 'mass':j.mass_nom})
                    
                    if self.isGoodBJet(j):
                        bjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})
                    else:
                        nonbjets.append({'pt':j.pt_nom, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx, 'mass':j.mass_nom})
                
                if self.isFwdJet(j):
                    fwdjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})

        cleanMaskW  = []
        isGoodFatJet   = []
        isFatJetfromW  = []
        isFatJetfromH = []
        
        goodfatjets = []
        Ws = []
        Hs = []
        bfatjets = []

        for f in fatjets:
            f.cleanMask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if p.isVeto:
                        if self.deltaR(f, p) < 0.4:
                            f.cleanMask = 0

            isGoodFatJet.append(1 if (self.isGoodFatJet(f) and f.cleanMask) else 0)
            isFatJetfromW.append(1 if (self.isFatJetfromW(f) and f.cleanMask) else 0)
            isFatJetfromH.append(1 if (self.isFatJetfromH(f) and f.cleanMask) else 0)
            
            cleanMaskW.append(f.cleanMask)

            # Fill the other collections
            if f.cleanMask:

                if self.isGoodFatJet(f):
                    goodfatjets.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
                    
                    if self.isGoodBFatJet(f):
                        bfatjets.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})

                    if self.isFatJetfromH(f):
                        Hs.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
                    elif self.isFatJetfromW(f):
                        Ws.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
        
        #JJs
        jj = []
        dphijmet = []
        maxptjj = -9999
        
        #make loop to calculate jj with highest pt and calculate variables
        for i in nonbjets:
            for j in nonbjets:
                if i != j:
                    if ( (i['pt'] + j['pt']) > maxptjj):
                        maxptjj = j['pt'] + i['pt']
                        leadjet = i
                        subjet = j
                        dphilead = self.deltaPhi(leadjet['phi'], met_phi)
                        dphisub = self.deltaPhi(subjet['phi'], met_phi)
                        mctjj = self.MCT(leadjet,subjet)
                        mjj = self.Mjj(leadjet,subjet)
                        ptjj = self.Ptjj(leadjet,subjet)
                        etajj = self.Etajj(leadjet,subjet)
                        phijj = self.Phijj(leadjet,subjet)

        #if there are two jets, append variables
        if maxptjj != -9999:
            dphijmet.append(dphilead)
            dphijmet.append(dphisub)
            mindphijmet = min(dphijmet)
            jj.append({'pt': ptjj, 'eta': etajj, 'phi': phijj, 'mass': mjj, "mct": mctjj, "mindphijmet": mindphijmet})

        #BBs

        bb = []
        dphibbmet = []
        mtbmet = []
        maxptbb = -9999
        
        for b in bjets:
            for j in bjets:
                if b != j:
                    if ( (b['pt'] + j['pt']) > maxptbb):
                        maxptbb = b['pt'] + j['pt']
                        leadjet = b
                        subjet = j
                        mtblead = self.MT(leadjet,met_pt,met_phi)
                        mtbsub = self.MT(subjet, met_pt,met_phi)
                        mct = self.MCT(leadjet,subjet)
                        Mbb = self.Mjj(leadjet,subjet)
                        ptbb = self.Ptjj(leadjet,subjet)
                        etabb = self.Etajj(leadjet,subjet)
                        phibb = self.Phijj(leadjet,subjet)
                        dphibb = self.deltaPhi(phibb, met_phi)
               

        #if there are two jets, append variables
        if maxptbb != -9999:
            mtbmet.append(mtblead)
            mtbmet.append(mtbsub)
            minmtbmet = min(mtbmet)
            bb.append({'pt': ptbb, 'eta': etabb, 'phi': phibb, 'mass': Mbb, "mct": mct, "dphibbmet": dphibb, "minmtbmet": minmtbmet})


        #BBJJs

        bbjj = []

        if (maxptjj != -9999) and (maxptbb != -9999): 
            dphibbjj = self.deltaPhi(phibb,phijj)
            bbjj.append({"dphibbjj": dphibbjj})
              

        #WH
        WH = []
        maxptwh = -9999
        for H in Hs:
            for W in Ws:
                if( (H['pt'] + W['pt']) > maxptwh):
                    if H != W:
                        maxptwh = H['pt'] + W['pt']
                        leadjet = H
                        leadojet = W
                        #print "%f"%(maxptwh)
                        whpt = self.Ptjj(leadjet, leadojet)
                        whptmet = whpt/met_pt
                        dphiWH = self.deltaPhi(leadjet['phi'],leadojet['phi'])
                        
        if maxptwh != -9999:
            WH.append({"WHptMET": whptmet, "dphiWH": dphiWH})

        #MT
        maxwpt = -9999
        maxhpt = -9999
        maxfjpt = -9999
        nW = 0
        nH = 0
        for W in Ws:
            nW += 1
            if W['pt'] > maxwpt:
                maxwpt = W['pt']
                leadW = W
                mtwMET = self.MT(leadW, met_pt, met_phi)
        for H in Hs:
            nH += 1
            if H['pt'] > maxhpt:
                maxhpt = H['pt']
                leadH = j 
                mthMET = self.MT(leadH, met_pt, met_phi)
        for j in goodfatjets:
            if nW == 0 and nH == 0:
                if j['pt'] > maxfjpt:
                    maxfjpt = j['pt']
                    leadfj = j
                    mtfjMET = self.MT(leadfj, met_pt, met_phi)
        
        leadW = []
        leadH = []
        leadFJ = []
        if maxwpt != -9999 and nW != 0:
            leadW.append({"mtwMET": mtwMET})
        if maxhpt != -9999 and nH != 0:
            leadH.append({"mthMET": mthMET})
        if maxfjpt != -9999 and nW == 0 and nH == 0:
            leadFJ.append({"mtfjMET": mtfjMET})


        #MATCHING AK8 AND AK4
        #HIGGS
        '''jmatchH = []
        bmatchH = []
        for H in Hs:
            for j in jets_out:
                if (self.deltaR(H,j) < 0.4):
                    jmatchH.append({'pt':j['pt'], 'eta':j['eta'], 'phi':j['phi']})
            for b in bjets:
                if (self.deltaR(H,b) < 0.4):
                    bmatchH.append({'pt':b['pt'], 'eta':b['eta'], 'phi':b['phi']})'''
                       

        # make sure the jets are properly sorted. they _should_ be sorted, but this can change once we reapply JECs if necessary
        bjets       = sorted(bjets, key = lambda i: i['pt'], reverse=True)
        nonbjets    = sorted(nonbjets, key = lambda i: i['pt'], reverse=True)
        fwdjets     = sorted(fwdjets, key = lambda i: i['pt'], reverse=True) # all jets, including forward ones
        jets_out    = sorted(jets_out, key = lambda i: i['pt'], reverse=True)
        Ws          = sorted(Ws, key = lambda i: i['pt'], reverse=True)
        Hs          = sorted(Hs, key = lambda i: i['pt'], reverse=True)

        # calculate MT, Mlb, Mlb_max, M_jj_b1, M_jj_b2 etc
        Mlb_00 = -99 # leading lepton, first b-jet
        Mlb_01 = -99 # leading lepton, second b-jet
        Mlb_10 = -99 # trailing lepton, first b-jet
        Mlb_11 = -99 # trailing lepton, second b-jet

        #if len(jets)>0:
        #    if len(leptons)>0:
        #        Mlb_00 = 
        #
        #if len(leptons)>1

            


        self.out.fillBranch("Muon_isTight",     isTightMuon)
        self.out.fillBranch("Muon_isVeto",      isVetoMuon)
        self.out.fillBranch("Electron_isTight", isTightElectron)
        self.out.fillBranch("Electron_isVeto",  isVetoElectron)
        self.out.fillBranch("Jet_leptonClean",   cleanMaskV)
        self.out.fillBranch("Jet_isGoodJet",    isGoodJet)
        #self.out.fillBranch("Jet_isGoodJetAll", isGoodJetAll)
        self.out.fillBranch("Jet_isGoodBJet",   isGoodBJet)
        self.out.fillBranch("Jet_WIdx",   WIdx)
        self.out.fillBranch("FatJet_leptonClean",   cleanMaskW)
        self.out.fillBranch("FatJet_isGoodJet",     isGoodFatJet)
        self.out.fillBranch("FatJet_fromW",         isFatJetfromW)
        self.out.fillBranch("FatJet_fromH",         isFatJetfromH)
        self.out.fillBranch("nGoodBTag",        sum(isGoodBJet))
        self.out.fillBranch("nGoodJet",         sum(isGoodJet))

        self.out.fillBranch("Tau_isVeto",           isVetoTau)
        self.out.fillBranch("nVetoTau",             sum(isVetoTau))
        self.out.fillBranch("IsoTrack_isVeto",      isVetoIsoTrack)
        self.out.fillBranch("nVetoIsoTrack",        sum(isVetoIsoTrack))

        
        '''jmatchH_pd = pd.DataFrame(jmatchH)
        self.out.fillBranch("njmH",             len(jmatchH_pd))
        if len(jmatchH_pd)>0:
            self.out.fillBranch("ak4matchedH_pt",    jmatchH_pd['pt'].tolist())
            self.out.fillBranch("ak4matchedH_eta",   jmatchH_pd['eta'].tolist())
            self.out.fillBranch("ak4matchedH_phi",   jmatchH_pd['phi'].tolist())
    
        bmatchH_pd = pd.DataFrame(bmatchH)
        self.out.fillBranch("nbmH",             len(bmatchH_pd))
        if len(bmatchH_pd)>0:
            self.out.fillBranch("bjetmatchedH_pt",    bmatchH_pd['pt'].tolist())
            self.out.fillBranch("bjetmatchedH_eta",   bmatchH_pd['eta'].tolist())
            self.out.fillBranch("bjetmatchedH_phi",   bmatchH_pd['phi'].tolist())'''

        goodjets_pd = pd.DataFrame(jets_out)
        self.out.fillBranch("nGoodJet",          len(goodjets_pd) )
        if len(goodjets_pd)>0:
            self.out.fillBranch("GoodJet_pt",        goodjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("GoodJet_eta",       goodjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("GoodJet_phi",       goodjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        bjets_pd = pd.DataFrame(bjets)
        self.out.fillBranch("nBJet",          len(bjets_pd) )
        if len(bjets_pd)>0:
            self.out.fillBranch("BJet_pt",        bjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("BJet_eta",       bjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("BJet_phi",       bjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        bfatjets_pd = pd.DataFrame(bfatjets)
        self.out.fillBranch("nBFatJet",          len(bfatjets_pd) )
        if len(bfatjets_pd)>0:
            self.out.fillBranch("BFatJet_pt",     bfatjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("BFatJet_eta",    bfatjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("BFatJet_phi",    bfatjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        bb_pd = pd.DataFrame(bb)
        self.out.fillBranch("nbb",          len(bb_pd) )
        if len(bb_pd)>0:
            self.out.fillBranch("bb_pt",        bb_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("bb_eta",       bb_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("bb_phi",       bb_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("bb_mass",       bb_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("bb_mct",         bb_pd.sort_values(by='pt', ascending=False)['mct'].tolist())

        jj_pd = pd.DataFrame(jj)
        self.out.fillBranch("njj",          len(jj_pd) )
        if len(jj_pd)>0:
            self.out.fillBranch("jj_pt",        jj_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("jj_eta",       jj_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("jj_phi",       jj_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("jj_mass",       jj_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("jj_mct",         jj_pd.sort_values(by='pt', ascending=False)['mct'].tolist())

        WH_pd = pd.DataFrame(WH)
        self.out.fillBranch("nWH",              len(WH_pd))
        if len(WH_pd) > 0:
            self.out.fillBranch("WHptMET",      WH_pd['WHptMET'].tolist() )

        goodfatjets_pd = pd.DataFrame(goodfatjets)
        self.out.fillBranch("nGoodFatJet",          len(goodfatjets_pd) )
        if len(goodfatjets_pd)>0:
            self.out.fillBranch("GoodFatJet_pt",        goodfatjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("GoodFatJet_eta",       goodfatjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist())
            self.out.fillBranch("GoodFatJet_phi",       goodfatjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist())

        Ws_pd = pd.DataFrame(Ws)
        self.out.fillBranch("nW",          len(Ws_pd) )
        if len(Ws_pd)>0:
            self.out.fillBranch("W_pt",       Ws_pd.sort_values(by='pt', ascending=False)['pt'].tolist())
            self.out.fillBranch("W_eta",      Ws_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("W_phi",      Ws_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        Hs_pd = pd.DataFrame(Hs)
        self.out.fillBranch("nH",          len(Hs_pd) )
        if len(Hs_pd)>0:
            self.out.fillBranch("H_pt",       Hs_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("H_eta",      Hs_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("H_phi",      Hs_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )


        # make pandas dataframe out of list
        leptons_pd = pd.DataFrame(leptons)

        self.out.fillBranch("nVetoLepton",          len(vleptons) )
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
