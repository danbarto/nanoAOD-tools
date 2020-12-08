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
#from PhysicsTools.NanoAODTools.postprocessing.framework.mt2Calculator import mt2Calculator


def hasBit(value,bit):
  """Check if i'th bit is set to 1, i.e. binary of 2^(i-1),
  from the right to the left, starting from position i=0."""
  # https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#GenPart
  # Gen status flags, stored bitwise, are:
  #    0: isPrompt,                          8: fromHardProcess,
  #    1: isDecayedLeptonHadron,             9: isHardProcessTauDecayProduct,
  #    2: isTauDecayProduct,                10: isDirectHardProcessTauDecayProduct,
  #    3: isPromptTauDecayProduct,          11: fromHardProcessBeforeFSR,
  #    4: isDirectTauDecayProduct,          12: isFirstCopy,
  #    5: isDirectPromptTauDecayProduct,    13: isLastCopy,
  #    6: isDirectHadronDecayProduct,       14: isLastCopyBeforeFSR
  #    7: isHardProcess,
  ###return bin(value)[-bit-1]=='1'
  ###return format(value,'b').zfill(bit+1)[-bit-1]=='1'
  return (value & (1 << bit))>0

class GenAnalyzer(Module):

    def __init__(self,isW,isWExt):
        self.isW = isW
        self.isWExt = isWExt
        print isW
        print isWExt
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        ## Define a first minimum set of objects needed for the analysis
        #FIXME objects should have cross-reference to full collection

        self.out.branch("GenW_pt", "F",         lenVar="nGenW")
        self.out.branch("GenW_eta", "F",        lenVar="nGenW")
        self.out.branch("GenW_phi", "F",        lenVar="nGenW")
        self.out.branch("GenW_pdgId", "I",      lenVar="nGeGW")

        self.out.branch("GenH_pt", "F",         lenVar="nGenH")
        self.out.branch("GenH_eta", "F",        lenVar="nGenH")
        self.out.branch("GenH_phi", "F",        lenVar="nGenH")
        self.out.branch("GenH_pdgId", "I",      lenVar="nGenH")

        self.out.branch("GenBfromH_pt", "F",       lenVar="nGenBfromH")
        self.out.branch("GenBfromH_eta", "F",      lenVar="nGenBfromH")
        self.out.branch("GenBfromH_phi", "F",      lenVar="nGenBfromH")
        self.out.branch("GenBfromH_pdgId", "I",    lenVar="nGenBfromH")
        self.out.branch("GenBfromH_fromH", "I",    lenVar="nGenBfromH")

        self.out.branch("GenQfromW_pt", "F",       lenVar="nGenQfromW")
        self.out.branch("GenQfromW_eta", "F",      lenVar="nGenQfromW")
        self.out.branch("GenQfromW_phi", "F",      lenVar="nGenQfromW")
        self.out.branch("GenQfromW_pdgId", "I",    lenVar="nGenQfromW")
        self.out.branch("GenQfromW_fromW", "I",    lenVar="nGenQfromW")

        #Re-adding genLepton and nLep info
        self.out.branch("GenL_pt", "F",       lenVar="nGenL")
        self.out.branch("GenL_eta", "F",      lenVar="nGenL")
        self.out.branch("GenL_phi", "F",      lenVar="nGenL")
        self.out.branch("GenL_pdgId", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromTop", "I",  lenVar="nGenL")
        self.out.branch("GenL_fromTau", "I",  lenVar="nGenL")
        self.out.branch("GenL_fromZ", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromW", "I",    lenVar="nGenL")

        self.out.branch("nGenLepFromTop",     "I")
        self.out.branch("nGenLepFromTau",    "I")
        self.out.branch("nGenLepFromW",    "I")
        self.out.branch("nGenLepFromZ",    "I")
        self.out.branch("nGenTau",    "I")

        #etc variables
        self.out.branch("stitch", "B")

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

    def hasAncestor(self, p, ancestorPdg, genParts):
        motherIdx = p.genPartIdxMother
        while motherIdx>0:
            if (abs(genParts[motherIdx].pdgId) == ancestorPdg): return True
            motherIdx = genParts[motherIdx].genPartIdxMother
        return False
    
    def MCT2(self, b1, b2):
      return 2*b1.pt*b2.pt*(1+math.cos(self.deltaPhi(b1.phi,b2.phi)))

    def MCT(self, b1, b2):
      return math.sqrt(self.MCT2(b1,b2))

    def Mbb(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).M()

    def Ptjj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).Pt()

    def getAncestorIdx(self, p, ancestorPdg, genParts):
        motherIdx = p.genPartIdxMother
        while motherIdx>=0:
            if (abs(genParts[motherIdx].pdgId) == ancestorPdg): return motherIdx
            motherIdx = genParts[motherIdx].genPartIdxMother
        return -1

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

#        print 'event %d *dr phil voice* shut the hell up bitch go kill yourself'%( event.event)

        #MET
        met_pt  = event.MET_pt
        met_phi = event.MET_phi
        
        #Particles
        GenParts    = Collection(event, "GenPart")
        GenAK8Jets  = Collection(event, "GenJetAK8")
        GenJets     = Collection(event, "GenJet")

        tops = [ p for p in GenParts if (abs(p.pdgId)==6 and hasBit(p.statusFlags,13) ) ] # last copy ts

        Ws = [ p for p in GenParts if (abs(p.pdgId)==24 and hasBit(p.statusFlags,13) ) ] # last copy Ws

        Hs = [ p for p in GenParts if (abs(p.pdgId)==25 and hasBit(p.statusFlags,13) )] # last copy Hs

        leptons = [ p for p in GenParts if ((abs(p.pdgId)==11 or abs(p.pdgId)==13) and hasBit(p.statusFlags,13) and (hasBit(p.statusFlags,0) or hasBit(p.statusFlags, 2)) ) ]
        neutrinos = [ p for p in GenParts if ((abs(p.pdgId)==12 or abs(p.pdgId)==14 or abs(p.pdgId)==16) and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0) and p.status==1 ) ]
        taus    = [ p for p in GenParts if ((abs(p.pdgId)==15) and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0) ) ]

        quarks  = [ p for p in GenParts if ((abs(p.pdgId)<5)  and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0)  and self.hasAncestor(p, 24, GenParts) and (GenParts[p.genPartIdxMother].pdgId==p.pdgId or abs(GenParts[p.genPartIdxMother].pdgId)==24)) ] # status 71, 52, 51?

        bs = [ p for p in GenParts if (abs(p.pdgId) == 5 and hasBit(p.statusFlags, 0) and hasBit(p.statusFlags, 13)) and self.hasAncestor(p, 25, GenParts) ]  #hard scatter bs


        for lep in leptons:
            lep.fromTop = ( 1 if self.hasAncestor(lep, 6, GenParts) else 0 )
            lep.fromTau = ( 1 if self.hasAncestor(lep, 15, GenParts) else 0 )
            lep.fromZ = ( 1 if self.hasAncestor(lep, 23, GenParts) else 0 )
            lep.fromW = ( 1 if self.hasAncestor(lep, 24, GenParts) else 0 )

        stitch = True
        if self.isW:
            sumPt = 0
            for nu in neutrinos:
                sumPt = sumPt + nu.pt
            if self.isWExt and sumPt <200: stitch = False
            if not self.isWExt and sumPt >200: stitch = False


        #--------------------------------------FILL---BRANCHES-----------------------------------------#
        
        #some jet stuff


        self.out.fillBranch("nGenLepFromTop", sum( [ l.fromTop for l in leptons ] ) )
        self.out.fillBranch("nGenLepFromTau", sum( [ l.fromTau for l in leptons ] ) )
        self.out.fillBranch("nGenLepFromW",   sum( [ l.fromW for l in leptons ] ) )
        self.out.fillBranch("nGenLepFromZ",   sum( [ l.fromZ for l in leptons ] ) )
        self.out.fillBranch("nGenTau",        len(taus) )

        self.out.fillBranch("nGenL",          len(leptons) )
        if len(leptons)>0:
            self.out.fillBranch("GenL_pt",        [ l.pt for l in leptons ])
            self.out.fillBranch("GenL_eta",       [ l.eta for l in leptons ])
            self.out.fillBranch("GenL_phi",       [ l.phi for l in leptons ])
            self.out.fillBranch("GenL_pdgId",     [ l.pdgId for l in leptons ])
            self.out.fillBranch("GenL_fromTop",   [ l.fromTop for l in leptons ])
            self.out.fillBranch("GenL_fromTau",   [ l.fromTau for l in leptons ])
            self.out.fillBranch("GenL_fromW",     [ l.fromW for l in leptons ])
            self.out.fillBranch("GenL_fromZ",     [ l.fromZ for l in leptons ])
            

        self.out.fillBranch("nGenW",          len(Ws) )
        if len(Ws)>0:
            self.out.fillBranch("GenW_pt",        [ p.pt    for p in Ws ])
            self.out.fillBranch("GenW_eta",       [ p.eta   for p in Ws ])
            self.out.fillBranch("GenW_phi",       [ p.phi   for p in Ws ])
            self.out.fillBranch("GenW_pdgId",     [ p.pdgId for p in Ws ])
        
        self.out.fillBranch("nGenH",          len(Hs) )
        if len(Hs)>0:
            self.out.fillBranch("GenH_pt",        [ p.pt    for p in Hs ])
            self.out.fillBranch("GenH_eta",       [ p.eta   for p in Hs ])
            self.out.fillBranch("GenH_phi",       [ p.phi   for p in Hs ])
            self.out.fillBranch("GenH_pdgId",     [ p.pdgId for p in Hs ])

        self.out.fillBranch("nGenBfromH",          len(bs) )
        if len(bs)>0:
            self.out.fillBranch("GenBfromH_pt",        [ p.pt    for p in bs ])
            self.out.fillBranch("GenBfromH_eta",       [ p.eta   for p in bs ])
            self.out.fillBranch("GenBfromH_phi",       [ p.phi   for p in bs ])
            self.out.fillBranch("GenBfromH_pdgId",     [ p.pdgId for p in bs ])

        self.out.fillBranch("nGenQfromW",          len(quarks) )
        if len(quarks)>0:
            self.out.fillBranch("GenQfromW_pt",        [ p.pt    for p in quarks ])
            self.out.fillBranch("GenQfromW_eta",       [ p.eta   for p in quarks ])
            self.out.fillBranch("GenQfromW_phi",       [ p.phi   for p in quarks ])
            self.out.fillBranch("GenQfromW_pdgId",     [ p.pdgId for p in quarks ])

        self.out.fillBranch("stitch",     stitch  )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genAnalyzer = lambda  isW, isWExt : GenAnalyzer( isW=isW, isWExt=isWExt )
