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

    def __init__(self,isW,isWExt,isT,isTExt,year):
        self.isW = isW
        self.isWExt = isWExt
        self.isT = isT
        self.isTExt = isTExt
        self.year = year
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        ## Define a first minimum set of objects needed for the analysis
        #FIXME objects should have cross-reference to full collection

        self.out.branch("GenW_pt", "F",         lenVar="ngenW")
        self.out.branch("GenW_eta", "F",        lenVar="ngenW")
        self.out.branch("GenW_phi", "F",        lenVar="ngenW")
        self.out.branch("GenW_pdgId", "I",    lenVar="ngenW")

        self.out.branch("GenH_pt", "F",         lenVar="ngenH")
        self.out.branch("GenH_eta", "F",        lenVar="ngenH")
        self.out.branch("GenH_phi", "F",        lenVar="ngenH")
        self.out.branch("GenH_pdgId", "I",    lenVar="ngenH")

        self.out.branch("Genb_pt", "F",         lenVar="ngenb")
        self.out.branch("Genb_eta", "F",        lenVar="ngenb")
        self.out.branch("Genb_phi", "F",        lenVar="ngenb")
        self.out.branch("Genb_pdgId", "I",    lenVar="ngenb")
        self.out.branch("Genb_fromH", "I",    lenVar="ngenb")

        self.out.branch("Genj_pt", "F",         lenVar="ngenj")
        self.out.branch("Genj_eta", "F",        lenVar="ngenj")
        self.out.branch("Genj_phi", "F",        lenVar="ngenj")
        self.out.branch("Genj_pdgId", "I",    lenVar="ngen")
        self.out.branch("Genj_fromW", "I",    lenVar="ngenj")

        #Re-adding genLepton and nLep info
        self.out.branch("GenL_pt", "F",         lenVar="nGenL")
        self.out.branch("GenL_eta", "F",        lenVar="nGenL")
        self.out.branch("GenL_phi", "F",        lenVar="nGenL")
        self.out.branch("GenL_pdgId", "I",      lenVar="nGenL")
        self.out.branch("GenL_fromTop", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromTau", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromZ", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromW", "I",    lenVar="nGenL")

        self.out.branch("GenJet_fromW", "I",    lenVar="nGenJet")
        self.out.branch("GenJet_genPartWIdx", "I",    lenVar="nGenJet")
        self.out.branch("GenJet_WIdx", "I",    lenVar="nGenJet")

        self.out.branch("Spectator_pt", "F",         lenVar="nSpectator")
        self.out.branch("Spectator_eta", "F",        lenVar="nSpectator")
        self.out.branch("Spectator_phi", "F",        lenVar="nSpectator")
        self.out.branch("Spectator_pdgId", "I",      lenVar="nSpectator")

        self.out.branch("Scatter_pt", "F",         lenVar="nScatter")
        self.out.branch("Scatter_eta", "F",        lenVar="nScatter")
        self.out.branch("Scatter_phi", "F",        lenVar="nScatter")
        self.out.branch("Scatter_pdgId", "I",      lenVar="nScatter")

        self.out.branch("W_pt", "F",         lenVar="nW")
        self.out.branch("W_eta", "F",        lenVar="nW")
        self.out.branch("W_phi", "F",        lenVar="nW")
        self.out.branch("W_mass", "F",    lenVar="nW")
        self.out.branch("W_pdgId", "I",    lenVar="nW")
        self.out.branch("W_genPartIdx", "I",    lenVar="nW")

        self.out.branch("Top_pt", "F",         lenVar="nTop")
        self.out.branch("Top_eta", "F",        lenVar="nTop")
        self.out.branch("Top_phi", "F",        lenVar="nTop")
        self.out.branch("Top_pdgId", "I",        lenVar="nTop")

        # Counter for good b-tags
        self.out.branch("nLepFromTop",     "I")
        self.out.branch("nLepFromTau",    "I")
        self.out.branch("nLepFromW",    "I")
        self.out.branch("nGenTau",    "I")
        self.out.branch("nLepFromZ",    "I")

        #miranda time
        self.out.branch("GenDeltaRbb", "F",   lenVar="ngenbb")
        self.out.branch("GenMCT", "F",        lenVar="ngenbb")
        self.out.branch("GenMbb", "F",         lenVar="ngenbb")

        self.out.branch("GenDeltaRjj", "F",   lenVar="ngenjj")
           

#        self.out.branch("GenMT2_WH", "F",     lenVar="ngenWH")
#        self.out.branch("GenMT2_bbjj", "F",     lenVar="ngenmt2")
#        self.out.branch("GenMT2_bjjb", "F",     lenVar="ngenmt2")

        self.out.branch("Genleadb_pt", "F")
        self.out.branch("Genleadnonb_pt", "F")

        self.out.branch("Genak8j_pt", "F",     lenVar="ngenak8jet")
        self.out.branch("Genak8j_eta", "F",     lenVar="ngenak8jet")
        self.out.branch("Genak8j_phi", "F",     lenVar="ngenak8jet")

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
        genMET_pt  = event.GenMET_pt

        #Particles
        GenParts    = Collection(event, "GenPart")
        GenAK8Jets  = Collection(event, "GenJetAK8")
        GenJets     = Collection(event, "GenJet")

        tops = [ p for p in GenParts if (abs(p.pdgId)==6 and hasBit(p.statusFlags,13) ) ] # last copy ts
        #print len(tops)

        #Ws = [ (p, idx) for idx, p in enumerate(GenParts) if (abs(p.pdgId)==24 and hasBit(p.statusFlags,13) ) ] # last copy Ws
        Ws = []
        for idx, p in enumerate(GenParts):
            if abs(p.pdgId)==24 and hasBit(p.statusFlags,13):
                p.idx = idx
                Ws.append(p)
            

        leptons = [ p for p in GenParts if ((abs(p.pdgId)==11 or abs(p.pdgId)==13) and hasBit(p.statusFlags,13) and (hasBit(p.statusFlags,0) or hasBit(p.statusFlags, 2)) ) ]
        neutrinos = [ p for p in GenParts if ((abs(p.pdgId)==12 or abs(p.pdgId)==14 or abs(p.pdgId)==16) and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0) and p.status==1 ) ]
        taus    = [ p for p in GenParts if ((abs(p.pdgId)==15) and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0) ) ]

        quarks  = [ p for p in GenParts if ((abs(p.pdgId)<5)  and hasBit(p.statusFlags,13) and hasBit(p.statusFlags,0)  and self.hasAncestor(p, 24, GenParts) and (GenParts[p.genPartIdxMother].pdgId==p.pdgId or abs(GenParts[p.genPartIdxMother].pdgId)==24)) ] # status 71, 52, 51?
        #print '##', len(quarks), len(Ws), event.event, event.run, event.luminosityBlock
        #for q in quarks:
        #    print q.pdgId, q.status
        for quark in quarks:
            idx = self.getAncestorIdx(quark, 24, GenParts)
            quark.WIdx = idx


        js = [ p for p in GenParts if ((abs(p.pdgId) == 1 or abs(p.pdgId) == 2 or abs(p.pdgId) == 3 or abs(p.pdgId) == 4 or abs(p.pdgId) == 5 or abs(p.pdgId) == 6) and hasBit(p.statusFlags, 7))] 

        nonbjs = [ p for p in GenParts if ((abs(p.pdgId) == 1 or abs(p.pdgId) == 2 or abs(p.pdgId) == 3 or abs(p.pdgId) == 4 or abs(p.pdgId) == 6) and hasBit(p.statusFlags, 7))] 


        Hs = [ p for p in GenParts if (abs(p.pdgId)==25 and hasBit(p.statusFlags,13) )] # last copy Hs

        bs = [ p for p in GenParts if (abs(p.pdgId) == 5 and hasBit(p.statusFlags, 7))]  #hard scatter bs


        for lep in leptons:
            lep.fromTop = ( 1 if self.hasAncestor(lep, 6, GenParts) else 0 )
            lep.fromTau = ( 1 if self.hasAncestor(lep, 15, GenParts) else 0 )
            lep.fromZ = ( 1 if self.hasAncestor(lep, 23, GenParts) else 0 )
            lep.fromW = ( 1 if self.hasAncestor(lep, 24, GenParts) else 0 )

        for j in js:
          j.fromW = (1 if self.hasAncestor(j, 24, GenParts) else 0 )

        stitch = True
        if self.isW:
            sumPt = 0
            for nu in neutrinos:
                sumPt = sumPt + nu.pt
            if self.isWExt and sumPt <200: stitch = False
            if not self.isWExt and sumPt >200: stitch = False
        elif self.isT:
            if self.year == 2018:
                if self.isTExt and genMET_pt<80: stitch = False
                if not self.isTExt and genMET_pt>80: stitch = False
            else:
                if self.isTExt and genMET_pt<150: stitch = False
                if not self.isTExt and genMET_pt>150: stitch = False

        for jet in GenJets:
            jet.genPartWIdx = -1 # this is the reference to the W in the GenPart collection
            jet.WIdx = -1 # this is the reference to the W in the W collection
            matched = False
            for quark in quarks:
                if self.deltaR(jet, quark)<0.4:
                    matched = True
                    jet.genPartWIdx = quark.WIdx
                    break
            if matched:
                for idx, W in enumerate(Ws):
                    if jet.genPartWIdx == W.idx:
                        jet.WIdx = idx
            jet.fromW = matched

        #print len(tops), len(Ws), len(scatter), len(spectator), len(leptons)

        leadbj_pt = []

        for b in bs:
          b.fromH = ( 1 if self.hasAncestor(b, 25, GenParts) else 0 )
          leadbj_pt.append(b.pt)

        leadnonbj_pt = []

        for n in nonbjs:
          leadnonbj_pt.append(n.pt)
          

        deltaRbb = []
        deltaRjj = []
        Mct = []
        mbb = []

        ibb = 0
        
        for b,s in itertools.combinations(bs, 2):
          if (b.fromH == 1 and s.fromH == 1 and (b.genPartIdxMother == s.genPartIdxMother)):
            ibb += 1
            dr = self.deltaR(b, s)
            mct = self.MCT(b,s)
            Mbb = self.Mbb(b,s)
          else:
            continue
          deltaRbb.append(dr)
          Mct.append(mct)
          mbb.append(Mbb)
 
         
        ijj = 0

        for b,s in itertools.combinations(js, 2):
          if (b.fromW == 1 and s.fromW == 1 and (b.genPartIdxMother == s.genPartIdxMother)):
            ijj += 1
            dr = self.deltaR(b, s)
          else:
            continue
          deltaRjj.append(dr)

        WHptMET = []

        if (len(Hs) > 0) and (len(Ws) > 0):
          for H in Hs:
            for W in Ws:
              whpt = self.Ptjj(H,W)
              whptmet = whpt / met_pt
              WHptMET.append(whptmet)

        #MT2WH
            
#        mt2WH = []

#        mt2Calculator.setMet(met_pt, met_phi)

#        if (len(Ws) > 0) and (len(Hs) > 0):
#          for H in Hs:
#            for W in Ws:
#              mt2Calculator.setJet1(H.pt, H.eta, H.phi)
#              mt2Calculator.setJet2(W.pt, W.eta, W.phi)
#              mt2WH.append(mt2Calculator.mt2jj())
        
        #MT2bbjj

#        mt2bbjj =[]
#        mt2bjjb = []
        
#        imt2 = 0

#        if len(bs) > 1 and len(js) > 1:
#          for b,s in itertools.combinations(bs, 2):
#            for j,z in itertools.combinations(js, 2):
#              if (b.fromH == 1 and s.fromH == 1 and j.fromW == 1 and z.fromW == 1 and (b.genPartIdxMother == s.genPartIdxMother) and (j.genPartIdxMother == z.genPartIdxMother)):
#                imt2 += 1
#                mt2Calculator.setBJet1(b.pt, b.eta, b.phi)
#                mt2Calculator.setBJet2(s.pt, s.eta, s.phi)
#                mt2Calculator.setJet1(j.pt, j.eta, j.phi)
#                mt2Calculator.setJet2(z.pt, z.eta, z.phi)
#                mt2bbjj.append(mt2Calculator.mt2bbjj())
#                mt2bjjb.append(mt2Calculator.mt2bjjb())


        #NEW MT2

        #create list of jets with no overlap
        
        goodjets = []

        for j in GenJets:
          overlap = False
          for g in GenAK8Jets:
            deltar = self.deltaR(j, g)
            if (deltar < 0.4):
              overlap = True
              break
          if overlap:
            continue
          goodjets.append(j)

        for g in GenAK8Jets:
          goodjets.append(g)

        #sort list of jets by highest pt
        goodjets.sort(reverse = True, key = lambda x : x.pt)


        #--------------------------------------FILL---BRANCHES-----------------------------------------#
        
        #some jet stuff

        self.out.fillBranch("ngenbb",          ibb)
        if ibb > 0:
          self.out.fillBranch("GenDeltaRbb",     deltaRbb)
          self.out.fillBranch("GenMCT",          Mct)
          self.out.fillBranch("GenMbb",           mbb)

        self.out.fillBranch("ngenjj",          ijj)
        if ijj > 0:
          self.out.fillBranch("GenDeltaRjj",     deltaRjj)


        if len(leadbj_pt)>0:
          self.out.fillBranch("Genleadb_pt",    max(leadbj_pt))
        if len(leadnonbj_pt)>0:
          self.out.fillBranch("Genleadnonb_pt", max(leadnonbj_pt))

        self.out.fillBranch("nWH",             len(Ws)*len(Hs))
        if len(Ws) * len(Hs) > 0:
          self.out.fillBranch("WHptMET",     WHptMET)
        #MT2
#        self.out.fillBranch("ngenWH",          len(Hs)*len(Ws) )
#        if (len(Hs)*len(Ws) >0):
#          self.out.fillBranch("GenMT2_WH",       mt2WH)

#        self.out.fillBranch("ngenmt2",          imt2)
#        if imt2 > 0:
#          self.out.fillBranch("GenMT2_bbjj",       mt2bbjj)
#          self.out.fillBranch("GenMT2_bjjb",       mt2bjjb)

        self.out.fillBranch("nLepFromTop", sum( [ l.fromTop for l in leptons ] ) )
        self.out.fillBranch("nLepFromTau", sum( [ l.fromTau for l in leptons ] ) )
        self.out.fillBranch("nLepFromW",   sum( [ l.fromW for l in leptons ] ) )
        self.out.fillBranch("nLepFromZ",   sum( [ l.fromZ for l in leptons ] ) )

        self.out.fillBranch("nGenTau",        len(taus) )

        if len(GenJets)>0:
            self.out.fillBranch("GenJet_fromW", [ j.fromW for j in GenJets])
            self.out.fillBranch("GenJet_genPartWIdx", [ j.genPartWIdx for j in GenJets])
            self.out.fillBranch("GenJet_WIdx", [ j.WIdx for j in GenJets])

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
            

        self.out.fillBranch("ngenW",          len(Ws) )
        if len(Ws)>0:
            self.out.fillBranch("GenW_pt",        [ p.pt    for p in Ws ])
            self.out.fillBranch("GenW_eta",       [ p.eta   for p in Ws ])
            self.out.fillBranch("GenW_phi",       [ p.phi   for p in Ws ])
            self.out.fillBranch("GenW_pdgId",     [ p.pdgId for p in Ws ])
        
        self.out.fillBranch("ngenH",          len(Hs) )
        if len(Hs)>0:
            self.out.fillBranch("GenH_pt",        [ p.pt    for p in Hs ])
            self.out.fillBranch("GenH_eta",       [ p.eta   for p in Hs ])
            self.out.fillBranch("GenH_phi",       [ p.phi   for p in Hs ])
            self.out.fillBranch("GenH_pdgId",     [ p.pdgId for p in Hs ])

        self.out.fillBranch("ngenb",          len(bs) )
        if len(bs)>0:
            self.out.fillBranch("Genb_pt",        [ p.pt    for p in bs ])
            self.out.fillBranch("Genb_eta",       [ p.eta   for p in bs ])
            self.out.fillBranch("Genb_phi",       [ p.phi   for p in bs ])
            self.out.fillBranch("Genb_pdgId",     [ p.pdgId for p in bs ])
            self.out.fillBranch("Genb_fromH",     [ p.fromH for p in bs ])

        self.out.fillBranch("ngenj",          len(js) )
        if len(js)>0:
            self.out.fillBranch("Genj_pt",        [ p.pt    for p in js ])
            self.out.fillBranch("Genj_eta",       [ p.eta   for p in js ])
            self.out.fillBranch("Genj_phi",       [ p.phi   for p in js ])
            self.out.fillBranch("Genj_pdgId",     [ p.pdgId for p in js ])
            self.out.fillBranch("Genj_fromW",     [ p.fromW for p in js ])

        self.out.fillBranch("ngenak8jet",        len(GenAK8Jets) )
        if len(GenAK8Jets) > 0:
            self.out.fillBranch("Genak8j_pt",       [j.pt   for j in GenAK8Jets])
            self.out.fillBranch("Genak8j_eta",        [j.eta   for j in GenAK8Jets])
            self.out.fillBranch("Genak8j_phi",        [j.phi   for j in GenAK8Jets])

            self.out.fillBranch("W_pt",        [ p.pt    for p in Ws ])
            self.out.fillBranch("W_eta",       [ p.eta   for p in Ws ])
            self.out.fillBranch("W_phi",       [ p.phi   for p in Ws ])
            self.out.fillBranch("W_mass",      [ p.mass  for p in Ws ])
            self.out.fillBranch("W_pdgId",     [ p.pdgId for p in Ws ])
            self.out.fillBranch("W_genPartIdx",   [ p.idx for p in Ws ])

        self.out.fillBranch("nTop",          len(tops) )
        if len(tops)>0:
            self.out.fillBranch("Top_pt",        [ p.pt    for p in tops ])
            self.out.fillBranch("Top_eta",       [ p.eta   for p in tops ])
            self.out.fillBranch("Top_phi",       [ p.phi   for p in tops ])
            self.out.fillBranch("Top_pdgId",     [ p.pdgId for p in tops ])

        self.out.fillBranch("stitch",     stitch  )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genAnalyzer = lambda  isW, isWExt, isT, isTExt, year : GenAnalyzer( isW=isW, isWExt=isWExt, isT=isT, isTExt=isTExt, year=year )
