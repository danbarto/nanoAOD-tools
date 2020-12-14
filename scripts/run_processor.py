#!/usr/bin/env python

import sys
import re

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.GenAnalyzer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.lumiWeightProducer import *

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2018, puWeight_2017, puWeight_2016

files       = sys.argv[1].split(',')
lumiWeight  = float(sys.argv[2])
isData      = int(sys.argv[3]) == 1
year        = int(sys.argv[4])
era         = sys.argv[5]
isFastSim   = int(sys.argv[6]) == 1

if files[0].startswith('/store/'):
    print "Had to add a prefix"
    #files = [ 'root://xrootd.t2.ucsd.edu:2040/' + f for f in files ] #FIXME xcache not working all the time..
    files = [ 'root://cmsxrootd.fnal.gov/' + f for f in files ]

#json support to be added

print "Sumweight:", lumiWeight
print "Data:", isData
print "Year:", year
print "FastSim", isFastSim 
print "Files:", files

# corrector for Type-1 MET
met = createJMECorrector(isMC=(not isData), dataYear=year, runPeriod=era, jesUncert="Total", jetType = "AK4PFchs", applySmearing = False, isFastSim = isFastSim )

# corrector for AK8 jets
AK8_type = "AK8PFchs" if int(year) == 2016 and isFastSim else "AK8PFPuppi" # there are no puppi JECs for 2016 rn. #FIXME needs to be checked!
print "Using the following jet type for AK8:", AK8_type
ak8 = createJMECorrector(isMC=(not isData), dataYear=year, runPeriod=era, jesUncert="Total", jetType = AK8_type, applySmearing = True, isFastSim = isFastSim )

modules = [\
    lumiWeightProd(lumiWeight, isData),
    met(),
    ak8(),
    ]

if not isData:
    isW = False
    isWExt = False
    for f in files:
        print f
        if re.search("W[1-4]Jets", f): isW=True
        if re.search("NuPt", f): isWExt=True
        print isW, isWExt
    modules += [genAnalyzer(isW, isWExt)]
    modules += [puWeight_2018()]

modules += [\
    selector(year, isData),
    ]

if isData:
    if year==2018:
        jsonInput='PhysicsTools/NanoAODTools/data/lumi/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
    elif year==2017:
        jsonInput='PhysicsTools/NanoAODTools/data/lumi/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    elif year==2016:
        jsonInput='PhysicsTools/NanoAODTools/data/lumi/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
    else:
        jsonInput=None # FIXME
else:
    jsonInput=None

# apply PV requirement	
cut  = 'PV_ndof>4 && sqrt(PV_x*PV_x+PV_y*PV_y)<=2 && abs(PV_z)<=24'
# loose skim	
#cut += '&& (MET_pt>200 || (nElectron+nMuon)>=2)' # this is a too loose lepton skim
cut += '&& (MET_pt>200 || (Sum$(Electron_pt>10 && abs(Electron_eta)<2.5 && Electron_miniPFRelIso_all<0.2) + Sum$(Muon_pt>10 && abs(Muon_eta)<2.5 && Muon_miniPFRelIso_all<0.2 && Muon_looseId))>=2)'
#cut += '&& MET_pt>200'
cut += '&& (Sum$(Jet_pt>30&&abs(Jet_eta<2.4))>=2 || Sum$(FatJet_pt>170&&abs(FatJet_eta<2.4))>=1)'

#p = PostProcessor('./', files, cut=cut, modules=modules,fwkJobReport=True, prefetch=True,\
p = PostProcessor('./', files, cut=cut, modules=modules,fwkJobReport=True, prefetch=True,\
#    branchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop_in.txt',\
    outputbranchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop.txt',
    jsonInput=jsonInput )

p.run()
