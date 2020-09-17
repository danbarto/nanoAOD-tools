#!/usr/bin/env python

import sys

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.GenAnalyzer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.lumiWeightProducer import *

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import *

jetmet = createJMECorrector(isMC=True, dataYear=2018, jesUncert="Total", jetType = "AK4PFchs", applySmearing = True, isFastSim = False )

files = sys.argv[1].split(',')

if files[0].startswith('/store/'):
    print "Had to add a prefix"
    files = [ 'root://cmsxrootd.fnal.gov/' + f for f in files ]

#json support to be added
print "Sumweight:", sys.argv[2]
print "Files:", files

modules = [\
    lumiWeightProd(float(sys.argv[2])),
    jetmet(),
    genAnalyzer(),
    selector2018(),
    ]


# apply PV requirement
cut  = 'PV_ndof>4 && sqrt(PV_x*PV_x+PV_y*PV_y)<=2 && abs(PV_z)<=24'
# loose skim
cut += '&& (Sum$(Electron_pt>30&&abs(Electron_eta)<2.4&&Electron_miniPFRelIso_all<0.1&&Electron_cutBased>=3)+Sum$(Muon_pt>25&&abs(Muon_eta)<2.4&&Muon_mediumId>0&&Muon_miniPFRelIso_all<0.1))>0'
cut += '&& ( (Sum$(Jet_pt>25&&abs(Jet_eta)<2.4)>=4) || (Sum$(Jet_pt>25&&abs(Jet_eta)<2.4)>=2 && (Sum$(Electron_pt>10&&abs(Electron_eta)<2.4)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.4&&Muon_mediumId>0))>=3) )'


p = PostProcessor('./', files, cut=cut, modules=modules,fwkJobReport=True, prefetch=True,\
#    branchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop_in.txt',\
    outputbranchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop.txt' )

p.run()