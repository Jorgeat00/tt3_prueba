#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import copy
import numpy as np
import awkward as ak
import coffea
import sys
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

from cafea.modules.GetValuesFromJsons import get_param, get_lumi
from cafea.analysis.objects import *
#from cafea.analysis.corrections import ApplyJetSystematicsRun3, GetPUSF_run3,  AttachMuonSFsRun3, AttachElecSFsRun3, AttachTrigSFsRun3, ApplyJetCorrectionsRun3
from cafea.analysis.corrections import TriggerSF, ElectronSF, MuonSF, GetBtagSF_R2, GetPUSF, ApplyJetSystematicsRun2,ApplyJetCorrectionsRun2, MuonsplitSF, ElectronsplitSF, ApplyRochesterCorrections
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path

#el cambio de norm se ve buscando CAMBIO

doSyst = False
doJES = False
'''
def AttachTrigSF(e0, m0, events):
  TrigSFe, TrigSFedo, TrigSFeup = GetTriggerSF5TeV(e0.pt, np.abs(e0.eta), 'e')
  TrigSFm, TrigSFmdo, TrigSFmup = GetTriggerSF5TeV(m0.pt, np.abs(m0.eta), 'm')
  TrigSFe   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFe, 1.)), nan=1)
  TrigSFedo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFedo, 1.)), nan=1)
  TrigSFeup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFeup, 1.)), nan=1)
  TrigSFm   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFm, 1.)), nan=1)
  TrigSFmdo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmdo, 1.)), nan=1)
  TrigSFmup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmup, 1.)), nan=1)
  events['trigger_sf']    = TrigSFe*TrigSFm
  events['trigger_sfUp'] = TrigSFeup*TrigSFmup
  events['trigger_sfDown'] = TrigSFedo*TrigSFmdo
'''
class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):

        self._samples = samples

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({
        'dummy'      : hist.Hist("Dummy", hist.Cat("sample", "sample"), hist.Bin("dummy", "Number of events", 1, 0, 1)),
        'PDF'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("PDF",     "Counts", 103, 0, 103)),
        'Scales'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Scales",  "Counts", 9, 0, 9)),
        'counts'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_now'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("counts_now",  "Counts", 1, 0, 10)),
        'l0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("l0pt",  "Leading lepton $p_{T}$ (GeV)", 10, 20, 120)),
        'l0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("l0eta", "Leading lepton $\eta$ ", 10, -2.5, 2.50)),
        'l1pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("l1pt",  "Subleading lepton $p_{T}$ (GeV)", 10, 20, 120)),
        'l1eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("l1eta", "Subleading lepton $\eta$ ", 10, -2.5, 2.50)),
        'ept'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("ept",  "Electron $p_{T}$ (GeV)", 10, 20, 120)),
        'eeta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("eeta", "Electron $\eta$ ", 10, -2.5, 2.50)),
        'mpt'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("mpt",  "Muon $p_{T}$ (GeV)", 10, 20, 120)),
        'meta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("meta", "Muon $\eta$ ", 10, -2.5, 2.50)),
        'j0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0pt",  "Leading jet $p_{T}$ (GeV)", 10, 0, 300)),
        'j0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0eta", "Leading jet $\eta$ ", 12, -2.5, 2.50)),
        'deltaphi'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("deltaphi","$\Delta\\varphi (e\mu)$ (rad/$\pi$)", 10, 0, 1)),
        'invmass'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 20, 0, 300)),
        'invmass2'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass2", "$m_{\ell\ell}$ (GeV) ", 30, 75, 105)),
        'invmass3'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass3", "$m_{\ell\ell}$ (GeV) ", 30, 85, 95)),
        'njets'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("njets",   "Jet multiplicity", 6, 0, 6)),
        "nbtagsl"    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("nbtagsl", "Loose btag multiplicity ", 5, 0, 5)),
        "nbtagsm"    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("nbtagsm", "Medium btag multiplicity ", 5, 0, 5)),
        'met'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("met",     "MET (GeV)", 10, 0, 200)),
        'ht'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ht",      "H$_{T}$ (GeV)", 10, 0, 400)),
        'nvtxPU'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("nvtxPU",      "Number of vertex", 40, 0, 80)),
        'ptll'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("ptll",  "$p_{T}^{\ell\ell}$ (GeV)", 26, 40, 300)),
        #'caloPU'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("caloPU",      "rhoCentralcaloPU", 40, 0, 80)),
        #'chargedPU'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("chargedPU",      "rhoCentralchargedPUPileUp", 40, 0, 80)),
        })

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    # Main function: run on a given dataset
    def process(self, events):
        # Dataset parameters
        dataset = events.metadata["dataset"]
        histAxisName = self._samples[dataset]["histAxisName"]
        year         = self._samples[dataset]["year"]
        
        xsec         = self._samples[dataset]["xsec"]
        sow          = self._samples[dataset]["nSumOfWeights"]
        isData       = self._samples[dataset]["isData"]
        isSystSample = ('mtop' in histAxisName) or ('hdamp' in histAxisName) or ('UE' in histAxisName)
        doPS         = (histAxisName in ['tt', 'ttPS', 'TTTo2L2Nu']) and events.PSWeight is not None and len(events.PSWeight[0])>=4
        #doPDFunc= False
        doPDFunc = "sumPDFWeights" in self._samples[dataset]
        #AVISO: para calcular las PDF los json tienen que estar actualizados
        
        #if histAxisName in ['tt', 'ttPS', 'TTTo2L2Nu']: doSyst = True
        #else: doSyst = False
        # Golden JSON !
        if year == "2016" or year == "2016APV":
            golden_json_path = cafea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
        elif year == "2017":
            golden_json_path = cafea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
        elif year == "2018":
            golden_json_path = cafea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = np.ones_like(events['event'], dtype=bool)
        if isData:
          lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)
        if doPDFunc:
          sowPDF       = self._samples[dataset]["sumPDFWeights"]
          sowScale     = self._samples[dataset]["sumScaleWeights"]
          PDFnorm = 1./np.array(sowPDF)
          Scalenorm = 1./np.array(sowScale)
          scaleweights      = events.LHEScaleWeight.to_numpy()
          scaleweights_bins = ak.local_index(events.LHEScaleWeight)
          pdfweights        = events.LHEPdfWeight.to_numpy()
          pdfweights_bins   = ak.local_index(events.LHEPdfWeight)
          scaleweights      = scaleweights * Scalenorm
          pdfweights        = pdfweights * PDFnorm
        dataset = events.metadata["dataset"].split('_')[0] 
        # Initialize objects
        met  = events.MET
        e    = events.Electron
        mu   = events.Muon
        jets = events.Jet

        # Pre-selection 
        #e["btagDeepB"] = ak.fill_none(e.matched_jet.btagDeepB, -99)
        #mu["btagDeepB"] = ak.fill_none(mu.matched_jet.btagDeepB, -99)

        # Muon selection
        #mu["isLoose"] = MuonLoose(mu.pt, mu.eta, mu.dxy, mu.dz, mu.sip3d, mu.mediumId, mu.btagDeepB, ptCut=20, etaCut=2.4)
        #mu["isMVA"]   = MuonMVA(mu.miniPFRelIso_all, mu.mvaTTH)
        mu["pt"] = ApplyRochesterCorrections(year, mu, isData)
        mu['isGood'] = isMuonPOGT(mu, ptCut=20) #etaCut=1.5
        mu['isExtra'] = isMuonPOGT(mu, ptCut=10) #forveto

        # Electron selection
        #e['isLoose'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaNoIso_WPL, 20, 2.4)
        #e['isMVA']   = ElecMVA(e.miniPFRelIso_all, e.mvaTTH)
        #if not hasattr(events, "fixedGridRhoFastjetAll"): events["fixedGridRhoFastjetAll"] = np.zeros_like(events, dtype=float)
        #AttachCutBasedTight(e, events.fixedGridRhoFastjetAll)
        e['isGood'] = isElectronTight(e, ptCut=20, etaCut=2.5) #estaba en 20 y 2.5 #prueba con 25 y 1.5
        e['isExtra'] = isElectronTight(e, ptCut=10, etaCut=2.5) #forveto #1.5
        # Build good collections
        #m_sel = mu[mu.isLoose & mu.isMVA]
        #e_sel = e[e.isLoose & e.isMVA]
        m_sel = mu[mu.isGood]
        e_sel = e[e.isGood]

        e0 = e_sel[ak.argmax(e_sel.pt,axis=-1,keepdims=True)]
        m0 = m_sel[ak.argmax(m_sel.pt,axis=-1,keepdims=True)]

        # Build loose collections
        m_extra = mu[mu.isExtra]
        e_extra = e[e.isExtra]
        #print('Num elecs = ', np.sum(ak.num(e_sel)), '\n\n')

        l_sel = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMCandidate')
        l_sel_extra = ak.with_name(ak.concatenate([e_extra, m_extra], axis=1), 'PtEtaPhiMCandidate')
        llpairs = ak.combinations(l_sel, 2, fields=["l0","l1"])
        mll = (llpairs.l0+llpairs.l1).mass # Invmass for leading two leps
        deltaphi = (llpairs.l0.delta_phi(llpairs.l1))/np.pi

        
        l_sel_padded = ak.pad_none(l_sel, 2)
        l0 = l_sel_padded[:,0]
        l1 = l_sel_padded[:,1]

        leadinglep = l_sel[ak.argmax(l_sel.pt, axis=-1, keepdims=True)]
        subleadinglep = l_sel[ak.argmin(l_sel.pt, axis=-1, keepdims=True)]
        leadingpt = ak.flatten(leadinglep.pt) #ak.pad_none(l_sel.pt, 1)
        subleadingpt = ak.flatten(subleadinglep.pt) #ak.pad_none(l_sel.pt, 1)
        leadingeta=ak.flatten(leadinglep.eta)
        subleadingeta=ak.flatten(subleadinglep.eta)

        ### Attach scale factors
        if not isData:
          MuonSF(m_sel,year)
          ElectronSF(e_sel,year)
          #MuonsplitSF(m_sel,year)
          #ElectronsplitSF(e_sel,year)
        
        l_sel = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMCandidate')
        if not isData:
          #AddSFsRun3(events, l_sel)
          AddSFsR2(events,l_sel)
          #AddSFssplitR2(events,l_sel)
        events['isem'] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 1) & (ak.num(l_sel_extra) <= 2) #ak.num(l_sel_extra) <= 2
        events['ismm'] = (ak.num(m_sel) == 2) & (ak.num(e_sel) == 0) & (ak.num(l_sel_extra) <= 2)
        events['isee'] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 2) & (ak.num(l_sel_extra) <= 2)
        events['isOS'] = (ak.prod(l_sel.charge, axis=1) == -1) & (ak.num(l_sel_extra) <= 2)
        events['isSS'] = (ak.prod(l_sel.charge, axis=1) ==  1) & (ak.num(l_sel_extra) <= 2)
        #GetTriggerSF(2018, events, l0, l1) # from top EFT
        #if not isData: AttachTrigSFsRun3(events, e0, m0)
        if not isData:
          TriggerSF(year, events, leadinglep,subleadinglep)
          PUSF=GetPUSF(events.Pileup.nTrueInt,year); PUSFUp=GetPUSF(events.Pileup.nTrueInt,year,var=1); PUSFDown=GetPUSF(events.Pileup.nTrueInt,year,var=-1)
        hout = self.accumulator.identity()
        
        ## Add systematics
        if (isData): genw = np.ones_like(events["event"])
        else:          genw = events["genWeight"]
        weights_dict = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        #weights_dict.add("norm",genw if isData else (xsec/sow)*genw)
        weights_dict.add("norm",genw if isData else (xsec/sow)*genw) #get_lumi(year.upper())*1000.0 (para el Run 2)
        if not isData: # Apply SFs
          weights_dict.add("trigSF", copy.deepcopy(events.trigger_sf), copy.deepcopy(events.trigger_sfUp), copy.deepcopy(events.trigger_sfDown))
          weights_dict.add("lepSF_muon", copy.deepcopy(events.sf_muon), copy.deepcopy(events.sf_hi_muon), copy.deepcopy(events.sf_lo_muon))
          weights_dict.add("lepSF_elec", copy.deepcopy(events.sf_elec), copy.deepcopy(events.sf_hi_elec), copy.deepcopy(events.sf_lo_elec))
          #estos son para la prueba
          #weights_dict.add("lepSF_ID_muon", copy.deepcopy(events.sf_ID_muon), copy.deepcopy(events.sf_hi_ID_muon), copy.deepcopy(events.sf_lo_ID_muon))
          #weights_dict.add("lepSF_tight_elec", copy.deepcopy(events.sf_tight_elec), copy.deepcopy(events.sf_hi_tight_elec), copy.deepcopy(events.sf_lo_tight_elec))
          #weights_dict.add("lepSF_ISO_muon", copy.deepcopy(events.sf_ISO_muon), copy.deepcopy(events.sf_hi_ISO_muon), copy.deepcopy(events.sf_lo_ISO_muon))
          #weights_dict.add("lepSF_RECO_elec", copy.deepcopy(events.sf_RECO_elec), copy.deepcopy(events.sf_hi_RECO_elec), copy.deepcopy(events.sf_lo_RECO_elec))
          #['lepSF_ID_muon', 'lepSF_tight_elec', 'lepSF_ISO_muon','lepSF_RECO_elec']
          #weights_dict.add("PU", copy.deepcopy(PUSF), copy.deepcopy(PUSFUp), copy.deepcopy(PUSFDown))
          weights_dict.add("PU", events.puWeight, events.puWeightUp,events.puWeightDown)
          weights_dict.add('PreFiring', events.L1PreFiringWeight.Nom,  events.L1PreFiringWeight.Up,  events.L1PreFiringWeight.Dn)
        
         # PS = ISR, FSR (on ttPS only)
        if doPS: 
          i_ISRdown = 0; i_FSRdown = 1; i_ISRup = 2; i_FSRup = 3
          ISRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRup)])
          ISRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRdown)])
          FSRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRup)])
          FSRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRdown)])
          weights_dict.add('ISR', np.ones_like(events["event"]), ISRUp, ISRDo)
          weights_dict.add('FSR', np.ones_like(events["event"]), FSRUp, FSRDo)
        
        '''
        obj_correction_syst_lst = [
            f'JER_{year}Up',f'JER_{year}Down', # Systs that affect the kinematics of objects
            'JES_FlavorQCDUp', 'JES_AbsoluteUp', 'JES_RelativeBalUp', 'JES_BBEC1Up', 'JES_RelativeSampleUp', 'JES_FlavorQCDDown', 'JES_AbsoluteDown', 'JES_RelativeBalDown', 'JES_BBEC1Down', 'JES_RelativeSampleDown']
        #obj_correction_syst_lst = [
            #'JER_{year}Up',f'JER_{year}Down', # Systs that affect the kinematics of objects
            #'JES_FlavorQCDUp', 'JES_RelativeBalUp', 'JES_HFUp', 'JES_BBEC1Up', 'JES_EC2Up', 'JES_AbsoluteUp', f'JES_BBEC1_{year[0:4]}Up', f'JES_EC2_{year[0:4]}Up', f'JES_Absolute_{year[0:4]}Up', f'JES_HF_{year[0:4]}Up', f'JES_RelativeSample_{year[0:4]}Up','JES_FlavorQCDDown', 'JES_RelativeBalDown', 'JES_HFDown', 'JES_BBEC1Down', 'JES_EC2Down', 'JES_AbsoluteDown', f'JES_BBEC1_{year[0:4]}Down', f'JES_EC2_{year[0:4]}Down', f'JES_Absolute_{year[0:4]}Down', f'JES_HF_{year[0:4]}Down', f'JES_RelativeSample_{year[0:4]}Down']
        
        '''
        met_raw=met #esto es para aplicar las correcciones
        #if doJES and not isData: systJEC_list = ['nominal','JERUp','JERDown',"JES_AbsoluteStatDown", "JES_AbsoluteScaleDown", "JES_AbsoluteFlavMapDown", "JES_AbsoluteMPFBiasDown", "JES_FragmentationDown", "JES_SinglePionECALDown", "JES_SinglePionHCALDown", "JES_FlavorQCDDown", "JES_TimePtEtaDown", "JES_RelativeJEREC1Down", "JES_RelativePtBBDown", "JES_RelativePtEC1Down", "JES_RelativeBalDown", "JES_RelativeSampleDown", "JES_RelativeFSRDown", "JES_RelativeStatFSRDown", "JES_RelativeStatECDown", "JES_PileUpDataMCDown", "JES_PileUpPtRefDown", "JES_PileUpPtBBDown", "JES_PileUpPtEC1Down", "JES_AbsoluteStatUp", "JES_AbsoluteScaleUp", "JES_AbsoluteFlavMapUp", "JES_AbsoluteMPFBiasUp", "JES_FragmentationUp", "JES_SinglePionECALUp", "JES_SinglePionHCALUp", "JES_FlavorQCDUp", "JES_TimePtEtaUp", "JES_RelativeJEREC1Up", "JES_RelativePtBBUp", "JES_RelativePtEC1Up", "JES_RelativeBalUp", "JES_RelativeSampleUp", "JES_RelativeFSRUp", "JES_RelativeStatFSRUp", "JES_RelativeStatECUp", "JES_PileUpDataMCUp", "JES_PileUpPtRefUp", "JES_PileUpPtBBUp", "JES_PileUpPtEC1Up"]#'JES_FlavorQCDUp', 'JES_SubTotalPileUpUp', 'JES_SubTotalRelativeUp', 'JES_SubTotalAbsoluteUp','JES_TimePtEtaUp','JES_FlavorQCDDown', 'JES_SubTotalPileUpDown', 'JES_SubTotalRelativeDown', 'JES_SubTotalAbsoluteDown','JES_TimePtEtaDown']
        #if doJES and not isData: systJEC_list=['nominal']+obj_correction_syst_lst
        #if doJES and not isData: systJEC_list =['nominal','JESDown', 'JESUp', 'JERDown','JERUp'] #!!!!!!!!!!!!! #orden 2
        if doJES and not isData: systJEC_list =['nominal','JERUp', 'JERDown', 'JESUp','JESDown'] #!!!!!!!!!!!!! #orden 1
        
        #if doJES and not isData: systJEC_list=['nominal','JERUp', 'JERDown','JES_FlavorQCDUp', 'JES_AbsoluteUp', 'JES_RelativeBalUp', 'JES_BBEC1Up', 'JES_RelativeSampleUp', 'JES_FlavorQCDDown', 'JES_AbsoluteDown', 'JES_RelativeBalDown', 'JES_BBEC1Down', 'JES_RelativeSampleDown']
        else: systJEC_list = ['nominal']
        if isData: systJEC_list = ['nominal']
        for syst_var in systJEC_list:
          #hay que hacer una copia de los pesos para que no se cambie el valor "nominal" con cada iteracion
          weights_dict_syst= copy.deepcopy(weights_dict)
          # Jet cleaning, before any jet selection
          #vetos_tocleanjets = ak.with_name( l_sel, "PtEtaPhiMCandidate")
          #tmp = ak.cartesian([ak.local_index(jets.pt), vetos_tocleanjets.jetIdx], nested=True)
          #cleanedJets = jets[~ak.any(tmp.slot0 == tmp.slot1, axis=-1)] # this line should go before *any selection*, otherwise lep.jetIdx is not aligned with the jet index
          j_isclean=ak.all(jets.metric_table(l_sel) > 0.4, axis=2)
          cleanedJets=jets[j_isclean]
          #print('analisis 261')
          #print(cleanedJets.pt)          # Selecting jets and cleaning them
          #print('------------')
          jetptname = "pt_nom" if hasattr(cleanedJets, "pt_nom") else "pt"
          
          
          # Without JEC
          if doJES:
            # Jet energy corrections
            cleanedJets["pt_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.pt
            cleanedJets["mass_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.mass
            if not isData: cleanedJets["pt_gen"] = ak.values_astype(ak.fill_none(cleanedJets.matched_gen.pt, 0), np.float32)
            cleanedJets["rho"] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, cleanedJets.pt)[0]
            events_cache = events.caches[0]
            if not isData:
              #cleanedJets = ApplyJetCorrectionsRun2(year, corr_type='jets').build(cleanedJets, lazy_cache=events_cache) #topEFT
              #cleanedJets=ApplyJetSystematicsRun2(year,cleanedJets,syst_var)
              #met=ApplyJetCorrectionsRun2(year, corr_type='met').build(met_raw, cleanedJets, lazy_cache=events_cache)
              cleanedJets = ApplyJetCorrectionsRun2('jets',year).build(cleanedJets, lazy_cache=events_cache)
              met = ApplyJetCorrectionsRun2('met',year).build(met_raw,cleanedJets, lazy_cache=events_cache)
              cleanedJets=ApplyJetSystematicsRun2(cleanedJets,syst_var)

          else: 
            cleanedJets["pt"]=(1 - cleanedJets.rawFactor)*cleanedJets.pt
            cleanedJets["E"]=(1 - cleanedJets.rawFactor)*cleanedJets.E
            cleanedJets["mass"]=(1 - cleanedJets.rawFactor)*cleanedJets.mass
          
          nvtxPU = events.PV.npvsGood
          #caloPU = events.Rho.fixedGridRhoFastjetCentralCalo
          #chargedPU = events.Rho.fixedGridRhoFastjetCentralChargedPileUp

          ################################ Jet selection
          jetptcut = 30
          metcut = 40
          print(syst_var)
          print('-------------')
          cleanedJets["isGood"] = isTightJet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, jetPtCut=jetptcut)
          goodJets = cleanedJets[cleanedJets.isGood]

          # Count jets
          njets = ak.num(goodJets)
          ht = ak.sum(goodJets.pt,axis=-1)
          j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]
          #btagwpl = get_param("btag_wp_loose_UL18")
          #isBtagJetsLoose = (goodJets.btagDeepFlavB > btagwpl)
          isBtagJetsLoose = isbjetloose(goodJets.btagDeepFlavB,year)
          isNotBtagJetsLoose = np.invert(isBtagJetsLoose)
          nbtagsl = ak.num(goodJets[isBtagJetsLoose])
          #btagwpm = get_param("btag_wp_medium_UL18")
          isBtagJetsMedium = isbjet(goodJets.btagDeepFlavB,year)
          isNotBtagJetsMedium = np.invert(isBtagJetsMedium)
          nbtagsm = ak.num(goodJets[isBtagJetsMedium])       
          
          trig = trgPassNoOverlap(events,isData,dataset,year)  
          METfilters = PassMETfilters(events,isData)
          if not isData:
            btagSF, btagSFUp, btagSFDo = GetBtagSF_R2(goodJets.pt, goodJets.eta, goodJets.hadronFlavour,year, isbjet(goodJets.btagDeepFlavB,year), True)
            weights_dict_syst.add("btagging", copy.deepcopy(btagSF), copy.deepcopy(btagSFUp), copy.deepcopy(btagSFDo))
            #weights_dict.add("btagging", copy.deepcopy(btagSF), copy.deepcopy(btagSFUp), copy.deepcopy(btagSFDo))
          
          # We need weights for: normalization, lepSF, triggerSF, pileup, btagSF...
          systList = ["norm"] #CAMBIO ESTO POR NOMINAL
          #systList = ["nominal"]
          #systJets = ['JESUp', 'JESDo'] if doJES else []
          #systJets= systJEC_list if doJES else []
          #if not isData and not isSystSample: systList = systList + ["lepSFUp","lepSFDown", "trigSFUp", "trigSFDown", "PUUp", "PUDown"]+systJets
          #if not isData and not isSystSample: systList = systList + ["eleceffUp","eleceffDown", "muoneffUp", "muoneffDown", "trigSFUp", "trigSFDown", "PUUp", "PUDown"]+systJets

          #if not isData and not isSystSample: systList = systList + [ "lepSF_elecUp","lepSF_elecDown","lepSF_muonUp","lepSF_muonDown","PUUp", "PUDown","trigSFUp", "trigSFDown"]# 
          if not isData and not isSystSample: systList = systList + [ "lepSF_elecUp","lepSF_elecDown","lepSF_muonUp","lepSF_muonDown","trigSFUp", "trigSFDown","btaggingUp","btaggingDown", "PUUp", "PUDown", "PreFiringUp", "PreFiringDown"]# 'lepSF_ID_muonUp', 'lepSF_tight_elecUp', 'lepSF_ISO_muonUp','lepSF_RECO_elecUp','lepSF_ID_muonDown', 'lepSF_tight_elecDown', 'lepSF_ISO_muonDown','lepSF_RECO_elecDown'
          #if not isData and not isSystSample: systList = systList + ["trigSFUp", "trigSFDown", "lepSF_elecUp","lepSF_elecDown","lepSF_muonUp","lepSF_muonDown", "PUUp", "PUDown", "PreFiringUp", "PreFiringDown"] #,"btaggingUp","btaggingDown"
          #if not isData and not isSystSample: systList = systList + ["trigSFUp", "trigSFDown","lepSF_muonUp","lepSF_muonDown"]
          if doPS: systList += ['ISRUp', 'ISRDown', 'FSRUp', 'FSRDown']

          if not doSyst: systList = ["norm"] #AQUI TAMBIEN LO CAMBIO
          #if not doSyst: systList = ["nominal"]
          
          if not isData and doSyst and syst_var != 'nominal': systList = [syst_var]
          # Counts
          counts = np.ones_like(events['event'], dtype=float)
          counts_now = np.ones_like(events['event'], dtype=float)
          
          # Initialize the out object, channels and levels

          channels = ['em', 'ee', 'mm'] 
          #levels = ['dilep','g2jets', 'offZ', 'metcut','g2jetsg1b','g2jetsg2b','offZg2jets'] # 'g1jets'
          levels=['dilep','g2jetsg1b', 'g2jetsg1bloose']
          # Add selections...

          #Adding secuancial preselection for debugging
          printevents = False
          if printevents == True:
             np.set_printoptions(threshold=sys.maxsize)
             printarray = np.array([events.event,events.luminosityBlock,events.run,trig,events.isem,events.isee,events.ismm,]) 
             #print(printarray.transpose())
             selections = PackedSelection(dtype='uint64')
             print("counts per selec level")
             selections.add("lumimask", lumi_mask)
             selections.add("trigger", trig)
             selections.add("metfilter", METfilters)
             cutlum = selections.all(*["lumimask"])
             print("lumimask",len(counts[cutlum]))
             cuttrig = selections.all(*["lumimask","trigger"])
             print("trigger",len(counts[cuttrig]))
             cutmetfilter = selections.all(*["lumimask","trigger","metfilter"])
             print("metfilter",len(counts[cutmetfilter]))
             selections.add("OS", ( (events.isOS)))
             selections.add("em", ( (events.isem)))
             cutos = selections.all(*["lumimask","trigger","em","OS"])
             print("em_os",len(counts[cutos]))
             mllvalues = np.where(ak.num(mll)==0, [[0]], mll)
             mllvalues = np.where(ak.num(mllvalues)>1, [[0]], mllvalues)
             mllvalues = ak.flatten(mllvalues, axis=1)
             selections.add("mll", ( (mllvalues>20)))
             selections.add("ptl1l2", ( (leadingpt>35) & (subleadingpt>35)))
             cutpt = selections.all(*["lumimask","trigger","em","OS","ptl1l2"])
             cutmll = selections.all(*["lumimask","trigger","em","OS","ptl1l2","mll"])
             print("pt",len(counts[cutpt]))
             print("mll",len(counts[cutmll]))
          selections = PackedSelection(dtype='uint64')
          selections.add("em", ( (events.isem)&(trig)&(METfilters)))
          selections.add("ee", ( (events.isee)&(trig)&(METfilters)))
          selections.add("mm", ( (events.ismm)&(trig)&(METfilters)))
          selections.add("OS", ( (events.isOS)))
          selections.add("SS", ( (events.isSS)))
          selections.add("dilep",  (njets >= 0)&(leadingpt>25)&(subleadingpt>20)&(lumi_mask))#25, 20
          #selections.add('eta', ((abs(leadingeta)<1.7) & (abs(subleadingeta)<1.7)))
          selections.add("g1jets", (njets >= 1))
          selections.add("g2jets", (njets >= 2))
          selections.add("g2jetsg1b", (njets >= 2)&(nbtagsm>=1))
          selections.add("g2jetsg1bloose", (njets >= 2)&(nbtagsl>=1))
          selections.add("g2jetsg2b", (njets >= 2)&(nbtagsm>=2))
          selections.add("0jet", (njets == 0))
          selections.add("1jet", (njets == 1))
          selections.add("2jet", (njets == 2))
          selections.add("3jet", (njets == 3)) 
          selections.add("g4jet", (njets >= 4))
          mllvalues = np.where(ak.num(mll)==0, [[0]], mll)
          mllvalues = np.where(ak.num(mllvalues)>1, [[0]], mllvalues)
          mllvalues = ak.flatten(mllvalues, axis=1)
          selections.add("offZ",   ( np.abs(mllvalues-90) > 15)&(njets >= 2))
          #selections.add('met', (met.pt >=metcut))
          selections.add("metcut", (met.pt >= metcut)&( np.abs(mllvalues-90) > 15)&(njets >= 2))
          selections.add("mll", ( (mllvalues>20)))
          #printarray = np.array(events.event[cut])
          #print(printarray)
          selections.add("offZg2jets", ( np.abs(mllvalues-90) > 15)&(njets >= 2))
          #print(ak.num(l_sel[selections.all(*['em'])]))
          print(njets)
          print('--------------')
          ##### Loop over the hists we want to fill
          #for syst in systList:
          syst = "norm"
          for syst in systList:
           njets_var = njets
           ht_var = ht

           for ch in channels:
            '''
            if syst == "norm":
              for lev in ['0jet', '1jet', '2jet', '3jet', 'g4jet','g2jetsg1b','g2jetsg2b']:
                cuts = [ch] + [lev] + ['mll', 'dilep'] + ['OS']
                cut   = selections.all(*cuts)
                weights = weights_dict.weight(None)
                weight = weights[cut]
                mll_flat = mllvalues[cut]
                #hout['invmass'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_flat, syst=syst, weight=weight)
                #hout['invmass2'].fill(sample=histAxisName, channel=ch, level=lev, invmass2=mll_flat, syst=syst, weight=weight)
                #hout['invmass3'].fill(sample=histAxisName, channel=ch, level=lev, invmass3=mll_flat, syst=syst, weight=weight)
            '''
            for lev in levels:
              cuts = [ch] + [lev] + ['mll', 'dilep']
              #nocut=[ch] + [lev]+['OS']  
              cutsOS = cuts + ['OS']
              cutsSS = cuts + ['SS']
              cut   = selections.all(*cutsOS)
              cutSS = selections.all(*cutsSS)
              #nocutsel=selections.all(*nocut)
              if (syst in ["norm","nominal"]) or (syst in systJEC_list): #MAS CAMBIO
              #if (syst in ["nominal"]) or (syst in systJEC_list):
                weights = weights_dict_syst.weight(None)
                #weights = weights_dict.weight(None)
              else:
                weights = weights_dict_syst.weight(syst)
                #weights = weights_dict.weight(syst)

              # Fill all the variables
              weightSS = weights[cutSS]
              weight = weights[cut]
              #weight = np.ones_like(weights[cut])
              mll_flat = mllvalues[cut]
              ptll = ak.flatten(llpairs.l0.pt[cut]+llpairs.l1.pt[cut])
              #deltaphi_cut = deltaphi[cut]
              lep0pt = ak.flatten(llpairs.l0.pt[cut])
              lep0eta = ak.flatten(llpairs.l0.eta[cut])
              jet0pt  = ak.flatten(j0.pt[cut])
              '''
              if syst_var=='nominal': print('analisis linea 450 nominal')
              if syst_var=='nominal': print(weight)
              if syst_var=='JESUp': print('analisis linea 452 JESUp')
              if syst_var=='JESUp': print(weight)
              if syst_var=='JESDown': print('analisis linea 454 JESDown')
              if syst_var=='JESDown': print(weight)
              
              if syst_var=='JERUp': print('analisis linea 471 JERUp')
              if syst_var=='JERUp': print(weight)
              if syst_var=='JERDown': print('analisis linea 473 JERDown')
              if syst_var=='JERDown': print(weight)
              print('----------------')
              '''
              print(l0.eta[cut])
              print(l1.eta[cut])
              print('------')
              jet0eta = ak.flatten(j0.eta[cut])
              hout['counts_now'].fill(sample=histAxisName, channel=ch, level=lev, counts_now=counts_now[cut],  syst=syst, sign='OS', weight=np.ones_like(weight))
              hout['counts_now'].fill(sample=histAxisName, channel=ch, level=lev, counts_now=counts_now[cut],  syst=syst, sign='SS', weight=np.ones_like(weight))
              hout['counts'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut],  syst=syst, sign='OS', weight=weight)
              hout['counts'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cutSS], syst=syst, sign='SS', weight=weightSS)
              hout['njets'].fill(sample=histAxisName, channel=ch, level=lev, njets=njets[cut], syst=syst, sign='OS', weight=weight)
              hout['njets'].fill(sample=histAxisName, channel=ch, level=lev, njets=njets[cutSS], syst=syst, sign='SS', weight=weightSS)
              hout['nbtagsl'].fill(sample=histAxisName, channel=ch, level=lev, nbtagsl=nbtagsl[cut], syst=syst, sign='OS', weight=weight)
              hout['nbtagsl'].fill(sample=histAxisName, channel=ch, level=lev, nbtagsl=nbtagsl[cutSS], syst=syst, sign='SS', weight=weightSS)
              hout['nbtagsm'].fill(sample=histAxisName, channel=ch, level=lev, nbtagsm=nbtagsm[cut], syst=syst, sign='OS', weight=weight)
              hout['nbtagsm'].fill(sample=histAxisName, channel=ch, level=lev, nbtagsm=nbtagsm[cutSS], syst=syst, sign='SS', weight=weightSS)
              hout['ht'].fill(sample=histAxisName, channel=ch, level=lev, ht=ht[cut], syst=syst, weight=weight)
              hout['deltaphi'].fill(sample=histAxisName, channel=ch, level=lev, deltaphi=ak.flatten(deltaphi[cut]), syst=syst, weight=weight)
              hout['met'].fill(sample=histAxisName, channel=ch, level=lev, met=met.pt[cut], syst=syst, weight=weight)
              hout['l0pt'] .fill(sample=histAxisName, channel=ch, level=lev, l0pt=ak.flatten(llpairs.l0.pt[cut]), syst=syst, sign='OS', weight=weight)
              hout['l0eta'].fill(sample=histAxisName, channel=ch, level=lev, l0eta=ak.flatten(llpairs.l0.eta[cut]), syst=syst, sign='OS', weight=weight)
              hout['l1pt'] .fill(sample=histAxisName, channel=ch, level=lev, l1pt=ak.flatten(llpairs.l1.pt[cut]), syst=syst, sign='OS', weight=weight)
              hout['l1eta'].fill(sample=histAxisName, channel=ch, level=lev, l1eta=ak.flatten(llpairs.l1.pt[cut]), syst=syst, sign='OS', weight=weight)
              hout['ptll'].fill(sample=histAxisName, channel=ch, level=lev, ptll=ak.flatten(llpairs.l0.pt[cut]+llpairs.l1.pt[cut]), syst=syst, sign='OS', weight=weight)
              hout['l0pt'] .fill(sample=histAxisName, channel=ch, level=lev, l0pt=ak.flatten(llpairs.l0.pt[cutSS]), syst=syst, sign='SS', weight=weightSS)
              hout['l0eta'].fill(sample=histAxisName, channel=ch, level=lev, l0eta=ak.flatten(llpairs.l0.eta[cutSS]), syst=syst, sign='SS', weight=weightSS)
              hout['l1pt'].fill(sample=histAxisName, channel=ch, level=lev, l1pt=ak.flatten(llpairs.l1.pt[cutSS]), syst=syst, sign='SS', weight=weightSS)
              hout['ptll'].fill(sample=histAxisName, channel=ch, level=lev, ptll=ak.flatten(llpairs.l0.pt[cutSS]+llpairs.l1.pt[cutSS]), syst=syst, sign='SS', weight=weightSS)
              hout['l1eta'].fill(sample=histAxisName, channel=ch, level=lev, l1eta=ak.flatten(llpairs.l1.pt[cutSS]), syst=syst, sign='SS', weight=weightSS)
              hout['invmass'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_flat, syst=syst, weight=weight)
              hout['invmass2'].fill(sample=histAxisName, channel=ch, level=lev, invmass2=mll_flat, syst=syst, weight=weight)
              hout['invmass3'].fill(sample=histAxisName, channel=ch, level=lev, invmass3=mll_flat, syst=syst, weight=weight)
              hout['nvtxPU'].fill(sample=histAxisName, channel=ch, level=lev, nvtxPU=nvtxPU[cut], syst=syst, weight=weight)
              #hout['caloPU'].fill(sample=histAxisName, channel=ch, level=lev, caloPU=caloPU[cut], syst=syst, weight=weight)
              #hout['chargedPU'].fill(sample=histAxisName, channel=ch, level=lev, chargedPU=chargedPU[cut], syst=syst, weight=weight)

              if lev != 'dilep':
                #jet1pt  = ak.flatten(j1.pt)
                #jet1eta = ak.flatten(j1.eta)
                hout['j0pt'].fill(sample=histAxisName, channel=ch, level=lev, j0pt=jet0pt, syst=syst, weight=weight)
                hout['j0eta'].fill(sample=histAxisName, channel=ch, level=lev, j0eta=jet0eta, syst=syst, weight=weight)
              if ch == 'em':
                e = e_sel; m = m_sel
                ept  = ak.flatten(e.pt [cut])
                eeta = ak.flatten(e.eta[cut])
                mpt  = ak.flatten(m.pt [cut])
                meta = ak.flatten(m.eta[cut])
                hout['ept' ].fill(sample=histAxisName, channel=ch, level=lev, ept =ak.flatten(e.pt [cut  ]), sign='OS', syst=syst, weight=weight)
                hout['eeta'].fill(sample=histAxisName, channel=ch, level=lev, eeta=ak.flatten(e.eta[cut  ]), sign='OS', syst=syst, weight=weight)
                hout['mpt' ].fill(sample=histAxisName, channel=ch, level=lev, mpt =ak.flatten(m.pt [cut  ]), sign='OS', syst=syst, weight=weight)
                hout['meta'].fill(sample=histAxisName, channel=ch, level=lev, meta=ak.flatten(m.eta[cut  ]), sign='OS', syst=syst, weight=weight)
                hout['ept' ].fill(sample=histAxisName, channel=ch, level=lev, ept =ak.flatten(e.pt [cutSS]), sign='SS', syst=syst, weight=weightSS)
                hout['eeta'].fill(sample=histAxisName, channel=ch, level=lev, eeta=ak.flatten(e.eta[cutSS]), sign='SS', syst=syst, weight=weightSS)
                hout['mpt' ].fill(sample=histAxisName, channel=ch, level=lev, mpt =ak.flatten(m.pt [cutSS]), sign='SS', syst=syst, weight=weightSS)
                hout['meta'].fill(sample=histAxisName, channel=ch, level=lev, meta=ak.flatten(m.eta[cutSS]), sign='SS', syst=syst, weight=weightSS)

              # Fill scale and pdf uncertainties
              
              if doPDFunc:
                scale_w = np.transpose(scaleweights[cut])*(weight)
                pdf_w   = np.transpose(pdfweights  [cut])*(weight)
                hout['Scales'].fill(sample=histAxisName, channel=ch, level=lev, Scales=ak.flatten(scaleweights_bins[cut]), syst="norm", weight=ak.flatten(scale_w))
                hout['PDF']   .fill(sample=histAxisName, channel=ch, level=lev, PDF   =ak.flatten(pdfweights_bins[cut]),   syst="norm", weight=ak.flatten(pdf_w))
              
        return hout
  
    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    # Load the .coffea files
    outpath= './coffeaFiles/'
    samples     = load(outpath+'samples.coffea')
    topprocessor = AnalysisProcessor(samples)

