from config import *
from cafea.plotter.xsec import *
from cafea.modules.fileReader import *
from PDFscaleUncertainties import *

names = {
  'lepSF_muon' : 'Muon efficiences',
  'lepSF_elec' : 'Electron efficiences',
  'btagging': 'b--tagging efficiency',
  #'eleceff' : 'Electron efficiences',
  #'muoneff' : 'Muon efficiences',
  'trigSF' : 'Trigger efficiencies',
  'JES' : 'Jet Energy Scale',
  'JER' : 'Jet Energy Resolution',
  'UE' : 'Underlying event',
  'hdamp' : 'ME/PS matching ($h_\mathrm{damp}$)',
  'mtop' : 'Top mass',
  'ISR' : 'Initial-state radiation',
  'FSR' : 'Final-state radiation',
  'DY' : 'Drell--Yan',
  'PU' : 'Pileup reweighting',
  'semilep' : r'$t\bar{t} \rightarrow 1 \ell$',
  #'JES_FlavorQCD': 'JES FlavorQCD',
  #'JES_SubTotalPileUp': 'JES Pileup',
  #'JES_SubTotalRelative':'JES Relative',
  #'JES_SubTotalAbsolute': 'JES Absolute',
  #'JES_TimePtEta': 'JES Time, pt and eta'
}

#path = 'histos/run3/5jun2022/'
### Fix categories
categories = {'channel':'em', 'level': 'g2jetsg1b', 'sign':'OS'} #al poner el syst se contraen todas las incertidumbres y no se pueden calcular
categoriesPDF = {'channel':'em', 'level': 'g2jets'}#, 'syst':'norm'} 


processDic = {
  'tt': 'TTTo2L2Nu',# ',
  #'tW': 'tbarW, tW',
  'tW' : 'tW_noFullHad, tbarW_noFullHad',
  'semilep':'TTToSemiLeptonic',
  'WJets':'WJetsToLNu_MLM', #WJetsToLNu_MLM
  #'DY': 'DYJetsToLL_M50, DYJetsToLL_M10to50', 
  'DY' : 'DYJetsToLL_M_10to50, DYJetsToLL_M_50_MLM',
  #'Diboson' : 'WW, WZ, ZZ',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  'VV' : 'WWTo2L2Nu, WZ,ZZ, ZZZ_ext1, WWW_ext1, WWZ_ext1, WZZ_ext',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  'ttX' : 'TTZToLL_M_1to10, TTZToLLNuNu_M_10,TTWJetsToLNu,TTWJetsToQQ, TTG',
  #'data' : 'MuonEG, DoubleMuon, SingleMuon,EGamma',
  'data': 'MuonEG, DoubleMuon, SingleMuon, DoubleEG, SingleElectron, EGamma',
}

bkglist    = ['tt', 'tW','semilep', 'WJets', 'DY', 'VV', 'ttX']
### Create plotter
p = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var='counts')


pathvar=path+'/variaciones/'
### Add hdamp and tune uncertainties
hdampup,hdampdo = GetModSystHistos(pathvar, 'TTTo2L2Nu_hdamp', 'hdamp', var='counts')
#mtopup,mtopdo = GetModSystHistos(pathvar, 'TTTo2L2Nu_mtop', 'mtop', var='counts') #no hay que incluir esta al parecer
tuneup , tunedo = GetModSystHistos(pathvar, 'TTTo2L2Nu_TuneCP5', 'UE', var='counts')

p.AddExtraBkgHist([hdampup, hdampdo], add=True)
#p.AddExtraBkgHist([mtopup, mtopdo], add=True)
p.AddExtraBkgHist([tuneup, tunedo], add=True)


### Create xsec object
experimental = ['lepSF_muon', 'lepSF_elec','trigSF','btagging','JES','JER', 'PU', 'PreFiring']# 'PU','PreFiring',,'JER','JES_FlavorQCD','JES_SubTotalPileUp','JES_SubTotalRelative','JES_SubTotalAbsolute','JES_TimePtEta']
modeling = ['hdamp','UE','ISR','FSR'] # ['UE', 'hdamp', 'ISR', 'FSR']'mtop'
#modeling=[]
x = xsec('tt', 0.025, {'tW':0.10,'semilep':0.2,'WJets':0.3, 'DY':0.2, 'VV':0.3, 'ttX': 0.3}, plotter=p, verbose=4, thxsec=833.9, experimental=experimental, modeling=modeling, categories=categories)
#x = xsec('tt', 0.06, {'tW':0.15,'semilep':0.2,'WJets':0.3, 'DY':0.2, 'Diboson':0.3}, plotter=p, verbose=4, thxsec=833.9, categories=categories)
x.SetNames(names)
x.ComputeCrossSection()
#x1.ComputeCrossSection()
#x.ComputeXsecUncertainties()


pdf   = Get1bPDFUnc(  path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=True)
scale = Get1binScaleUnc(path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=True)
x.AddModUnc('PDF$+\\alpha_{S}$', pdf, isRelative=True)
x.AddModUnc('$\mu_R, \mu_F$ scales', scale, isRelative=True)
#x.AddModUnc('mtop', 0.005972, isRelative=True)

#jecs = GetJECSystHistos(path, 'variations/TTTo2L2Nu_withJEC', var='counts', categories=categories)

#x.AddExpUnc('Jet energy scale (external)', jecs, isRelative=True)
#x.AddExpUnc('Jet energy scale', 0.02, isRelative=True)
#x.ComputeXsecUncertainties()

x.ComputeXsecUncertainties()
'''
path_gen = '/nfs/fanae/user/juanr/cafea/histos/TTTo2L2Nu_gen.pkl.gz'
hcounts = GetHisto(path_gen, "counts", categories={"sample":"TTTo2L2Nu", "channel":"all", "level":"all", "syst":"norm", "sign":"all"})
genEvents= list(hcounts.values().values())[0][0]
hfidu = GetHisto(path_gen, "counts", categories={"sample":"TTTo2L2Nu", "channel":"em", "level":"g2jets", "syst":"norm", "sign":"OS"})
fiduEvents = list(hfidu.values().values())[0][0]
BRem = 0.031938
BRSF = 0.015969
BRdilep = 0.108*0.108*9

x.SetGenEvents(genEvents)
x.SetFiduEvents(fiduEvents)
x.SetBR(BRem)
x.SetBRsample(BRdilep)

acc, accunc = x.GetAcceptance()
eff, effunc = x.GetEfficiency()
x.ComputeFiducial()
'''
x.GetYieldsTable()
x.GetUncTable(form='%1.2f')

