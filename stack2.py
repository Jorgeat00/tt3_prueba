from config import *
from PDFscaleUncertainties import *

def Draw(var, categories, output=None, label='', outpath='temp/', doRatio=True):
  
  plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)

  
  if not CheckHistoCategories(plt.hists[var], categories):
    print("Nope")
    return
  plt.SetRatio(doRatio)
  plt.SetOutpath(outpath)
  plt.plotData = doData
  plt.SetLumi(lumi, "fb$^{-1}$", '13 TeV')
  plt.SetRatioRange(0.8,1.2)
  #plt.SetLogY()
  if var in ['counts', 'l0pt','ept', 'mpt', 'l0eta', 'eeta', 'meta', 'njets','nbtagsl','nbtagsm','ptll']:
  #if var in ['counts','nbtagsm']:
    categories['sign'] = 'OS'
  plt.SetCategories(categories)
  plt.SetDataName('data')
  label = (GetChLab(categories['channel']) if isinstance(categories['channel'], str) else GetChLab(categories['channel'][0]) ) + GetLevLab(categories['level'])
  #AddLabel(self, x, #y, text, options={}):
  plt.SetRegion(label)
  plt.SetOutput(output)
  
  b0 = None; bN = None
  if   var in ['deltaphi','invmass','ptll']:
  #if   var in ['invmass']:
    b0 = 2

  if b0 is not None:
    plt.SetRebin(var, b0, bN, includeLower=True, includeUpper=True)
  
  
  #PrintHisto(plt.hists[var]) #esto printea las variables que tiene el histograma, cada una de las ramas

  ############## Uncertainties
  # Muon and electron efficiencies --> From SF, already included in the histograms, automatic
  # FSR, ISR -> From weights, already included in the histograms, automatic
  # Pileup -> Not yet...

  # JES -> From JES tt sample without JEC (but, for the moment, flat)
  jecs = 0.0#GetJECSystHistos(path, 'variations/TTTo2L2Nu_withJEC', var='counts', categories=categories)[0]

  # hdamp -> Flat on tt
  '''
  categoriesh=categories
  del categoriesh['syst']
  pathvar=path+'/variaciones/'
  hdampup,hdampdo = GetModSystHistos(pathvar, 'TTTo2L2Nu_hdamp', 'hdamp', var='counts')
  #print(hdampup.values())
  #print('--------')
  for c in categoriesh:
    hdampup = hdampup.integrate(c, categoriesh[c])
    hdampdo = hdampdo.integrate(c, categoriesh[c])
  #print(hdampup.values())
  hdampup = hdampup.integrate('syst', 'hdampUp')
  hdampdo = hdampdo.integrate('syst', 'hdampDown')
  '''
  '''
  #print('hdamp', hdampMean)
  # For the moment, hard-coded:
  '''
  ######## correctedPuppiJets // primera version (sin JECS)
  #hdamp = 0.01335580063899193##0.009859073495631691#0.0082 # 0.82%
  hdamp=0.0052
  #mtop=0.005972

  # PDF and scales -> From weights, flat on tt
  #categoriesPDF = {'channel':'em', 'level': 'g2jetsg1b'}
  #pdf   = Get1bPDFUnc(  path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
  #scale = Get1binScaleUnc(path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)

  #totFlatTTunc = np.sqrt(hdamp**2 + pdf**2 + scale**2) #jecs**2
  totFlatTTunc=0.01
  # Background uncertainties
  #plt.SetNormUncDict({'tt':totFlatTTunc, 'tW':0.15,'semilep':0.2, 'WJets':0.3, 'DY':0.2, 'Diboson':0.3})
  plt.SetNormUncDict({'tt':totFlatTTunc, 'tW':0.10,'Non-prompt':0.3, 'DY':0.2, 'VV':0.2, 'ttX': 0.3})
  plt.SetSystematics(syst=['lepSF_muon', 'lepSF_elec','PU','trigSF','ISR' ,'JES','JER','FSR','btagging','PreFiring'])#'ISR', 'FSR',
  #JES y JER hacen que nbtagsm sea asimetrico, 
  #plt.SetSystematics(syst=['lepSF_muon', 'lepSF_elec','PU','trigSF', 'ISR','FSR', 'btagging','PreFiring'])#'ISR', 'FSR',
  #plt.SetSystematics(syst=['stat']) #'norm'
  #plt.SetSystematics(syst=['lepSF_muon', 'lepSF_elec','PU','trigSF', 'ISR','FSR', 'btagging','PreFiring',f'JER_{year}','JES_FlavorQCD', 'JES_Absolute', 'JES_RelativeBal', 'JES_BBEC1', 'JES_RelativeSample'])
  plt.Stack(var, xtit='', ytit='Sucesos', dosyst=True)
  #print(plt.Stack())

  #plt.PrintYields(var='counts', save=True, bkgprocess={'tW','semilep','WJets','DY','Diboson'}, signalprocess={'tt'}) #aqui deberia salir el numero de sucesos pero no sale nada

outpath = outpatho
#outpath = '/nfs/fanae/user/andreatf/www/private/ttrun3/withLepSF_withoutJECPU_metfiltersOverlap_correctLepSF_recoMuonSF_PU_triggerSF_withMllfixed_puppiJetsCorrectedRecommended_newlLepJetCleaning/' 
outpath='/nfs/fanae/user/jorgeat/Documentos/tfm/cafea/prueba/'+year 
#outpath='/nfs/fanae/user/jorgeat/Documentos/tfm/cafea/pruebaSF/trigger/'+year 
'''
def Print2lplots():
  for c in ['em','ee','mm']:#,'ee','mm']:#','em', 'mm']:
    for l in ['dilep','g2jetsg1b']:#,'g2jetsg1b','g2jetsg2b']: #g2jets, g2jetsg2b
      outp = outpath+'/'+l+'/'
      cat = {'channel':c, 'level':l}#,'syst':'norm'}
      for var in ['invmass', 'invmass2', 'j0pt']: #, 'invmass2','met'
        #for var in ['invmass']:
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, output=outname, outpath=outp)
      #for var in ['counts', 'l0pt','ept', 'mpt', 'l0eta', 'eeta', 'meta', 'njets','nbtagsl','nbtagsm','nvtxPU','deltaphi','ptll']:
      for var in ['counts', 'l0pt', 'l0eta', 'njets','nbtagsm']:
      #for var in ['counts', 'l0pt','l0eta', ]:
      #for var in ['counts','nbtagsm']:
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp)


if not var is None:
  var='invmass'
  ch='em'; level='g2jetsg1b'
  categories = { 'channel': ch, 'level' : level}#, 'syst':'norm'}
  outp = outpath+'/'+level+'/'
  outname = "%s_%s_%s"%(var, ch, level)
  Draw(var, categories, outname, outpath=outp)

else:
#Draw('njets', {'channel':['em'], 'level':'dilep','syst':'norm'}, output='njetsmalnrm', outpath=outpath)
  Print2lplots()

'''
var='nbtagsm'
ch='em'; level='g2jets'
categories = { 'channel': ch, 'level' : level}#, 'syst':'norm'}
outp = outpath+'/'+level+'/'
outname = "%s_%s_%s"%(var, ch, level) #_sinJES_JER
Draw(var, categories, outname, outpath=outp)

#h=GetHisto('histosprueba/SingleMuon.pkl.gz') #esto funciona
#PrintHisto(h[var])
#print(h)

