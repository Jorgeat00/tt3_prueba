from config import *
from PDFscaleUncertainties import *

from cafea.modules.CreateDatacardFromRootfile import Datacard
from cafea.plotter.plotter import GetHisto

from rootval_cafea import vartotdo, vartotup,vartot #para las incertidumbres 'shape'
if  not '/' in inputFile: outpath = './'
else: outpath = inputFile[:inputFile.rfind('/')]
if not path.endswith('/'): path += '/'


#channels = ['e', 'm']
#levels = ['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
channels=['em','mm','ee']
levels=['dilep','g2jetsg1b','g2jetsg2b']
verbose=1
def GetChanLevFromName(fname):
  # medianDRjj_e_g5j1b.root
  inputs = fname[fname.rfind('/')+1:].replace('.root', '').split('_')
  chan = None; lev = None
  for i in inputs:
    if i in channels: chan = i
    if i in levels: lev = i
  if chan is None or lev is None:
    print("WARNING: could not get channel or level from file name: %s"%fname)
  return chan, lev

def GetModUnc(path, chan, lev):
  nbin = channels.index(chan) *len(levels) + levels.index(lev)
  if os.path.isfile(path + '/masterhistos/master.pkl.gz'):
    # print('Loading masterhistos from %s'%path)
    '''
    pathvar=path+'/variaciones/'
    ### Add hdamp and tune uncertainties
    hdampup,hdampdo = GetModSystHistos(pathvar, 'TTTo2L2Nu_hdamp', 'hdamp', var='counts')
    #mtopup,mtopdo = GetModSystHistos(pathvar, 'TTTo2L2Nu_mtop', 'mtop', var='counts') #no hay que incluir esta al parecer
    tuneup , tunedo = GetModSystHistos(pathvar, 'TTTo2L2Nu_TuneCP5', 'UE', var='counts')
    '''
    histo = GetHisto(path + 'masterhistos/master.pkl.gz', 'master').integrate('process', 'tt')
    nominal = histo.integrate('syst', 'norm').values()[()][nbin]
    # PDF and scale
    pdfUp = histo.integrate('syst', 'PDFUp').values()[()][nbin]
    pdfDown = histo.integrate('syst', 'PDFDown').values()[()][nbin]
    pdf = (abs(pdfUp-nominal) + abs(pdfDown-nominal))/2/nominal
    scaleUp = histo.integrate('syst', 'ScaleUp').values()[()][nbin]
    scaleDown = histo.integrate('syst', 'ScaleDown').values()[()][nbin]
    scales = (abs(scaleUp-nominal) + abs(scaleDown-nominal))/2/nominal
    # hdamp and UE
    tot = sum(histo.integrate('syst', 'norm').values()[()])
    hdampUp = sum(histo.integrate('syst', 'hdampUp').values()[()])
    hdampDown = sum(histo.integrate('syst', 'hdampDown').values()[()])
    hdamp = max(abs(hdampUp-tot), abs(hdampDown-tot))/tot
    UEUp = sum(histo.integrate('syst', 'UEUp').values()[()])
    UEDown = sum(histo.integrate('syst', 'UEDown').values()[()])
    UE = max(abs(UEUp-tot),abs(UEDown-tot))/tot
    return round(pdf,3), round(scales,3), round(hdamp,3), round(UE,3)
  else:
    print("WARNING: please provide master histograms to take modeling uncertaintes... for now, returning hardcoded values")
    pdf = 0.007
    scales = 0.002
    hdamp = 0.007
    UE = 0.005
  return pdf, scales, hdamp, UE

def GetIncshape(maxs):
	histogramas=['DY','VV','Non-prompt', 'tW', 'ttX']
	var=['lepSF_muon','lepSF_elec','trigSF','btagging','JES','JER','PU','PreFiring']
	DY=[]
	VV=[]
	Nonprompt=[]
	tW=[]
	ttX=[]
	tt=[maxs[-10:]]
	total=[DY,VV,Nonprompt,tW,ttX]
	for i in range(len(histogramas)):
		a0=len(var)*i
		a=len(var)*(i+1)
		total[i].append(maxs[a0:a])
	
	return tt, DY, VV, Nonprompt, tW, ttX
	

def GetIncs(proc,var):
	#var=['lepSF_muon','lepSF_elec','trigSF','btagging','JES','JER','PU','PreFiring']
	muestras=GetIncshape(var)
	#print(muestras)
	nombres={'DY':1, 'VV':2, 'Non-prompt': 3, 'tW': 4, 'ttX': 5, 'tt': 0}
	muestra=muestras[nombres[proc]][0]
	
	lepSF_muon=round(float(muestra[0]),3)
	lepSF_elec=round(float(muestra[1]),3)
	trigSF=round(float(muestra[2]),3)
	btagging=round(float(muestra[3]),3)
	JES=round(float(muestra[4]),3)
	JER=round(float(muestra[5]),3)
	PU=round(float(muestra[6]),3)
	PreFiring=round(float(muestra[7]),3)
	if len(muestra) >8:
		ISR=round(float(muestra[8]),3)
		FSR=round(float(muestra[9]),3)
		return lepSF_muon, lepSF_elec, trigSF, btagging, JES, JER, PU, PreFiring, ISR, FSR
	else:
		return lepSF_muon, lepSF_elec, trigSF, btagging, JES, JER, PU, PreFiring

def CreateDatacard(fname, outpath=outpath, oname=output):
  chan, lev = GetChanLevFromName(fname)
  if oname is None:
    oname = fname[fname.rfind('/')+1:] if '/' in fname else fname
    if oname.endswith('.root'): oname = oname[:-5]
    if '/' in oname: oname[oname.rfind('/')+1:]
  oname = 'dat_'+oname
  if not oname.endswith('.txt'): oname += '.txt'
  
  lumiUnc = 0.025
  #bkg =  ['tW', 'WJets', 'QCD', 'DY']
  #norm = [0.095, 0.2, 0.3, 0.2]
  #bkg=['tW','Non-prompt', 'DY', 'VV', 'ttX']
  bkg=['DY','VV','Non-prompt', 'tW', 'ttX']
  #norm=[0,0,0,0,0]
  #norm=[0.10,0.3, 0.2,0.3,0.3, 0.3]
  norm= [0.2, 0.3, 0.3, 0.1, 0.3]
  signal = 'tt'
  #systList = ['muonSF', 'elecSF', 'btagSF', 'FSR', 'ISR', 'JES', 'prefire']# 'hdamp', 'UE', 'trigSF', 'Scales', 'PDF', 'Prefire']
  systList=['lepSF_muon','lepSF_elec','trigSF','btagging','JES','JER','PU','PreFiring']#,'ISR','FSR']
  systList2=['ISR','FSR']
  #systList=['JES','JER']
  #systList=['lepSF_muon','lepSF_elec']
  #systList=['JES','JER']
  sList=[]
  #systList += ['MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res']
  verbose=1
  d = Datacard(fname, signal, bkg, lumiUnc, norm, sList, nSpaces=12, verbose=verbose) #sList--> systList
  #d.AddExtraUnc('prefiring', 0.014, ['tt', 'tW', 'WJets', 'DY'])
  
  #pdf   = Get1bPDFUnc(  fname, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
  #scale = Get1binScaleUnc(fname, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
  
  pdf, scales, hdamp, UE = GetModUnc(path, chan, lev)
  valpdf=str(1+pdf)
  if pdf <=0.0001: valpdf='-'
  valscales=str(1+scales)
  if scales <=0.0001: valscales='-'
  valhdamp=str(1+hdamp)
  if hdamp <=0.0001: valhdamp='-'
  valUE=str(1+UE)
  if UE <=0.0001: valUE='-'
  
  d.AddExtraUncTH('PDF', valpdf, signal)
  d.AddExtraUncTH('Scales', valscales, signal)
  #if chan == 'e': d.AddExtraUnc('trigElecSF', 0.02, signal)
  #else          : d.AddExtraUnc('trigMuonSF', 0.01, signal)
  d.AddExtraUncTH('hdamp', valhdamp, signal)
  d.AddExtraUncTH('UE', valUE, signal)
  
 
  '''
  for j in range(len(systList)):
    #print(systList[j], i, GetIncs(i)[j])
    vtot=[]
    for i in ['tt']+bkg:
      valup=1+GetIncs(i,vartotup)[j]
      valdo=1-GetIncs(i,vartotdo)[j]
      if valup<0:
        valup=1.001
        #valdo=0.98
      if valdo<0:
        valdo=0.999
        #valup=1.02
      if valup>=1.1 and valdo <=0.9:
        valup=1.1
        valdo=0.9
      val=str(valdo)+'/'+str(valup)

      if valdo >=0.999 and valup <=1.001:
        val='-'
      #val=str(valdo)+'/'+str(valup)
      vtot.append(val)
    d.AddExtraUnc(systList[j],vtot , ['tt']+bkg)
  '''
  valup=1+GetIncs('tt',vartotup)[8]
  valdo=1-GetIncs('tt',vartotdo)[8]
  valISR=str(valdo)+'/'+str(valup)
  if valdo >=0.9999 and valup <=1.0001:
    valISR='-'
  valup=1+GetIncs('tt',vartotup)[9]
  valdo=1-GetIncs('tt',vartotdo)[9]
  valFSR=str(valdo)+'/'+str(valup)
  if valdo >=0.9999 and valup <=1.0001:
    valFSR='-'
  d.AddExtraUncTH('FSR', valFSR,signal)
  d.AddExtraUncTH('ISR', valISR,signal)
  #d.AddExtraUnc('Prefire', 0.01, signal)
  
  #['lepSF_muon','lepSF_elec','trigSF','btagging','JES','JER','PU','PreFiring']
  #      0            1          2         3        4     5     6       7 
  for j in range(len(systList)):
    #print(systList[j], i, GetIncs(i)[j])
    vtot=[]
    if j==0 or j==1 or j==2 or j==4 or j==5:
      for i in ['tt']+bkg:
        valup=1+GetIncs(i,vartotup)[j]
        valdo=1-GetIncs(i,vartotdo)[j]
        if valup<0:
          valup=1.001
          #valdo=0.98
        if valdo<0:
          valdo=0.999
          #valup=1.02
        if valup>=1.1 and valdo <=0.9:
          valup=1.1
          valdo=0.9
        val=str(valdo)+'/'+str(valup)

        if valdo >=0.999 and valup <=1.001:
          val='-'
        #val=str(valdo)+'/'+str(valup)
        vtot.append(val)
    else: 
      for i in ['tt']+bkg:
        vtot.append('-')
    print(vtot)
    d.AddExtraUnc(systList[j],vtot , ['tt']+bkg)
  
  d.SetOutPath(outpath)
  print('outpath = %s'%outpath)
  d.Save(oname)

if inputFile == '': 
  print('Please provide a root file to create datacards from using --inputFile /path/to/inputs/');
  exit(1)

if os.path.isdir(inputFile):
  for d in os.listdir(inputFile):
    if not d.endswith('.root'): continue
    fname = os.path.join(inputFile, d)
    CreateDatacard(fname)
else:
  CreateDatacard(inputFile)
