from config import *

prdic = {'Prompt': 'TTTo2L2Nu, tbarW_noFullHad, tW_noFullHad, WWTo2L2Nu, WZ,ZZ, ZZZ_ext1, WWW_ext1, WWZ_ext1, WZZ_ext1, DYJetsToLL_M_10to50_MLM, DYJetsToLL_M_50_MLMTTZToQQ,TTZToLL_M_1to10, TTZToLLNuNu_M_10,TTWJetsToLNu,TTWJetsToQQ, TTG', 'Nonprompt': 'WJetsToLNu_MLM, TTToSemiLeptonic','data': 'MuonEG,DoubleMuon,SingleMuon,DoubleEG, EGamma'}

p = plotter(path, prDic=prdic, bkgList=bkglist, colors=colordic, lumi=lumi)
p.SetDataName('data')
outpath='/nfs/fanae/user/jorgeat/Documentos/tfm/cafea/prueba/'+year
level='g2jetsg1b'

def GetYieldsNP(plotter, sign='OS', process='Nonprompt', level='g2jetsg1b'):
  cat = {'level':level, 'sign':sign, 'channel':'em', 'syst':'norm'}
  if process == 'data':
    h = plotter.GetHistogram('counts', process=None, categories=cat) #antes era process=None
    #h.scale(lumi)
    #print(h.values())
    h,_ = plotter.GetData('counts', h)
  else:
    h = plotter.GetHistogram('counts', process=process, categories=cat)
    h.scale(lumi)
  #print(h.values())
  y, e = h.values(sumw2=True)[()]
  #y=h.values(sumw2=True)
  #print(y)
  if process == 'data': e = np.sqrt(y)
  return y, e

def Nonprompt(plt, level='g2jetsg1b', save=False):
  NdataSS, NdataSS_e = GetYieldsNP(p, 'SS', 'data', level)
  NpSS   , NpSS_e    = GetYieldsNP(p, 'SS', 'Prompt', level)
  NnpOS, NnpOS_e = GetYieldsNP(p, 'OS', 'Nonprompt', level)
  NnpSS, NnpSS_e = GetYieldsNP(p, 'SS', 'Nonprompt', level)
  NSS, NSS_e = (NdataSS - NpSS, NdataSS_e + NpSS_e)
  fact   = NnpOS / NnpSS
  fact_e = fact*np.sqrt( (NnpOS_e/NnpOS)**2 + (NnpSS_e/NnpSS)**2 )

  NOS   = NSS*fact
  NOS_e = NOS * np.sqrt( (NSS_e/NSS)**2 + (fact_e/fact)**2 )
  
  SF=NOS/NnpOS
  SF_e=SF* np.sqrt( (NOS_e/NOS)**2 + (NnpOS_e/NnpOS)**2 )
  
  if save:
    t = OutText(outpath, 'Nonprompt_'+level, 'new', 'tex')
    t.bar()
    t.SetTexAlign('l c')
    t.line("Source" + t.vsep() + "$\mathrm{e}^\pm\mu^\mp$"); t.sep()
    t.line("Prompt SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NpSS, NpSS_e));
    t.line("Data SS" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NdataSS, NdataSS_e));
    t.line("Data - prompt (MC) SS" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NSS, NSS_e)); t.sep()
    
    t.line("Nonprompt SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NnpSS, NnpSS_e) );
    t.line("Nonprompt OS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NnpOS, NnpOS_e) );
    t.line("Ratio OS/SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(fact, fact_e) ); t.sep()
    t.line("SF (NOS/NnpOS)"+t.vsep()+"%1.2f $\pm$ %1.2f"%(SF,SF_e))
    t.line("Nonprompt estimate" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NOS, NOS_e) ); t.sep()
    t.write()

  return NOS, NOS_e

def DrawNonprompt(p, var='l0pt', level=level, chan=ch):
  #prdic =  #{'tt': 'TTTo2L2Nu, tbarW, tW', 'DY': 'WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M50, DY_M10to50', 'Nonprompt': 'WJetsToLNu, TTToSemiLep'}
  #colordic = {'tt' : '#ff6782', 'DY':'#00b4bd', 'Nonprompt': '#47ce33'}
  #prdic = {'Prompt': 'TTTo2L2Nu, tbarW, tW, WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M50, DY_M10to50', 'Nonprompt': 'WJetsToLNu, TTToSemiLep'}
  prdic = {'Prompt': 'TTTo2L2Nu, tbarW_noFullHad, tW_noFullHad, WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M_10to50_MLM, DYJetsToLL_M_50_MLM', 'Nonprompt': 'WJetsToLNu, TTToSemiLeptonic'}#,'data': 'MuonEG,DoubleMuon,SingleMuon, EGamma'}
  colordic = {'Prompt' : '#ff6782', 'Nonprompt': '#47ce33'}
  bkgList = ['Prompt','Nonprompt']
  p = plotter(path, prDic=prdic, bkgList=bkgList, colors=colordic, lumi=lumi, var=var)
  p.SetDataName('data')
  categories = {'level':level, 'sign':'SS', 'channel':chan, 'syst':'norm'}
  p.SetRatio(True)
  p.SetOutpath(outpath)
  p.plotData = True
  p.SetLumi(lumi, "fb$^{-1}$", '13 TeV')
  p.SetCategories(categories)
  label = (GetChLab(categories['channel']) if isinstance(categories['channel'], str) else GetChLab(categories['channel'][0]) ) + GetLevLab(categories['level'])  + ', same-sign'
  p.SetRebin(var, 2)
  p.SetRegion(label)
  p.SetOutput(output)
  #p.SetSystematics(syst=['lepSF', 'JES', 'trigSF', 'FSR', 'ISR'])#, 'lepSF']) # FSR, ISR, JES, lepSF, trigSF
  p.Stack(var, xtit='', ytit='', dosyst=True)



if __name__ == '__main__':
  Nonprompt(p, level=level, save=True)
  DrawNonprompt(p, 'njets')


