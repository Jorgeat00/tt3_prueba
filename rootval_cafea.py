# -*- coding: utf-8 -*-
#from cafea.plotter.plotter import GetXYfromH1D
import uproot3
from config import *

#hay que hacer esto función para que coja bien los valores


per=year
abspath='/nfs/fanae/user/jorgeat/Documentos/tfm/cafea/'
#path=abspath+'histosprueba/'+per+'/combineFiles/'
path=abspath+'pruebacombine/'
histogramas=['DY','VV','Non-prompt', 'tW', 'ttX','tt']
incertidumbres=['lepSF_muon','lepSF_elec','trigSF','btagging','JES','JER','PU','PreFiring']#'PDF','Scales','hdamp','UE']
#incertt=['ISR','FSR']
incertt=[]

datos=[]
canal='em'
nivel='g2jetsg1b'

treename='counts_'+canal+'_'+nivel
print(per)
#da=r.TFile(path+treename+'_'+per+'.root')
#da=uproot3.open(path+'prueba_merged.root')
da=uproot3.open(path+treename+'_'+per+'.root')
#print(da)
#print(da['tt'].values)


datrootup=[]
datrootdo=[]
datroot=[]
for i in histogramas:
	if i=='tt':
		incertidumbres += incertt
	else:
		incertidumbres=incertidumbres
	for j in incertidumbres:
		varup=i+'_'+j+'Up'
		vardo=i+'_'+j+'Down'
		var=i
		#print(var)
		d=da[var].values
		dup=da[varup].values
		ddo=da[vardo].values
		datrootup.append(float(dup))
		datrootdo.append(float(ddo))
		datroot.append(d)
'''
print('variaciones')
print('nom')
print(datroot)
print('up')
print(datrootup)
print('down')
print(datrootdo)
'''
vartotup=[]
vartotdo=[]
vartot=[]
for i in range(len(datroot)):
	nom=datroot[i]
	varup=datrootup[i]
	vardo=datrootdo[i]
	difup=abs(varup-nom)/nom
	difdo=abs(vardo-nom)/nom
	vartotup.append(difup)
	vartotdo.append(difdo)
	varmax=max(difup,difdo)
	vartot.append(varmax)
	

'''
for i in range(len(vartotdo)):
	print(abs(vartotup[i]-vartotdo[i]))
'''
'''
print('variaciones finales')
print('vartotup')
print('-------------')
print(vartotup)
print('vartotdo')
print('----------------------')
print(vartotdo)
'''






'''
scales=uproot3.recreate(path+treename+'_'+per+'.root')
with scales as file:
	var='DY'
	file[var].values=file[var].values*2
		
		
print('ya esta')

'''
