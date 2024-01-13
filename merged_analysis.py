# -*- coding: utf-8 -*-
import os, sys, subprocess
#funciona cambiando la variable ejecutar con las muestras que se quieren y el year
#muestras de 2016
muones=['DoubleMuon_Run2016F','DoubleMuon_Run2016G', 'DoubleMuon_Run2016H']
electrones=['DoubleEG_Run2016F','DoubleEG_Run2016G','DoubleEG_Run2016H']
muel=['MuonEG_Run2016F','MuonEG_Run2016G','MuonEG_Run2016H']
muon=['SingleMuon_Run2016F','SingleMuon_Run2016G','SingleMuon_Run2016H']
elec=['SingleElectron_Run2016F','SingleElectron_Run2016G','SingleElectron_Run2016H']
ejecutar16=[muones,electrones,muel,muon,elec]

#muestras de 2016apv
muonesapv=['DoubleMuon_Run2016B','DoubleMuon_Run2016C', 'DoubleMuon_Run2016D','DoubleMuon_Run2016E','DoubleMuon_Run2016F']
electronesapv=['DoubleEG_Run2016B','DoubleEG_Run2016C','DoubleEG_Run2016D','DoubleEG_Run2016E','DoubleEG_Run2016F']
muelapv=['MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F']
muonapv=['SingleMuon_Run2016B','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F']
elecapv=['SingleElectron_Run2016B','SingleElectron_Run2016C','SingleElectron_Run2016D','SingleElectron_Run2016E','SingleElectron_Run2016F']
ejecutarapv=[muonesapv,electronesapv,muelapv,muonapv,elecapv]

#muestras de 2017
muones17=['DoubleMuon_Run2017B','DoubleMuon_Run2017C', 'DoubleMuon_Run2017D','DoubleMuon_Run2017E','DoubleMuon_Run2017F']
electrones17=['DoubleEG_Run2017B','DoubleEG_Run2017C','DoubleEG_Run2017D','DoubleEG_Run2017E','DoubleEG_Run2017F']
muel17=['MuonEG_Run2017B','MuonEG_Run2017C','MuonEG_Run2017D','MuonEG_Run2017E','MuonEG_Run2017F']
muon17=['SingleMuon_Run2017B','SingleMuon_Run2017C','SingleMuon_Run2017D','SingleMuon_Run2017E','SingleMuon_Run2017F']
elec17=['SingleElectron_Run2017B','SingleElectron_Run2017C','SingleElectron_Run2017D','SingleElectron_Run2017E','SingleElectron_Run2017F']
ejecutar17=[muones17,electrones17,muel17,muon17,elec17]

#muestras de 2018
muones18=['DoubleMuon_Run2018A','DoubleMuon_Run2018B','DoubleMuon_Run2018C','DoubleMuon_Run2018D']
elec18=['EGamma_Run2018A','EGamma_Run2018B','EGamma_Run2018C','EGamma_Run2018D']
muon18=['SingleMuon_Run2018A','SingleMuon_Run2018B','SingleMuon_Run2018C','SingleMuon_Run2018D']
muel18=['MuonEG_Run2018A','MuonEG_Run2018B','MuonEG_Run2018C','MuonEG_Run2018D']
ejecutar18=[muones18,elec18,muon18,muel18]


#Montecarlo (igual para todos)
MC1=['TTTo2L2Nu', 'tW_noFullHad', 'tbarW_noFullHad','WJetsToLNu_MLM']
MC2=['DYJetsToLL_M_10to50_MLM', 'DYJetsToLL_M_50_MLM']
MC3=['TTToSemiLeptonic','WWTo2L2Nu'] #'WZTo3LNu'
MC4=['WZ','ZZ'] #'WW','ZZTo4L','ZZTo2L2Nu'
MC5=['TTG', 'TTZToLL_M_1to10', 'TTZToLLNuNu_M_10','TTWJetsToLNu','TTWJetsToQQ', 'TTZToQQ']
MC6=['WWW_ext1', 'WWZ_ext1', 'WZZ_ext1','ZZZ_ext1'] #voy a quitar 'WWW','WWZ', 'WZZ', 'ZZZ'

MC_s=['tW_noFullHad', 'tbarW_noFullHad','WJetsToLNu_MLM','DYJetsToLL_M_10to50_MLM','TTToSemiLeptonic','WWTo2L2Nu','WZTo3LNu','WW','WZ','ZZ']
MC_l=['TTTo2L2Nu','DYJetsToLL_M_50_MLM','ZZTo2L2Nu']

#Montecarlo para las incertidumbres
MC_teor=['TTTo2L2Nu_TuneCP5down','TTTo2L2Nu_TuneCP5up','TTTo2L2Nu_hdampDOWN','TTTo2L2Nu_hdampUP']#'TTTo2L2Nu_mtop169p5','TTTo2L2Nu_mtop175p5',

#ejecutar=[MC1,MC2]
#ejecutar=[MC1,MC2,MC3]
#ejecutar=[MC3]
#ejecutar=[MC1,MC2,MC3,MC4,MC5,MC6]#MC1,MC2,
#ejecutar=[['TTTo2L2Nu']]
#ejecutar=[['tW_noFullHad']]
#ejecutar=[MC1,MC2,MC3,MC4]#+ejecutar16
#ejecutar=[MC5,MC6]
#ejecutar=ejecutar18
ejecutar=ejecutarapv
#ejecutar=[MC_teor]
ncores=16 #32
batch='-j'
#batch=''
#json='json/2016/TTTo2L2Nu_2016.json'
per='2016apv'
var=False
if ejecutar==[MC_teor]: var=True
else: var=False

def correr(lista,ncores,per,batch):
	'''
	función para correr a la vez varios jsons en función del periodo y la muestra
	'''
	for i in lista:
		if var==True:
			json='json/'+per+'/'+'variaciones/'+ i+'_'+per+'.json' 
			#output='histosprueba/'+per+'/variaciones/'
			output='histosprueba/SF/trigger/'+per+'/variaciones/'
			#output='histosprueba/Run2/variaciones/'
		else:
			json='json/'+per+'/'+ i+'_'+per+'.json' 
			#output='histosprueba/'+per
			output='histosprueba/SF/trigger/'+per
			#output='histosprueba/Run2/'
		if 'Run' in i:
			if per=='2016apv':
				nombre=i+'apv'
			else:
				nombre=i
		else:
			nombre=i+'_'+ per
		os.system('python analysis/tt3_prueba/run.py %s -n %s %s -p %s -o %s' %(json, ncores, batch,output,nombre))
for i in ejecutar:
	correr(i,ncores,per,batch)
