import os
from config import *

if var is None: var = "counts" #l0pt

outpath = path + "/combineFiles/"
if not os.path.exists(outpath):
    os.makedirs(outpath)

#pathcomb = "/nfs/fanae/user/juanr/CMSSW_10_2_13/src/combinescripts/tt5TeVljets/"+datatoday + "/" 
pathcomb="/nfs/fanae/user/jorgeat/Documentos/tfm/cafea/pruebacombine/final2/"

if not os.path.exists(pathcomb):
    os.makedirs(pathcomb)

#levels = ['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
#channels = ['e', 'm']
levels=['g2jetsg1b']#['dilep','g2jetsg1b','g2jetsg2b']
channels=['em']#, 'mm', 'ee']
iteration = 0
total = len(channels)*len(levels)
verbose=1
for ch in channels:
    for level in levels:
        print("\r[{:<100}] {:.2f} % ".format('#' * int(float(iteration)/total*100), float(iteration)/total*100),end='')
        var='counts'
        outname = "%s_%s_%s_%s.root"%(var, ch, level,year)
        #if not os.path.exists(f"{pathcomb+outname}"): #quitar esto cuando se cambien los root
        command = "python analysis/tt3_prueba/SaveRootfile.py -p %s -v %s -l %s -c %s --data"%(path, var, level, ch)
        if verbose: print("Running: %s"%(command))
        os.system(command)

        # Move the file to the combine folder
        mvcommand = "cp %s %s" %(outpath+outname, pathcomb)
        if verbose: print("Running: %s"%(mvcommand))
        os.system(mvcommand)

        # Create the datacard
        cardcommand = f"python analysis/tt3_prueba/CreateDatacard.py -p {path} --inputFile {pathcomb+outname}"
        if verbose: print("Running: %s"%(cardcommand))
        os.system(cardcommand)
        if verbose >= 1: print("  ")
        iteration+= 1
print("\r[{:<100}] {:.2f} % ".format('#' * int(float(iteration)/total*100), float(iteration)/total*100))
print('Datacards created in ', pathcomb)
