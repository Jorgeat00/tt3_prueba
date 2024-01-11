#'codigo' para definir el periodo y la luminosidad de las muestras
import json
#2016, 2016apv,2017,2018
per='2016apv'
def muestras(var):
	return per


def get_lumi(year):
    lumi_json = "lumi.json"
    year=year.upper()
    with open(lumi_json) as f_lumi:
       lumi = json.load(f_lumi)
       lumi = lumi[year]
    return lumi

