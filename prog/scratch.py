import netCDF4 as nc
import numpy as np
import xarray as xr

fn = "/home/guest/work/michele/data/ERA5_Z500_day_djfm_northem_500hPa_2015-2019.nc" #nome del file di partenza
data = nc.Dataset(fn) #oggetto dataset di netCDF4
#per avere info sui tipi di dati presento da console apro il dataset. Es arrivato a questo punto chiamo data dalla riga di comando di python
zg_info = data["zg"] #se assegno così ho le info sui data
print(zg_info)
zg = data["zg"][:] #se assegno così ho il dato stesso che in questo caso è un array 4-d di numpy
print("\n\n")
print(zg[100,0,20,70])

#SINTASSI per l'altezza di geopotenziale: zg[time,plev(?),latitudine,longitudine]
#for i in range(0,100):
#	print(zg[100,0,20,i])

#ALTERNATIVA: utilizzare xarray

data = xr.load_dataset(fn)
zg = data["zg"] #posso chiamare l'etichetta
#print(zg)
print("\n\n")
#print(zg("lat"))
print(zg[:,0,20,70])
#zg = zg*2
#print(zg[100,0,20,70]) #stampa molte info ma posso anche effettuare operazioni algebriche
#for i in range(0,100):
#	print(zg[100,0,20,i])

