import numpy as np
import xarray as xr

fn = "/home/guest/work/michele/data/ERA5_Z500_day_djfm_northem_500hPa_2015-2019.nc" #nome del file di partenza

#____NOTE FILTRI APPLICATI____
#Il dataset proviene dalla rianalisi ERA5 e contiene l'altezza di geopotenziale a 500hPa, nei mesi
#di dicembre,gennaio,febbraio e marzo, dal 2015 al 2019 nell'emisfero nord. La griglia in questo casi è di 2.5°

#____IMPORTO I DATI____
#in un oggetto xdataset e memorizzo zg in un xarray
data = xr.load_dataset(fn)
zg = data["zg"]

#____DEF punctual IB____
def pIB(ZG):

	#inizializzo le tuple con i valori dell'altezza di geopotenziale
	ZG_zero = ZG.loc[:,:,30:75,:].values
	ZG_sud = ZG.loc[:,:,15:60,:].values
	ZG_nord = ZG.loc[:,:,45:90,:].values

	#eseguo calcolo GHGS/GHGN per le tuple
	GHGS = (+ ZG_zero - ZG_sud)/15
	GHGN = (- ZG_zero + ZG_nord)/15

	#cerco con la funzione where di xarrays quando GHGS/GHGN soddisfano le
	#condizioni e moltiplico valore per valore le matrici (no prod matriciale)
	BoolGHGS = xr.where(GHGS > 0, 1 , 0)
	BoolGHGN = xr.where(GHGN < -10, 1, 0)
	TuplepIB = BoolGHGN * BoolGHGS

	#Uso il costruttore di un XArray per dare un output nello stesso formato
	#dell'input
	times = ZG.coords["time"].values
	plev = ZG.coords["plev"].values
	lon = ZG.coords["lon"].values
	lat = np.arange(30,75 + 1,2.5) #per non sbagliare la dimensione aggiungo "+1"
	pIB = xr.DataArray(TuplepIB,coords=[times,plev,lat,lon],dims = ZG.dims)
	return pIB


#VERIFICA RISULTATI
Rev = pIB(zg)
print(Rev)
RevBool = Rev.values
Counter = sum(sum(sum(sum(RevBool))))
print("\nNumero di IB puntuali tra il 2015 e il 2019, nell'emisfero nord, nei mesi di dicembre, gennaio, febbraio e marzo: ", Counter)
#print(Rev.shape)
#print(zg.shape)
