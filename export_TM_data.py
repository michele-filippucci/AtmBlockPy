from lib.BlockTools import BlockTools
import lib.BlockTools as BT

input="/home/guest/work/michele/data/ERA5/Z500/"+\
     "ERA5_Z500_day_djfm_r144x73_500hPa_northem_1979-2008.nc"
output="/home/guest/work/michele/prog/plots/dataset.dat"

tm = BlockTools()
print("Reading data:")
print(tm.read(input))
print("OK")
print("Creating dataset:")
print(tm.TM(output))
