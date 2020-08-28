from netCDF4 import Dataset
fh=Dataset("subSets/subSet_wrfout_d04_20140611_15:00-16:00.nc","r")
dbz=fh.variables.get("dbz")

