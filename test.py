import netCDF4
import os

view = netCDF4.Dataset(os.path.join(".", "XAERDT_L3_MEASURES_QD_HH.A2022365.2030.001.2024316153136.nc"), "r+")

print(view.variables.keys())
