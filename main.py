import json
import os

# import netCDF4
import numpy as np
import pandas
import requests
from netCDF4 import Dataset

json_loc = os.path.join(".", "site_map.json")
net_loc = os.path.join(
    ".", "XAERDT_L3_MEASURES_QD_HH.A2022365.2030.001.2024316153136.nc"
)

test_loc = os.path.join(
    ".", "XAERDT_L3_MEASURES_QD_HH.A2022365.2030.001.2024316153136.nc"
)


def make_df():
    site_list_url = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
    need2update = False
    site_df = None
    site_list = requests.get(site_list_url)
    site_list = site_list.text
    loc = os.path.join(".", "aeronet_locations_v3.txt")
    if site_list is not None and os.path.exists(loc):
        with open(loc, "r") as file:
            if file.read() == site_list:
                need2update = False
            else:
                need2update = True
            file.close()
    else:
        need2update = True

    if need2update:
        with open(loc, "w") as file:
            file.write(site_list)
        file.close()

    return pandas.read_csv(loc, skiprows=1, delimiter=",")


def getNGP(lat, lon, site_lat, site_lon):
    R = 6371000  # radius of the earth in meters
    lat1 = np.radians(site_lat)
    lat2 = np.radians(lat)
    delta_lat = np.radians(lat - site_lat)
    delta_lon = np.radians(lon - site_lon)

    # print((np.sin(delta_lat / 2)))
    # print((np.sin(delta_lat / 2)))
    # print((np.cos(lat1)))
    # print((np.cos(lat2)))
    # print((np.sin(delta_lon / 2)))
    # print(np.sin(delta_lon / 2))
    a = (np.sin(delta_lat / 2)) * (np.sin(delta_lat / 2)) + (np.cos(lat1)) * (
        np.cos(lat2)
    ) * (np.sin(delta_lon / 2)) * (np.sin(delta_lon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = R * c
    return np.unravel_index(d.argmin(), d.shape)


def process(site_df):
    # sites = None
    # site_df = pandas.read_csv(loc, skiprows=1, delimiter=",") if site_df is None else site_df
    geods = Dataset(net_loc, "r")
    # print(geods.variables.keys())
    # var_list = [var for geods.variables.keys() if "lat" in geod]

    # for variable in geods.variables.keys():
    #     print(variable, type(geods[variable]))
    lat = geods.variables["lat"][:][:]
    lon = geods.variables["lon"][:][:]

    print()

    for _, site in site_df.iterrows():

        for variable in geods.variables.keys():
            if "AOD" in variable:
                # print(type(geods[variable]))
                # print(geods[variable])
                # #             # print(geods[variable].keys())
                # #             # lon = geods[variable].variables["longitude"][:][:]
                # #             # lat = geods[variable].variables["latitude"][:][:]
                # #             # print(lat, lon)
                site_lat = site["Latitude(decimal_degrees)"]
                site_lon = site["Longitude(decimal_degrees)"]
                # gets (and then prints) the x,y location of the nearest point in data to entered location, accounting for no data values
                x, y = getNGP(lat, lon, site_lat, site_lon)
        #
        # # site_df.loc[site["Site_Name"], "Longitude(decimal_degrees)"] = x
        # # site_df.loc[site["Site_Name"], "Latitude(decimal_degrees)"] = y
        # print(x, y)


## NEXT Run Calculations against new lat, lngs in the data frame against netcdf

## SAVE for export


if __name__ == "__main__":
    site_df = make_df()
    process(site_df)
