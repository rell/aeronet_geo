"""

Author: Terrell Credle

Generates csv data from a given list of NetCDF files for geo data for all aeronet locations.

"""

import json
import os
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import numpy.ma as ma
import pandas
import requests
from netCDF4 import Dataset
from tqdm.auto import tqdm

progress_bar = None  # Global progress for tracking the processing of sites
AERONET_FILE_LOC = os.path.join(".", "aeronet_locations_v3.txt")
PREPROCESS_LOC = os.path.join(".", "aeronet_locations.csv")


# Creates an array of all NC files from a given list
def collect_nc4s():
    try:
        with open("filelist.txt", "r") as file:
            files = file.readlines()
            files = [file.strip() for file in files if file.strip()]
            return files
    except:
        print("no files list found")


# Testing short list of locations for developement
def make_df_test():
    AERONET_FILE_LOC = os.path.join(".", "aeronet_locations_v3_test.txt")
    return pandas.read_csv(AERONET_FILE_LOC, skiprows=1, delimiter=",")


# Builds site_df from aeronet location list
def make_df():
    site_list_url = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
    need2update = False
    site_df = None
    site_list = requests.get(site_list_url)
    site_list = site_list.text
    if site_list is not None and os.path.exists(AERONET_FILE_LOC):
        with open(AERONET_FILE_LOC, "r") as file:
            if file.read() == site_list:
                need2update = False
                print("Location list is already upto date, moving to processing.")
            else:
                need2update = True
            file.close()
    else:
        need2update = True

    if need2update:
        print("Updating location list\n")
        with open(AERONET_FILE_LOC, "w") as file:
            file.write(site_list)
        file.close()

    return (need2update, pandas.read_csv(AERONET_FILE_LOC, skiprows=1, delimiter=","))


def get_updated_df():
    return pandas.read_csv(PREPROCESS_LOC, delimiter=",")


def getNGP(lat, lon, site_lat, site_lon):
    """
    NGP with haversine formula from pawan's repo:
        https://github.com/pawanpgupta/DTAerosols/blob/main/read_viirs_at_a_location.ipynb
    """
    R = 6371000
    lat1 = np.radians(site_lat)
    lat2 = np.radians(lat)
    delta_lat = np.radians(lat - site_lat)
    delta_lon = np.radians(lon - site_lon)
    a = (np.sin(delta_lat / 2)) * (np.sin(delta_lat / 2)) + (np.cos(lat1)) * (
        np.cos(lat2)
    ) * (np.sin(delta_lon / 2)) * (np.sin(delta_lon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = R * c
    return np.unravel_index(d.argmin(), d.shape)


def build_aeronet_df(site, ds_path, new_lat, new_lon):
    try:
        geods = Dataset(ds_path, "r")
        index, site = site
        # TODO: possibly need to use netcdf  lat lng instead of norm lat lng
        site_lat = site["Latitude(decimal_degrees)"]
        site_lon = site["Longitude(decimal_degrees)"]
        site_data = {"Site_Name": site["Site_Name"]}
        x, y = getNGP(new_lat, new_lon, site_lat, site_lon)
        site_data["NGP_lat"] = x
        site_data["NGP_lon"] = y
        site_data["netCDF_lat"] = new_lat[x, y]
        site_data["netCDF_lon"] = new_lon[x, y]
        site_data["lat"] = site_lat
        site_data["lon"] = site_lon

        update_progress()
        return site_data

    except Exception as e:
        print(e)


# Just update global progress bar for tracking sites
def update_progress(n=1):
    global progress_bar
    if progress_bar is not None:
        progress_bar.update(n)


# Build variable list for variable exclusion to prevent errors
def get_geodslist(geods):
    list_of_vars = []
    for variable in geods.variables.keys():
        if geods.variables[variable][:].ndim >= 2:
            list_of_vars.append(variable)
    return list_of_vars


def process_aeronet_list(file, site_df):
    ds_path = os.path.join(".", file)
    geods = Dataset(ds_path, "r")

    # create new lat lon dims
    lat = geods.variables["lat"][:]
    lon = geods.variables["lon"][:]

    new_lat = np.zeros((720, 1440), dtype=float)
    new_lon = np.zeros((720, 1440), dtype=float)

    for i in range(0, 1440):
        new_lat[:, i] = lat

    for i in range(0, 720):
        new_lon[i, :] = lon

    process_site_partial = partial(
        build_aeronet_df,
        ds_path=ds_path,
        new_lat=new_lat,
        new_lon=new_lon,
    )

    with Pool(7) as pool:
        site_data_list = list(
            tqdm(
                pool.imap_unordered(
                    process_site_partial,
                    site_df.iterrows(),
                    chunksize=50,
                ),
                total=site_df.shape[1],
                desc="Processing Locations",
            )
        )

    df = pandas.DataFrame(site_data_list)
    df.to_csv(PREPROCESS_LOC, index=False)


# Grab all netCDF output of all vars for a given site
def process_site(site, ds_path, var_list):
    try:
        geods = Dataset(ds_path, "r")
        index, site = site

        site_data = {"Site_Name": site["Site_Name"]}
        for variable in var_list:
            try:
                # Set fill to prevent numpy from messing up fill netCDF fill value
                data = geods.variables[variable][
                    ..., site["netCDF_lat"], site["netCDF_lon"]
                ]
                ma.set_fill_value(data, -9999)
                cleaned_data = data.filled()
                site_data[variable] = ma.getdata(cleaned_data)
            except Exception as e:
                pass

        update_progress()
        return site_data

    except Exception as e:
        print(e)


def process(site_df, files):
    save_path = os.path.join(".", "csv")
    os.makedirs(save_path, exist_ok=True)

    for ds in tqdm(files, desc="Processing Files"):
        save_csv = ds.split(".nc")[0] + ".csv"
        ds_path = os.path.join(".", ds)
        geods = Dataset(ds_path, "r")
        list_of_vars = get_geodslist(geods)

        process_site_partial = partial(
            process_site,
            ds_path=ds_path,
            var_list=list_of_vars,
        )

        # Safe pool sizing to prevent crashing of server (Possible to increase to a factor of system thread count)
        with Pool(8) as pool:
            site_data_list = list(
                tqdm(
                    pool.imap_unordered(
                        process_site_partial,
                        site_df.iterrows(),
                        chunksize=50,
                    ),
                    total=site_df.shape[0],
                    desc="Processing Locations",
                )
            )

        # Fixes previous issue of attempting to directly add masked array to dataframe
        flattened_data = [
            {
                key: (
                    value[0]
                    if isinstance(value, np.ndarray)
                    and hasattr(value, "__len__")
                    and len(value) == 1
                    else value
                )
                for key, value in record.items()
            }
            for record in site_data_list
        ]

        df = pandas.DataFrame(flattened_data)

        completed_df = pandas.merge(site_df, df, on="Site_Name", how="right")
        completed_df.to_csv(os.path.join(save_path, save_csv), index=False)


if __name__ == "__main__":
    files = collect_nc4s()
    site_df = make_df()
    if site_df[0] or not os.path.exists(PREPROCESS_LOC):
        print(
            "aeronet preprocess csv is not generated or needs to be updated... generating now."
        )
        process_aeronet_list(files[0], site_df[1])
    else:
        print("aeronet preprocess csv is updated and does not need to be generated.")

    site_df = get_updated_df()
    # site_df = make_df_test()
    process(site_df, files)
