# -*- coding: utf-8 -*-
"""SunPy sample data files"""
from __future__ import absolute_import

from os import remove
import os.path
from zipfile import ZipFile
from urllib2 import URLError
from shutil import move

from astropy.utils.data import download_file

from sunpy.util.net import url_exists
from sunpy import config

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"


sampledata_dir = config.get("downloads", "sample_dir")

# urls to search for the sample data
_base_urls = (
    'http://data.sunpy.org/solarbextrapolation/sample_data/',
    'http://hesperia.gsfc.nasa.gov/~schriste/sunpy-sample-data/',
    'https://github.com/ehsteve/sunpy-sample-data/raw/master/')

# keys are file shortcuts
# values consist of filename as well as optional file extension if files are
# hosted compressed. This extension is removed after download.
_files = {
    "HMI_2011": ("2011-02-14__20-35-25__01_hmi.fits", ""),
    "AIA_2011": ("2011-02-14__20-35-25__02_aia.fits", ""),
    "Bxyz_2011": ("2011-02-14__20-35-25__03_Bxyz.npy", ""),
}

sample_files = {}
for key in _files:
    sample_files[key] = os.path.abspath(os.path.join(sampledata_dir, _files[key][0]))


def download_sample_data(progress=True, overwrite=True):
    """
    Download the sample data.

    Parameters
    ----------
    progress: bool
        Show a progress bar during download
    overwrite: bool
        If exist overwrites the downloaded sample data.

    Returns
    -------
    None
    """
    number_of_files_fetched = 0
    print("Downloading sample files to " + sampledata_dir)
    for file_name in _files.itervalues():
        if not overwrite:
            if os.path.isfile(os.path.join(sampledata_dir,
                                           file_name[0])):
                number_of_files_fetched += 1
                continue

        for base_url in _base_urls:
            full_file_name = file_name[0] + file_name[1]
            print(full_file_name)
            if url_exists(os.path.join(base_url, full_file_name)):
                f = download_file(os.path.join(base_url, full_file_name))
                real_name, ext = os.path.splitext(full_file_name)

                if file_name[1] == '.zip':
                    print("Unpacking: %s" % real_name)
                    with ZipFile(f, 'r') as zip_file:
                        zip_file.extract(real_name, sampledata_dir)
                    remove(f)
                else:
                    # move files to the data directory
                    move(f, os.path.join(sampledata_dir, file_name[0]))
                # increment the number of files obtained to check later
                number_of_files_fetched += 1
                break

    if number_of_files_fetched < len(_files.keys()):
        raise URLError("Could not download all samples files. Problem with accessing sample data servers.")
