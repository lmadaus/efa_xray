#!/usr/bin/env python
from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date

def get_metar_observations(siteids=None, obtypes=['air_temperature'], use_times=[datetime.utcnow()-timedelta(hours=6), datetime.utcnow()], bounds=None):
    """
    Returns METAR observations
    """
    from siphon.catalog import TDSCatalog
    from siphon.ncss import NCSS

    # copied from the browser url box
    metar_cat_url = 'http://thredds.ucar.edu/thredds/catalog/nws/metar/ncdecoded/catalog.html?dataset=nws/metar/ncdecoded/Metar_Station_Data_fc.cdmr'
    # parse the xml
    metar_cat = TDSCatalog(metar_cat_url)
    # what datasets are here? only one "dataset" in this catalog
    dataset = list(metar_cat.datasets.values())[0]
    print(dataset.name)
    ncss_url = dataset.access_urls['NetcdfSubset']

    ncss = NCSS(ncss_url)
    query = ncss.query()

    # Build the query
    query.time_range(use_times[0], use_times[-1])
    if bounds is not None:
        query.lonlat_box(*bounds)
    query.variables(*obtypes)
    query.accept("netcdf")

    data = ncss.get_data(query)

    print data

if __name__ == '__main__':
    get_metar_observations()


