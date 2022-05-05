import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': 'total_column_water_vapour',
        'year': '2021',
        'month': '03',
        'day': '19',
        'time': '03:00',
        'area': [
            26, 96, 23,
            99,
        ],
    },
    'tcwv.nc')