import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'surface_net_solar_radiation',
        'year': '2022',
        'month': [
            '04','05'
        ],
        'time': '00:00',
        'area': [
            45, 105, 25,
            125,
        ],
    },
    'E:\\ncdata_test\\ERA5_2022_monthly_NetRad_China.nc')