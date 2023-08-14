import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': '2m_temperature',
        'year': [
             '2023'
        ],
        'month': [
            '03',
        ],
        'time': [
            '00:00', '06:00', '12:00',
            '18:00',
        ],
        'grid': [
            0.1, 0.1,
         ],        
        'area': [
            45, 105, 25,
            125,
        ],
    },
    'E:\\ncdata_test\\ERA5_2023_monthly_Tem_China.nc')
