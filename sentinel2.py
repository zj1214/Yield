from sentinelsat import SentinelAPI
import os
from osgeo import gdal, ogr, gdalconst
import zipfile
import numpy as np

class Sentinel2():
    def __init__(self, user_name, password, shp, directory_path, start_date, end_date) :
        self.user_name = user_name  # 用户名
        self.password = password  # 密码
        self.shp = shp  # 下载范围
        self.directory_path = directory_path  # 保存目录
        self.start_date = start_date  # 开始日期
        self.end_date = end_date  # 截止日期

    def unzip(self, filepath):
        sen_path = filepath.split(".zip")[0]
        zip_file = zipfile.ZipFile(filepath)
        for f in zip_file.namelist():
            zip_file.extract(f, sen_path)  # 循环解压文件到指定目录
        zip_file.close()
        xml_path = ''
        for root, dirs, files in os.walk(sen_path):
            for file in files:
                file_path = os.path.join(root, file)
                if 'MTD_MSIL2A' in file_path:
                    xml_path = file_path
        print('已解压')
        return xml_path

    def write_img(self, filename, im_proj, im_geotrans, im_data):
        if 'int8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32

        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        else:
            im_bands, (im_height, im_width) = 1, im_data.shape

        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)
        dataset.SetGeoTransform(im_geotrans)
        dataset.SetProjection(im_proj)
        if im_bands == 1:
            dataset.GetRasterBand(1).WriteArray(im_data)
        else:
            for i in range(im_bands):
                dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
        del dataset

    def S2tif(self, filename):
        # 打开栅格数据集
        root_ds = gdal.Open(filename)
        ds_list = root_ds.GetSubDatasets()  # 获取子数据集。该数据以数据集形式存储且以子数据集形式组织
        visual_ds = gdal.Open(ds_list[0][0])  # 打开第1个数据子集的路径。ds_list有4个子集，内部前段是路径，后段是数据信息
        visual_arr = visual_ds.ReadAsArray()
        visual_ds2 = gdal.Open(ds_list[1][0])
        visual_arr2 = visual_ds2.ReadAsArray()

        # 创建.tif文件
        band_count = visual_ds.RasterCount
        xsize = visual_ds.RasterXSize
        ysize = visual_ds.RasterYSize
        out_tif_name = filename.split(".SAFE")[0] + ".tif"
        driver = gdal.GetDriverByName("GTiff")
        out_tif = driver.Create(out_tif_name, xsize, ysize, band_count + 1, gdal.GDT_Float32)
        out_tif.SetProjection(visual_ds.GetProjection())
        out_tif.SetGeoTransform(visual_ds.GetGeoTransform())

        # 获取B11波段
        xsize2 = visual_ds2.RasterXSize
        ysize2 = visual_ds2.RasterYSize
        B11_name = filename.split(".SAFE")[0] + "B11" + ".tif"
        driver2 = gdal.GetDriverByName("GTiff")
        out_tif2 = driver2.Create(B11_name, xsize2, ysize2, 1, gdal.GDT_Float32)
        out_tif2.SetProjection(visual_ds2.GetProjection())
        out_tif2.SetGeoTransform(visual_ds2.GetGeoTransform())
        for index, band in enumerate(visual_arr2):
            band = np.array([band])
            if index == 4:
                out_tif2.GetRasterBand(1).WriteArray(band[0])  # 将波段band2、3、4、8数据写入内存
        out_tif2.FlushCache()  # 最终将数据写入硬盘
        out_tif2 = None  # 关闭tif文件

        # B11波段重采样成10m
        inputrasfile = gdal.Open(B11_name, gdal.GA_ReadOnly)
        inputProj = inputrasfile.GetProjection()
        # 创建重采样输出文件
        B11_reference = filename.split(".SAFE")[0] + "B11_reference" + ".tif"
        driver3 = gdal.GetDriverByName('GTiff')
        output3 = driver3.Create(B11_reference, visual_ds.RasterXSize, visual_ds.RasterYSize, 1, gdal.GDT_Float32)
        output3.SetGeoTransform(visual_ds.GetGeoTransform())
        output3.SetProjection(visual_ds.GetProjection())
        # 参数说明 输入数据集、输出文件、输入投影、参考投影、重采样方法(最邻近内插\双线性内插\三次卷积等)、回调函数
        gdal.ReprojectImage(inputrasfile, output3, inputProj, visual_ds.GetProjection(), gdalconst.GRA_Bilinear, 0.0, 0.0)
        output3 = None

        for index, band in enumerate(visual_arr):
            band = np.array([band])
            for i in range(len(band[:])):
                # 数据写出
                out_tif.GetRasterBand(index + 1).WriteArray(band[i])
        B11_reference_tif = gdal.Open(B11_reference)
        B11_bandarray = B11_reference_tif.ReadAsArray()
        out_tif.GetRasterBand(5).WriteArray(B11_bandarray)
        out_tif.FlushCache()
        out_tif = None

        return (B11_name, B11_reference, out_tif_name)

    def getVIs(self, filepath):
        dataset = gdal.Open(filepath)
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        band_R = im_data[2, 0:im_height, 0:im_width] / 10000
        band_NIR = im_data[3, 0:im_height, 0:im_width] / 10000

        # EVI2
        L = np.ones((im_height, im_width))
        dataMin = (band_NIR - band_R) * 2.50
        dataSum = (band_NIR + 2.4 * L * band_R + L) * 1.00
        Mask = dataSum == 0
        dataSum[Mask] = np.nan
        EVI2 = dataMin * 1.00 / dataSum * 1.00
        outpath = filepath.split(".tif")[0] + '_EVI2.tif'
        self.write_img(outpath, im_proj, im_geotrans, EVI2)
        print('EVI2指数计算完成：' + str(outpath))

    def download_sentinel_data(self):
        shape = ogr.Open(self.shp)
        lyr = shape.GetLayer()
        count = lyr.GetFeatureCount()
        points = []
        for i in range(0, count):
            poly = lyr.GetFeature(i)
            geometry = poly.GetGeometryRef()
            points = points + geometry.ExportToWkt().replace('POLYGON ((', '').replace('))', '').split(",")
        lon = []
        lat = []
        for i in points:
            lon.append(float(i.split(" ")[0]))
            lat.append(float(i.split(" ")[1]))
        footprint = str(
            'POLYGON((' + str(min(lon)) + ' ' + str(min(lat)) + ',' + str(min(lon)) + ' ' + str(max(lat)) + ',' + str(
                max(lon)) + ' ' + str(max(lat)) + ',' + str(max(lon)) + ' ' + str(min(lat)) + ',' + str(
                min(lon)) + ' ' + str(min(lat)) + '))')

        # 创建SentinelAPI，请使用哥白尼数据开放获取中心自己的用户名及密码
        api = SentinelAPI(self.user_name, self.password, 'https://scihub.copernicus.eu/dhus/')
        products = api.query(footprint,
                             date=(self.start_date, self.end_date),
                             platformname='Sentinel-2',
                             producttype='S2MSI2A',
                             cloudcoverpercentage=(0, 30))
        print('搜索到影像：' + str(len(products)) + '幅')
        for product in products:
            product_info = api.get_product_odata(product)
            print('正在下载：' + product_info['title'])
            api.download(product, self.directory_path)

            file = self.directory_path + '\\' + str(product_info['title']) + '.zip'
            xml_path = self.unzip(file)
            removefile = self.S2tif(xml_path)
            if os.path.exists(str(removefile[0])):
                os.remove(str(removefile[0]))
                os.remove(str(removefile[1]))
            self.getVIs(removefile[2])

if __name__ == "__main__":
    user_name = 'z1214'
    password = 'mirrors03.'
    shp = "E:\潍柴雷沃估产数据\聊城SHP\聊城农场.shp"  # shp路径
    directory_path = "E:\潍柴雷沃估产数据\s2"  # 保存目录
    start_date = '20220428'
    end_date = '20220520'

    e = Sentinel2(user_name, password, shp, directory_path, start_date, end_date)
    e.download_sentinel_data()
