import os
import math
import random
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from osgeo import gdal, gdalnumeric, ogr, osr, gdal_array
from PIL import Image,ImageDraw
from scipy.optimize import minimize
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib as mpl
gdal.UseExceptions()
from datetime import datetime, timedelta

class Cliptif():
    def __init__(self, shp, img, out) :
        self.shp = shp # 输入shp
        self.img = img # 待裁剪tif
        self.out = out # 输出tif

    def WGStoUTM(self, lon, lat,EPSG):
        '''
        WGS坐标转UTM坐标
        '''
        source = osr.SpatialReference()
        source.ImportFromEPSG(4326) # WGS84
        target = osr.SpatialReference()
        target.ImportFromEPSG(EPSG) # UTM
        ct = osr.CoordinateTransformation(source, target)
        pt = ct.TransformPoint(lat, lon)
        return pt

    def shpExtent(self, dataset, points,EPSG):
        '''
        返回点集合在栅格内的UTM投影的坐标信息
        '''
        geoTrans = dataset.GetGeoTransform()
        # [0]横/[3]纵坐标  [1]水平/[4]垂直分辨率
        ulX_min = geoTrans[0]
        ulY_max = geoTrans[3]
        ulX_max = ulX_min + dataset.RasterXSize * geoTrans[1]
        ulY_min = ulY_max + dataset.RasterYSize * geoTrans[5]
        x = []
        y = []
        for i in points:
            point = self.WGStoUTM(float(i.split(" ")[0]), float(i.split(" ")[1]),EPSG) # UTM投影的点坐标
            # 如果点在栅格内
            if point[0] >= ulX_min and point[0] <= ulX_max:
                if point[1] >= ulY_min and point[1] <= ulY_max:
                    x.append(point[0])
                    y.append(point[1])
            # 如果点不在栅格内
            if point[0] < ulX_min:
                x.append(ulX_min)
            if point[0] > ulX_max:
                x.append(ulX_max)
            if point[1] < ulY_min:
                y.append(ulY_min)
            if point[1] > ulY_max:
                y.append(ulY_max)
        return (min(x), min(y), max(x), max(y))

    def world2Pixel(self, geoMatrix, x, y):
        '''
        返回输入xy信息在栅格数据中对应坐标位置
        '''
        # [0]横/[3]纵坐标  [1]水平/[5]垂直分辨率  [2]行/[4]列旋转
        ulX = geoMatrix[0]
        ulY = geoMatrix[3]
        xDist = geoMatrix[1]
        yDist = geoMatrix[5]
        pixel = int((x - ulX) / xDist)
        line = int((ulY - y) / abs(yDist))
        return (pixel, line)

    def image2Array(self, i):
        a = gdal_array.numpy.frombuffer(i.tobytes(), 'b')
        a.shape = i.im.size[1], i.im.size[0]
        return a

    def write_img(self, filename, im_proj, im_geotrans, im_data):
        '''
        创建栅格影像
        '''
        # 数据类型
        if 'int8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32
        # 判断波段信息
        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        else:
            im_bands, (im_height, im_width) = 1, im_data.shape
        
        # 创建tif模板
        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)
        dataset.SetGeoTransform(im_geotrans)
        dataset.SetProjection(im_proj)
        # 依次写入波段栅格信息
        if im_bands == 1:
            dataset.GetRasterBand(1).WriteArray(im_data)
        else:
            for i in range(im_bands):
                dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
        del dataset

    def main(self):
        srcImage = gdal.Open(self.img)
        geoTrans = srcImage.GetGeoTransform()
        geoProj = srcImage.GetProjection()
        EPSG = int(str(geoProj).split('"')[-2])
        shape = ogr.Open(self.shp)
        lyr = shape.GetLayer()
        count = lyr.GetFeatureCount()
        points = [] 
        for i in range(0, count):
            poly = lyr.GetFeature(i)
            geometry = poly.GetGeometryRef()
            points = points + geometry.ExportToWkt().replace('POLYGON ((', '').replace('))', '').split(",") # shp边界的点集合
        minX, minY, maxX, maxY = self.shpExtent(srcImage, points,EPSG) # 点集合的外接矩形的坐标点
        # print(minX, minY, maxX, maxY)
        ulX, ulY = self.world2Pixel(geoTrans, minX, maxY) # 对应栅格数据中起点坐标
        lrX, lrY = self.world2Pixel(geoTrans, maxX, minY) # 对应栅格数据中终点坐标
        pxWidth = abs(int(lrX - ulX))
        pxHeight = abs(int(lrY - ulY))
        # print(ulX, ulY, pxWidth, pxHeight)
        _clip = srcImage.ReadAsArray(ulX, ulY, pxWidth, pxHeight) # 读取shp范围内的栅格数据

        # 创建新的geoMatrix
        geoTrans = list(geoTrans)
        geoTrans[0] = minX
        geoTrans[3] = maxY
        # 创建掩膜图层对图像进行裁剪
        rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
        rasterize = ImageDraw.Draw(rasterPoly)
        for i in range(0, count):
            pixels = []
            poly = lyr.GetFeature(i) # shp的第i个要素
            geometry = poly.GetGeometryRef()
            pts = geometry.GetGeometryRef(0) # 该要素的点集合
            for p in range(pts.GetPointCount()):
                # 点在栅格上的坐标位置
                pixels.append(
                    self.world2Pixel(geoTrans, self.WGStoUTM(pts.GetX(p), pts.GetY(p),EPSG)[0], self.WGStoUTM(pts.GetX(p), pts.GetY(p),EPSG)[1]))
            rasterize.polygon(pixels, 0)
        mask = self.image2Array(rasterPoly) # 掩膜图像
        clip = gdalnumeric.choose(mask, (_clip, 0))
        clip[clip == 0.0] = np.nan
        self.write_img(self.out, geoProj, geoTrans, clip)
        gdal.ErrorReset()

class Stack():
    def __init__(self, nc_path, shp_path,eviimg,out_path) :
        self.nc_path = nc_path # nc文件路径
        self.shp_path = shp_path # 农场shp路径
        self.eviimg = eviimg # s2数据路径
        self.out_path = out_path # 波段合成后导出路径
        self.e = Cliptif(self.shp_path, self.eviimg, self.out_path)

    def getvalue(self):
        '''
        获取shp对应的气象数据 返回排序好的气象数据列表
        '''
        data = ogr.Open(self.shp_path)
        layer = data.GetLayer()
        extent = layer.GetExtent()
        # 获取shp中心经纬度
        lon = (extent[0] + extent[1]) / 2
        lat = (extent[2] + extent[3]) / 2

        pathlist = []
        for root, dirs, files in os.walk(self.nc_path):
            if len(files) != 0:
                for i in range(len(files)):
                    path = root + '\\' + files[i]
                    pathlist.append(path)

        value = []
        for path in pathlist:
            nc_obj = Dataset(path)
            # 查看nc文件中的变量
            nc_name = []
            for i in nc_obj.variables.keys():
                nc_name.append(i)
            # shp经纬度在nc文件中对应的位置lat_num,lon_num
            latitude = nc_obj.variables['latitude'][:]
            lat_num = np.where(abs(latitude - lat) <= 0.25 / 2)
            longitude = nc_obj.variables['longitude'][:]
            lon_num = np.where(abs(longitude - lon) <= 0.25 / 2)
            
            # 当下载气象数据为5变量时
            if len(nc_name) == 5:
                EPOCH = datetime(1900, 1, 1)  # 1900/1/1
                # 查看每个变量的信息，nc_name[4]变量名称
                for i in range(len(nc_obj.variables[nc_name[4]])):
                    time = EPOCH + timedelta(days=int(nc_obj.variables[nc_name[3]][:][i])/24)
                    month = str(time).split("-")[1]  # 获取当前月份
                    # 判断并获取当前气象数据的数组
                    if str(nc_obj.variables[nc_name[4]][i,0,0,0]).replace(".", "").isdigit():
                        nc = nc_obj.variables[nc_name[4]][i, 0, :]
                    else:
                        nc = nc_obj.variables[nc_name[4]][i, 1, :]
                    # 为气象数据赋值和变量名
                    if nc_name[4] == 'tp':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] * 1000), 'Pre' + month])
                    elif nc_name[4] == 'ssr':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] / 1000000), 'Rad' + month])
                    elif nc_name[4] == 't2m':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] - 273.15), 'Tem' + month])
            
            # 当下载气象数据为4变量时
            if len(nc_name) == 4:
                EPOCH = datetime(1900, 1, 1)  # 1900/1/1
                # 查看每个变量的信息
                for i in range(len(nc_obj.variables[nc_name[3]])):
                    time = EPOCH + timedelta(days=int(nc_obj.variables[nc_name[2]][:][i])/24)
                    month = str(time).split("-")[1]  # 获取当前月份
                    # 获取当前气象数据的数组
                    nc = nc_obj.variables[nc_name[3]][i, :]
                    # 为气象数据赋值和变量名
                    if nc_name[3] == 'tp':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] * 1000), 'Pre' + month])
                    elif nc_name[3] == 'ssr':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] / 1000000), 'Rad' + month])
                    elif nc_name[3] == 't2m':
                        value.append([round(nc[lat_num[0][0]][lon_num[0][0]] - 273.15), 'Tem' + month])

        return value

    def getmaxEVI(self):
        '''
        获取下载S2数据中 农场最大EVI对应的影像
        '''
        pathlist = []
        for root, dirs, files in os.walk(self.eviimg):
            if len(files) != 0:
                for i in range(len(files)):
                    if 'tif' in str(files[i]).split('.')[-1] and 'EVI' in str(files[i]):
                        path = root + '\\' + files[i]
                        pathlist.append(path)
        shape = ogr.Open(self.shp_path)
        lyr = shape.GetLayer()
        count = lyr.GetFeatureCount()

        EVI_all = []
        for path in pathlist:
            srcImage = gdal.Open(path)
            geoTrans = srcImage.GetGeoTransform()
            geoProj = srcImage.GetProjection()
            EPSG = int(str(geoProj).split('"')[-2])
            mean = 0
            for i in range(count):
                poly = lyr.GetFeature(i)
                geometry = poly.GetGeometryRef()
                points = geometry.ExportToWkt().replace('POLYGON ((', '').replace('))', '').split(",")
                minX, minY, maxX, maxY = self.e.shpExtent(srcImage, points, EPSG) # shp边界
                ulX, ulY = self.e.world2Pixel(geoTrans, minX, maxY) # shp起点坐标
                lrX, lrY = self.e.world2Pixel(geoTrans, maxX, minY) # shp终点坐标
                pxWidth = abs(int(lrX - ulX))
                pxHeight = abs(int(lrY - ulY))
                eviclip = srcImage.ReadAsArray(ulX, ulY, pxWidth, pxHeight) # 读取shp范围的栅格数据
                mean = mean + eviclip.mean()
            EVI_all.append(mean/count) # EVI均值 mean/count
        EVI_max = pathlist[EVI_all.index(max(EVI_all))] # 最大EVI影像路径
        
        eviout = self.out_path + '\\' + 'EVI2.tif'
        clip = Cliptif(self.shp_path,EVI_max,eviout)
        clip.main()
        return eviout

    def create(self):
        '''
        创建10波段影像 EVI+9气象
        '''
        maxEVI = self.getmaxEVI()
        dataset = gdal.Open(maxEVI)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        value = self.getvalue()
        pathlist = [maxEVI] # 创建列表，按波段顺序添加
        for i in range(len(value)):
            im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
            im_data.ravel()[np.flatnonzero(im_data)] = value[i][0]  # 将栅格数组中非0部分赋气象数据值
            outpath = self.out_path + '\\' + str(value[i][1]) +'.tif'
            self.e.write_img(outpath, im_proj, im_geotrans, im_data) # 创建气象数据的tif文件
            pathlist.append(outpath)
            outpath = None
        paths = sorted(pathlist)
        print(paths)
        
        # 创建EVI2_ERA5数据
        driver = gdal.GetDriverByName("GTiff")
        EVI2_ERA5 = self.out_path + '\\' + 'EVI2_ERA5.tif'
        output = driver.Create(EVI2_ERA5, im_width, im_height, 10, gdal.GDT_Float32)
        output.SetGeoTransform(im_geotrans)
        output.SetProjection(im_proj)
        # 按列表顺序为栅格数据添加波段
        for i in range(len(paths)):
            dataset = gdal.Open(paths[i])
            im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
            output.GetRasterBand(i+1).WriteArray(im_data)
        output.FlushCache()
        output = None

        return paths,EVI2_ERA5

class tifToPng():
    def __init__(self, rasterpath, save_path) :
        self.rasterpath = rasterpath  # 估产tif
        self.save_path = save_path  # 保存路径
        
    def getTif(self):
        dataSet = gdal.Open(self.rasterpath)
        if dataSet == None:
           print('数据导入有误')
        im_width = dataSet.RasterXSize  # 列数
        im_height = dataSet.RasterYSize  # 行数
        im_data = dataSet.ReadAsArray(0, 0, im_width, im_height)
        nodata = dataSet.GetRasterBand(1).GetNoDataValue()
        # 设置背景值
        if nodata is None:
            nodata = 0
        else:
            nodata = nodata
        im_data[im_data == nodata] = np.nan
        
        # 均值 标准差
        mean = np.nanmean(im_data)
        std = np.nanstd(im_data)
        min = np.nanmin(im_data)
        max = np.nanmax(im_data)
        
        min_num = np.floor((mean - min) / (1/2*std)) #向下取整
        max_num = np.floor((max - mean) / (1/2*std))
        reclass_list = [round(min, 0), round(max, 0)]
 
        min_num = 1 if min_num <= 3 else 2  # 最小值与均值之间允许存在的刻度个数
        max_num = 1 if max_num <= 3 else 2
        
        # 按照标准差计算刻度值
        for i in range(0, min_num):
            reclass_list.append(round(mean - std * (i*2 + 1)/2, 0))
        for j in range(0, max_num):
            reclass_list.append(round(mean + std * (j*2 + 1)/2, 0))
        
        # （mean-min）/std小于3/2时，从另一方多算一个5/2倍标准差
        if min_num == 1 and (max-mean)/std >=5/2:
            reclass_list.append(round(mean + std * 5/2, 0))
        if max_num == 1 and (mean-min)/std >=5/2:
            reclass_list.append(round(mean - std * 5/2, 0))
        
        if len(reclass_list) < 6:
            reclass_list = [round(min, 0), round(max, 0)]
            mean = (max-min)/5
            for i in range(0,4):
                reclass_list.append(round(min + (i+1)*mean, 0))
        list = sorted(reclass_list)
        # print(list)

        # 重分类
        data_class = np.full(shape=im_data.shape, fill_value=np.nan)
        data_class[im_data <= list[1]] = 0
        data_class[(im_data > list[1]) & (im_data <= list[2])] = 1
        data_class[(im_data > list[2]) & (im_data <= list[3])] = 2
        data_class[(im_data > list[3]) & (im_data <= list[4])] = 3
        data_class[im_data >= list[4]] = 4

        return data_class, list
    
    def drawPng(self):
        im_data, list = self.getTif()
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置全局中文字体为黑体
        # 自适应调整子图
        plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
        # 去除刻度
        plt.xticks([])
        plt.yticks([])
        # 去除边框
        plt.axis('off')

        plt.imshow(im_data, cmap='RdYlGn')
        cmp = mpl.cm.RdYlGn
        norm = mpl.colors.BoundaryNorm(list, cmp.N)
        # plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmp), label='产量(kg/亩)')

        # 保存图片
        path = self.save_path + "\PreYield_clip.png"
        plt.savefig(path)
        # plt.show()

        plt.close()
        plt.figure(figsize=(0.5,1.1),dpi=300,facecolor='black')
        # 去除边框
        plt.axis('off')
        plt.rcParams['font.size'] = 4
        clb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmp),fraction=0.5)
        clb.ax.set_title('产量(kg/亩)',fontsize = 5,pad = 3,color = 'w')#设置标题
        clb.ax.yaxis.set_ticks_position('left')#刻度放在左侧
        clb.ax.tick_params(color = 'w',length = 2,width = 0.5,labelcolor = 'w')#调整刻度样式
        path = self.save_path + "\PreYield_clip_colorbar.png"
        plt.savefig(path)

class HLM_corn():
    def __init__(self, rasterpath, shppath, savepath, coef_excel) :
        self.rasterpath = rasterpath  # 估产tif
        self.shppath = shppath  # 田块边界shp
        self.savepath = savepath  # 产量tif
        self.coef_excel = coef_excel # 估产系数表格
        self.e = Cliptif(self.shppath, self.rasterpath, self.savepath)

    def getyield(self,params):
        '''
        估产tif计算
        '''
        dataset = gdal.Open(self.rasterpath)
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        # band0:EVImax; band1:Pre03; band4:Rad03; band7:Tem03
        for i in range(0, len(im_data)):
            im_data[im_data < 0.0] = np.nan
            
        c = params[9]
        d = params[19]
        for i in range(0, 9):
            c = c + params[i] * im_data[i + 1]
            d = d + params[i + 10] * im_data[i + 1]
        Yhlm = c * im_data[0] + d
        savepath = self.savepath + '\\' + 'PreYield.tif'
        self.e.write_img(savepath, im_proj, im_geotrans, Yhlm)

        return savepath
    
    def yield_adjust(self,yield_mean,params):
        '''
        最终估产tif和png
        '''
        # 获取估产系数
        df = pd.read_excel(self.coef_excel)
        name  = str(self.shppath).split("\\")[-1].split(".")[0]
        row_id = df.loc[df.iloc[:,0] == name].index.values
        coef_list=[[df.iloc[row,1],df.iloc[row,2:].values.cumprod()[-1]/yield_mean] for row in row_id] # 地块id 最终系数
        
        driver = ogr.GetDriverByName("ESRI Shapefile")
        shape = driver.Open(self.shppath,1)
        lyr = shape.GetLayer()
        # shp添加估产系数字段
        feature_defn = lyr.GetLayerDefn()
        field_name = 'coef'
        fieldIndex = feature_defn.GetFieldIndex(field_name)
        if fieldIndex < 0:
            fieldDefn  = ogr.FieldDefn(field_name, ogr.OFTReal)
            fieldDefn.SetPrecision(5)
            lyr.CreateField(fieldDefn,1)
        
        # 为字段赋值
        for feature in lyr:
            feature_value = int(feature.GetField('id'))
            if feature_value in [i[0] for i in coef_list]:
                coef_value = [i[1] for i in coef_list if i[0]==feature_value][0]
            else:
                coef_value = float(1)
            feature.SetField(field_name,coef_value)
            lyr.SetFeature(feature)

        # 初步估产结果tif 作为模板
        img = self.getyield(params)
        dataset = gdal.Open(img)
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        geoTrans = dataset.GetGeoTransform()
        geoProj = dataset.GetProjection()
        
        # 根据模板tif属性信息创建对应标准的目标栅格
        coef_tif = self.savepath + '\\' + name + "_coef.tif"
        driver = gdal.GetDriverByName("GTiff")
        target_ds = driver.Create(coef_tif, im_width, im_height, 1, gdal.GDT_Float32)
        target_ds.SetProjection(geoProj)
        target_ds.SetGeoTransform(geoTrans)
        
        # 设置背景数值
        band = target_ds.GetRasterBand(1)
        NoData_value = 0
        band.SetNoDataValue(NoData_value)
        band.FlushCache()
        # 栅格化函数
        gdal.RasterizeLayer(target_ds, [1], lyr, options=["ATTRIBUTE=coef"])
        coef_data = band.ReadAsArray()
        
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_data[im_data == params[19]] = np.nan # 为非目标区域栅格 赋空值
        yield2 = np.multiply(coef_data,im_data) # 最终估产值
        yield2[yield2 == 0.] = np.nan
        savepath = self.savepath + '\\' + name + '_PreYield.tif' # 产量tif文件
        self.e.write_img(savepath, geoProj, geoTrans, yield2)
        # del dataset
        
        # tif to png
        png = tifToPng(savepath,self.savepath)
        png.drawPng()    
      
class HLM_wheat():
    def __init__(self, rasterpath, shppath, savepath) :
        self.rasterpath = rasterpath  # 估产tif
        self.shppath = shppath  # 估产tif
        self.savepath = savepath  # 产量tif
        self.e = Cliptif(self.shppath, self.rasterpath, self.savepath)

    def getyield(self,params):
        '''
        估产tif计算
        '''
        dataset = gdal.Open(self.rasterpath)
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        # band0:EVImax; band1:Pre03; band4:Rad03; band7:Tem03
        for i in range(0, len(im_data)):
            im_data[im_data < 0.0] = np.nan
            
        c = params[9]
        d = params[19]
        for i in range(0, 9):
            c = c + params[i] * im_data[i + 1]
            d = d + params[i + 10] * im_data[i + 1]
        Yhlm = c * im_data[0] + d
        savepath = self.savepath + '\\' + 'PreYield.tif'
        self.e.write_img(savepath, im_proj, im_geotrans, Yhlm)
        clippath = self.savepath + '\\' + 'PreYield_clip.tif'
        clip = Cliptif(self.shppath, savepath, clippath)
        clip.main()
        
        # tif to png
        png = tifToPng(clippath,self.savepath)
        png.drawPng()

if __name__ == "__main__":
    nc_path = r"E:\ncdata"  # 只包含三个nc文件的文件夹
    evi_path = r"E:\潍柴雷沃估产数据\s22"  # sentinel下载的路径

    shp_path = r"E:\1潍柴雷沃玉米估产\2022年估产\SHP\masi.shp"  # 农场边界数据
    save_path = r"E:\1潍柴雷沃玉米估产"  # 产量tif
    crop = 'corn' # corn or wheat

    stack = Stack(nc_path, shp_path, evi_path, save_path)
    tifpath = stack.create() # EVI_ERA5路径
    for i in tifpath[0]:
        os.remove(i)
    
    ##### 玉米估产系数 ########image.png
    yield_mean = 653.7635
    model_corn = np.array([-4.04887897, -0.45785285, -0.31328856,  1.70103846,  1.86853625,  0.23360395,
                       2.2995999,   1.25867253,  1.62653312,  1.77348436, -3.97690748,  0.35837612,
                       0.79250955,  4.79144255,  4.62676819,  2.0039969,   7.09561622,  4.14458283,
                       4.39637883,  1.9991423 ])
    
    ###### 小麦估产系数 ########
    model_wheat = np.array([2.31604279, 1.28980504, 1.74452386, 2.4879997, 2.59906958, 3.52681908,
                            1.87340677, 3.08237613, 3.74004541, 1.14715985, 2.62341881, 1.34778218,
                            2.4000839, 3.28883175, 3.07335727, 5.34414233, 2.92456618, 3.85071211,
                            3.52495191, 1.36363804])
    
    if crop == 'corn':
        coef_excel = r"E:\1潍柴雷沃玉米估产\农场调查.xlsx"  #农场调查
        hlm = HLM_corn(tifpath[1],shp_path,save_path,coef_excel)
        hlm.yield_adjust(yield_mean,model_corn)

    if crop == 'wheat':
        hlm = HLM_wheat(tifpath[1],shp_path,save_path)
        hlm.getyield(model_wheat)