=================================
fy3调用示例
=================================

FY3 MERSI L1数据辐射定标
-----------------------------------------
支持对FY3 MERSI 1KM、250M L1数据辐射定标

.. code-block:: python

    from fypy.fy3 import FY3Scene
    from fypy.tools import readhdf

    l1file  = 'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = 'FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'

    scene = FY3Scene(l1file)

    # 1-4波段
    band14 = scene.Calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')

    # 5 - 19波段
    band519 = scene.Calibration(l1file, '/Data/EV_1KM_RefSB')

    # 20 - 23
    band2023 = scene.Calibration(l1file, '/Data/EV_1KM_Emissive')

    # 24-25
    band2425 = scene.Calibration(l1file, '/Data/EV_250_Aggr.1KM_Emissive')

FY3 轨道数据产品投影
-----------------------------------------
支持对FY3 MERSI L1、L2等轨道产品（ORBT）数据进行投影（WGS84）

.. code-block:: python

    from fypy.fy3 import FY3Scene
    from fypy.tools import readhdf

    l1file  = 'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = 'FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'
    scene = FY3Scene(l1file)
    data = scene.calibrate(l1file, '/Data/EV_250_Aggr.1KM_RefSB')

    lat = readhdf(geofile, '/Geolocation/Latitude')
    lon = readhdf(geofile, '/Geolocation/Longitude')

    ds = scene.project(data, lat, lon,
                       resolution=None, vmin=0, vmax=10000, resampleAlg='near')

    scene.ds2tiff('test.tif', ds)


FY3 10度块产品拼接（支持GLL、HAM投影）
-----------------------------------------
FY3的10度块产品大部分是GLL（等经纬WGS84），由于NDVI产品是HAM投影，新增了对HAM投影

.. code-block:: python

    from fypy.fy3 import FY3Scene

    filelist  = glob.glob(r'*NDVI*.HDF')

    scene = FY3Scene()

    ds = scene.block10Merge(filelist, ProdID='NDVI', sdsname='1000M_10day_NDVI', fillvalue=-32768)
    scene.ds2tiff('test1.tif', ds)

    # 将Harmmer投影转换为WGS84
    harmds = scene.hammer2Wgs84(ds, xRes=0.01, yRes=0.01)
    scene.ds2tiff('test2.tif', harmds)

FY3 绘制真彩图
-----------------------------------------

.. code-block:: python

    from fypy.fy3 import FY3Scene

    l1file = r'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    scene = FY3Scene(filename=l1file)
    scene.load(l1file, ProdID='truecolor')
    # scene.show()
    scene.saveThematic(r'test.png')

