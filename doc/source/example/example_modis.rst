 =================================
modis调用示例
=================================

modis L1定标，L1/L2轨道数据投影
-----------------------------------------

.. code-block:: python

    from fypy.modis import MODISScene
    from fypy.tools import readhdf4
    scene = MODISScene()
    filename = r'MOD021KM.A2022001.0330.061.2022001133325.hdf'
    data = scene.calibrate(filename, 'EV_500_Aggr1km_RefSB')

    geofile = r'MOD03.A2022001.0330.061.2022001090339.hdf'
    lat = readhdf4(geofile, 'Latitude')
    lon = readhdf4(geofile, 'Longitude')

    ds = scene.project(data, lat, lon,
                       resolution=0.01, resampleAlg='near')
    scene.ds2tiff('test.tif', ds)

    ds = scene.quickProject(filename, 'EV_500_Aggr1km_RefSB')
    scene.ds2tiff(r'test4.tif', ds)

modis L3级数据产品拼接、投影、裁剪
-----------------------------------------

.. code-block:: python

    from fypy.modis import MODISScene

    filename = glob.glob(r'MOD13A2.A2021353.*.hdf')
    sdsname = '1 km 16 days NDVI'

    scene = MODISScene()
    ds = scene.converModisByGDAL(filename, sdsname, resolution=None)
    scene.ds2tiff(r'test.tif', ds)

