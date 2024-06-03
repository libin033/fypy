=================================
ahi调用示例
=================================

H8/H9 AHI HSD文件解析和块拼接成FULL DISK
-----------------------------------------

.. code-block:: python

    from fypy.ahi.ahiscene import AHIScene

    nowdate = datetime.datetime.strptime('20230115_0400', '%Y%m%d_%H%M')
    srcdir = r'./H8'
    dstdir = r'./H8/20230115'
    outname = os.path.join(dstdir, 'HS_H08_%s_FLDK.h5' %(nowdate.strftime('%Y%m%d_%H%M')))
    SatID = 'H09'
    overwrite = 1
    if not os.path.isdir(dstdir) :
        os.makedirs(dstdir)

    tempdir = os.path.join(dstdir, 'temp')
    if not os.path.isdir(tempdir) :
        os.makedirs(tempdir)

    for chid in np.arange(1, 4) :
        searchfile = os.path.join(srcdir,
                                  r'*HS_{satid}_{date}_B{chid}_FLDK_R*_S*.DAT.bz2'.format(
                                      satid=SatID,
                                      date=nowdate.strftime('%Y%m%d_%H%M'),
                                      chid='%02d' %(chid)))
        filelist = glob.glob(searchfile)
        if len(filelist) == 0:
            print('未匹配到文件【%s】' %(searchfile))
            continue

        scene = AHIScene(subpoint=140.7, resolution=0.01)
        ds = scene.hsdBlock(filelist, tempdir)

        ds = scene.clip(ds)

        scene.ds2tiff(r'./H8/20230115/test.tif', srcDS=ds)

        scene.ds2netcdf(r'./H8/20230115/test.nc', 'data', srcDS=ds)
        scene.ds2hdf(r'./H8/20230115/test.hdf', 'data', srcDS=ds)

H8/H9 AHI 数据经纬度、行列号互转
-----------------------------------------

.. code-block:: python

    from fypy.ahi.ahiscene import AHIScene
    lat = [30]
    lon = [104]

    scene = AHIScene(subpoint=140.7, resolution=0.01)
    # lat, lon转xy
    x, y = scene.latlon2xy(lat, lon)
    print(x, y)
    # xy转lon,lat
    lon, lat = scene.xy2latlon(x, y)

    print(lon, lat)

H8/H9 AHI 各类角度计算
-----------------------------------------

.. code-block:: python

    from fypy.ahi import AHIScene

    nowdate = datetime.datetime.strptime('20230115_0400', '%Y%m%d_%H%M')
    scene = AHIScene(subpoint=140.7, resolution=0.01)

    lat = np.array([20])
    lon = np.array([104.7])

    suna, sunz = scene.calSolarAngle(lat, lon, nowdate)
    sata, satz = scene.calSatAngle(lat, lon)

    # 相对方位角
    rela = scene.calRelativeAzimuth(sata, suna)

    # 耀斑角
    sungl = scene.calSunGL(satz, sunz, rela)


H8/H9 AHI 真彩图绘制
-----------------------------------------

.. code-block:: python

    from fypy.ahi import AHIScene