=================================
rsdp调用示例
=================================

FY3 MERSI大气校正
-----------------------------------------
对FY3 MERSI L1数据的可见光通道进行大气校正

.. code-block:: python

    import datetime
    from fypy.RSDP.AtmCorr_FY3D_MERSI import AtmCorr_FY3D_MERSI
    from fypy.tools import readhdf, readhdf_fileinfo, writehdf

    l1file = r'FY3D_MERSI_GBAL_L1_20220624_0620_1000M_MS.HDF'
    geofile = r'FY3D_MERSI_GBAL_L1_20220624_0620_GEO1K_MS.HDF'

    corrpro = AtmCorr_FY3D_MERSI()
    nowdate = datetime.datetime.strptime('20220624_0620', '%Y%m%d_%H%M')
    ref = corrpro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')
    dn = readhdf(l1file, '/Data/EV_250_Aggr.1KM_RefSB')
    lat = readhdf(geofile, '/Geolocation/Latitude')
    lon = readhdf(geofile, '/Geolocation/Longitude')
    satA = readhdf(geofile, '/Geolocation/SensorAzimuth')
    satZ = readhdf(geofile, '/Geolocation/SensorZenith')
    solA = readhdf(geofile, '/Geolocation/SolarAzimuth')
    solZ = readhdf(geofile, '/Geolocation/SolarZenith')
    overwrite = 1
    for bandid in range(1, 5) :
        corrband = corrpro.FLAASH(nowdate, ref[bandid-1]*0.01, lat, lon,
                                  sunz=21.0, suna=278.0, satz=0.0, sata=60.0,
                                  SatID='FY3D', InstID='MERSI', BandId=bandid)

        writehdf('test.hdf', 'B%02d' %(bandid), corrband, overwrite=overwrite)
        overwrite = 0

GF1、2、4、6大气校正
-----------------------------------------
对GF1、2、4、6 L1数据的可见光通道进行大气校正

.. code-block:: python

    import os
    import datetime
    from fypy.RSDP.AtmCorr_GF import AtmCorr_GF

    srctif = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-MSS1.tiff'
    srcrpb = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-MSS1.rpb'
    srcxml = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-MSS1.xml'

    # srctif = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-PAN1.tiff'
    # srcrpb = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-PAN1.rpb'
    # srcxml = r'./GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163/GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-PAN1.xml'

    basename = os.path.basename(srctif)
    namelist = basename.split('_')

    outdir = r'./'
    atmcorname = os.path.join(outdir, basename)
    orthoname = os.path.join(outdir, basename.replace('.tiff', '_ortho.tiff'))
    fillvalue = 65535

    mpro = AtmCorr_GF()

    SatID = namelist[0]
    InstID = namelist[1]
    nowdate = datetime.datetime.strptime(namelist[4], '%Y%m%d')

    # 大气校正
    mpro.FLAASH(atmcorname, srctif, srcxml, nowdate, SatID, InstID, InstType='MSS') # ImageType: PAN、MSS

    # 正射校正
    mpro.ortho_rectification(orthoname, atmcorname, srcrpb)

     # Note：此处需要先对MSS和PAN做大气校正，PAN的大气校正仅仅只做了辐射定标，未真正做大气校正
    # 多光谱波段
    mssfile = os.path.join(outdir, 'GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-MSS1_ortho.tiff')
    # 全色波段
    panfile = os.path.join(outdir, 'GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-PAN1_ortho.tiff')

    pansharpenname = os.path.join(outdir, 'GF1_PMS1_E91.1_N29.4_20180125_L1A0002958163-PAN1_pansharpen.tiff')
    dict_spe = {
        mssfile : [1, 2, 3, 4]
    }

    # 波段融合
    mpro.pansharpen(pansharpenname, panfile, dict_spectrac=dict_spe,
                    nodata_value=fillvalue, creation_options=["COMPRESS=LZW"])

Landsat大气校正
-----------------------------------------
对Landsat L1数据的可见光通道进行大气校正

.. code-block:: python

    import os
    import datetime
    from fypy.RSDP import AtmCorr_Landsat

    srcdir = r'./LC08_L1TP_129039_20221022_20221101_02_T1'
    outdir = r'./LC08_L1TP_129039_20221022_20221101_02_T1\test'

    metafile = os.path.join(srcdir, 'LC08_L1TP_129039_20221022_20221101_02_T1_MTL.txt')
    srcfile = os.path.join(srcdir, 'LC08_L1TP_129039_20221022_20221101_02_T1_B1.TIF')
    nowdate = datetime.datetime.strptime('20221022', '%Y%m%d')
    mpro = AtmCorr_Landsat(metafile)
    mpro.FLAASH(nowdate, srcfile=srcfile)


Sentinel大气校正
-----------------------------------------
对Sentinel L1数据的可见光通道进行大气校正

.. code-block:: python

    import os
    import datetime
    from fypy.RSDP import AtmCorr_Sentinel

    srcdir = r'./S2B_MSIL2A_20230411T033539_N0509_R061_T48RVS_20230411T065637'
    outdir = r'./S2B_MSIL2A_20230411T033539_N0509_R061_T48RVS_20230411T065637/test'

    srcfile = os.path.join(srcdir, 'IMG_DATA', 'R10m', 'T48RVS_20230411T033539_B02_10m.jp2')
    metafile = os.path.join(srcdir, 'MTD_TL.xml')

    outname = os.path.join(outdir, 'T48RVS_20230411T033539_B02_10m.tif')

    mpro = AtmCorr_Sentinel()
    metadata = mpro.GetMeta(metafile)
    mpro.FLAASH(outname, srcfile, metadata, SatID='S2B', InstID='MSI', BandId=2)


