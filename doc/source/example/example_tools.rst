=================================
tools调用示例
=================================


nc2tif
-----------------------------------------

.. code-block:: python

        from fypy.tools import nc2tif

        filename = r'./L3m_20210607__171505720_4_AV-OLB_ZSD_DAY_00.nc'
        nc2tif('./data/test.tif', filename, 'ZSD_mean')

gifpro
-----------------------------------------

.. code-block:: python


        from fypy.tools import creategif1, splitgif

        pathin = r'./images'

        fils = glob.glob(os.path.join(pathin, '*.PNG'))
        fils.sort()
        creategif1('./data/test.gif', filelist=fils, duration=1 / 2.0)

        splitgif('./data/', './data/test.gif')
