=================================
draw调用示例
=================================

绘制专题图（drawThematic）
-----------------------------------------

.. code-block:: python

    from fypy.draw import drawThematic

    from fypy.tools import readnc, readtiff
    filename = r'./data/test.tif'

    data, trans, prj = readtiff(filename)
    mpro = drawThematic(data, prodid='RSM', cbfile='./colorbar_RSM.txt')

    # 增加Text
    mpro.text(0.1, 0.1, 'RSM')

    # 保存为图片文件
    mpro.save('./data/test1.png')



绘制图例（colorbar）
-----------------------------------------

.. code-block:: python

    from fypy.draw import colorbar, getColorList

    cbfile = r'./colorbar_CLM.txt'
    cblist, cbtitle, cbdir = getColorList(cbfile)
    outname = os.path.join('./data/cb/', os.path.basename(cbfile).replace('.txt', '_v.png'))
    colorbar( outname, cblist, width=800, height=600, cbtitle=cbtitle, ticks=None, ticklabels=None,
              vmin=None, vmax=None,  orientation='v', fmt='%.0f', font=None,
              tickfontsize=45, titlesize=50, loc='center', split=False, encoding='gbk')

    outname = os.path.join('./data/cb/', os.path.basename(cbfile).replace('.txt', '_h.png'))
    colorbar( outname, cblist, width=1260, height=180, cbtitle=cbtitle, ticks=None, ticklabels=None,
              vmin=None, vmax=None,  orientation='h', fmt='%.0f', font=None,
              tickfontsize=45, titlesize=50, loc='center', split=False)

