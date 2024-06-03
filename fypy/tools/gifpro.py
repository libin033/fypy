# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : gifpro.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
对GIF文件操作：

-    实现对图像合成GIF文件

-    实现将GIF文件拆分为单幅影像

'''
import imageio
import os
from PIL import Image, ImageDraw


def creategif(outname, filelist, duration=0.1):
    '''
    使用 imageio 生成 GIF

    Parameters
    ----------
    outname: str
        输出的GIF文件名
    filelist: list
        输入需要创建GIF的图片列表
    duration: float
        时间间隔，控制GIF播放速度
    Returns
    -------

    '''

    frames = []

    if not outname.endswith(".gif") and not outname.endswith(".GIF"):
        outname  += ".gif"

    num = len(filelist)
    for item in filelist:
        num -= 1
        if not item.lower().endswith(".jpg") and not item.lower().endswith(".png"):
            print('文件后缀不是PNG或JPG格式文件，将跳过该文件【%s】' %(os.path.basename(item)))
            continue

        if not os.path.isfile(item):
            print("%s not exist, will be continue !!!!" % item)
            continue

        print("【%d】%s" % (num,item))
        frames.append(imageio.imread(item))

    imageio.mimsave(outname, frames, "GIF", duration = duration)
    # print("outname:",outname)

def creategif1(outname, filelist, duration=0.1):
    '''
    使用 Image 生成 GIF

    Parameters
    ----------
    outname: str
        输出GIF文件名
    filelist: list
        输入需要创建GIF的图片列表
    duration: float
        时间间隔，控制GIF播放速度
    Returns
    -------

    '''

    imgs = []
    for item in filelist:
        print(item)
        if not item.lower().endswith(".jpg") and not item.lower().endswith(".png"):
            print(not item.lower().endswith(".jpg") , not item.lower().endswith(".png"))
            continue

        if not os.path.isfile(item):
            print("%s not exist, will be continue !!!!" % item)
            continue

        temp = Image.open(item)
        imgs.append(temp)

    imgs[0].save(outname, save_all=True, append_images=imgs, duration=duration)
    print("outname:",outname)

    return outname

def splitgif(outdir, gifname, format='PNG', quality=95):
    '''
    拆解GIF文件，还原成单幅图像

    Parameters
    ----------
    outdir : str
        输出路径
    gifname: str
        输入gif文件名
    format: str
        输出文件的文件格式
    quality: int
         保存图像的质量，值的范围从1（最差）到95（最佳）。
         默认值为95，使用中应尽量避免高于95的值;
         100会禁用部分JPEG压缩算法，
         并导致大文件图像质量几乎没有任何增益

    Returns
    -------
        list
        拆解的文件列表
    '''

    img = Image.open(gifname)

    filelist = []
    try:
        for i in range(img.n_frames) :
            img.seek(i)

            if format in ['PNG','png'] :
                new_frame = Image.new('RGBA', img.size)
            else:
                new_frame = Image.new('RGB', img.size)
            new_frame.paste(img)

            outname = os.path.join(outdir, os.path.basename(gifname)+'_%d.%s' %(i, format))
            new_frame.save(outname, format=format, quality=quality)
            print('[%d]成功创建：%s' %(i, outname))

            filelist.append(outname)
    except EOFError:
        pass

    return filelist










