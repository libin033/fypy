#coding:utf-8

# step1:
#     python setup.py sdist bdist_wheel
#
# step2:
#     C:\Users\admin\AppData\Roaming\Python\Python38\Scripts\twine.exe upload dist/*
# or
#     python -m twine upload dist/*

import os
from setuptools import setup,find_packages,Extension

DIR = os.path.dirname(os.path.abspath(__file__))
# INSTALL_PACKAGES = open(os.path.join(DIR, 'requirements.txt')).read().splitlines()


# __version__ = '1.0.0'
# from fypy import __version__
# from fypy.version import version as __version__

from fypy.version import get_versions
__version__ = get_versions()['version']

mod_expy=Extension('',
                   sources=[],
                   libraries=[],
                   include_dirs=[],
                   library_dirs=[],
                   define_macros=[('MAJOR_VERSION','1'),('MINOR_VERSION','0')],
                   )

#区分各种平台的编译
# if platform.system()=='Darwin':
#     os.environ['ARCHFLAGS']='-arch x86_64'
#     extmodlist=[mod_expy,]
# elif platform.system()=='Linux':
#     os.environ['ARCHFLAGS']='-arch i386 -arch x86_64'
#     extmodlist=[mod_expy,]
# else:
#     raise RuntimeError('Unknown system()=%s'%repr(platform.system()))


name = 'fypy'
description = "fypy Library for FENGYUN Satellite "
readme = open('README.md', encoding="utf-8").read()
content_type='text/markdown'

setup(
    name = name, # 项目名
    version = __version__, # 如 0.0.1/0.0.1.dev1
    description = description,
    long_description = readme,
    long_description_content_type = content_type,
    # url='', # 如果有github之类的相关链接
    author = 'The fypy Team', # 作者
    # author_email='xxx@163.com', # 邮箱
    license = 'MIT',
    platforms = ["windows", 'linux'],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    url='https://github.com/libin033/fypy',
    keywords=['FengYun', 'Meteorology', 'RemoteSensing'], # 关键词之间空格相隔一般
    # 需要安装的依赖
    namespace_packages = ['fypy'],
    install_requires = [
        'numpy >= 1.2.0',
        'h5py >= 1.0.0',
        'netcdf4 >= 1.0.0',
        'tqdm >= 4.0.0',
        'gdal >= 2.0.0',
        'pyshp >= 2.1.0',
        'imageio >= 2.15.0',
    ],
    setup_requires = [ ],
    packages=find_packages(),
    package_data={
    },
    include_package_data = True, #
    entry_points={ }, #如果发布的库包括了命令行工具
    # ext_modules=[mod_expy,],
    package_dir={'fypy': 'fypy'},
    python_requires='>=3.7',
)