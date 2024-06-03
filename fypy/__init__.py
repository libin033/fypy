# -*- coding:utf-8 -*-
'''
@Project  : fypy
@File     : __init__.py
@Modify Time      @Author    @Version    
--------------    -------    --------    
2022/7/20 11:29      Lee       1.0         
@Description
------------------------------------
 
'''

__import__ ('pkg_resources').declare_namespace (__name__)

from . import version
from .version import __version__

# from .version import get_versions
# __version__ = get_versions()['version']
# del get_versions

from fypy.ahi import AHIScene
from fypy.fy3 import FY3Scene
from fypy.fy4 import FY4Scene
from fypy.modis import MODISScene

