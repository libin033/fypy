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

from .version import get_versions
__version__ = get_versions()['version']
del get_versions

# __version__ = "1.0.0"
# version_tuple = (1, 0, 0)
# try:
#     from fypy.version import version as __version__  # noqa
# except ModuleNotFoundError:
#     raise ModuleNotFoundError(
#         "No module named fypy.version. This could mean "
#         "you didn't install 'fypy' properly. Try reinstalling ('pip "
#         "install').")


# try:
#     __import__('pkg_resources').declare_namespace(__name__)
# except ImportError:
#     pass  # must not have setuptools
