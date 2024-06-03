# fypy  
```angular2html
安装命令：
python -m pip install fypy
python -m pip install --upgrade fypy
```

API 参考 [官方使用说明](https://fypy.readthedocs.io/zh/latest/)


<div align="center"> 
    <p>获取最新更新信息，请关注作者公众号</p>
    <img decoding="async" src="https://github.com/libin033/lb_toolkits/blob/master/lb_toolkits/parm/icon/logo_wxgzh.jpg?raw=true" align="middle"  width="30%">
</div>


## 依赖库

<table>
    <tr>
        <th> 库名 </th>
        <th> 版本 </th>
        <th> 库名 </th>
        <th> 版本 </th>
        <th> 库名 </th>
        <th> 版本 </th>
    </tr>
    <tr>
        <td> numpy </td>
        <td> 1.2.0 </td>
        <td> pyhdf </td>
        <td> 0.10.0 </td>
        <td> h5py </td>
        <td> 1.0.0 </td>
    </tr>
    <tr>
        <td> netcdf4 </td>
        <td> 1.0.0 </td>
        <td> tqdm </td>
        <td> 4.0.0 </td>
        <td> gdal </td>
        <td> 2.0.0 </td>
    </tr>
</table>

------------------------------------------------

## fypy
提供了针对 ``FY3`` 、 ``FY4`` 卫星数据的预处理；
* 预处理包含对极轨卫星FY3的5分钟块、轨道数据的投影、拼接裁剪；
* 对 ``FY3`` 的10度块数据进行投影、拼接、裁剪功能
* 对 ``FY4`` 卫星数据的经纬度转行列号，行列号转经纬度，
  以及对数据进行投影转换、裁剪等功能；



