# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : jsonpro.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
对json\ascii\二进制文件进行读写操作

'''
import sys
import os
import numpy as np
import json



def readjson(jsonname, buffering=None, encoding=None, **kw):
    '''
    读取JSON文件

    Parameters
    ----------
    jsonname: str
        JSON文件名
    kwargs

    Returns
    -------
        dict
        json文件内容以dict形式返回
    '''

    fp = open(jsonname, 'r', encoding=encoding)
    dict_info_in = json.load(fp, **kw)

    fp.close()

    return dict_info_in

def writejson(jsonname, dict_info, indent=4, chinese=False):

    if chinese :
        ensure_ascii = False
    else:
        ensure_ascii = True

    fp = open(jsonname, 'w')
    json.dump(dict_info, fp, indent=indent, ensure_ascii=ensure_ascii)
    fp.close()


def readbinary(filename, shape, dtype=np.float32, offset=0, encoding='utf-8'):
    '''
    读取二进制文件

    Parameters
    ----------
    filename: str
        输入文件名
    shape: list or tuple
        获取的维度大小
    dtype: numpy.dtype
        数据类型
    offset: int
        偏移量
    encoding: str
        编码类型
    Returns
    -------
        numpy.array
    '''

    count = 1
    if os.path.isfile(filename) :
        for i in shape :
            count *= i
        data = np.fromfile(filename, dtype=dtype, count= count,offset=offset)

    return data.reshape(shape)


def writebinary(filename, data, overwrite=1, offset=0, encoding='utf-8'):
    '''
    写入二进制文件

    Parameters
    ----------
    filename: str
        输出文件名
    shape: list or tuple
        获取的维度大小
    overwrite: int
        是否覆盖、新建文件, 非0为新建
    offset: int
        偏移量
    encoding: str
        编码类型
    Returns
    -------
        numpy.array
    '''

    data = np.array(data)

    if overwrite:
        fp = open(filename, 'wb')
    else:
        fp = open(filename, 'ab')

    fp.seek(offset)
    data.tofile(fp)

    fp.close()

def readascii(filename, dtype=float, comments='#', delimiter=None,
                     converters=None, skiprows=0, usecols=None, unpack=False,
                     ndmin=0, encoding='bytes', max_rows=None):
    '''
    Load data from a text file.

    Each row in the text file must have the same number of values.

    Parameters
    ----------
    fname : file, str, or pathlib.Path
        File, filename, or generator to read.  If the filename extension is
        ``.gz`` or ``.bz2``, the file is first decompressed. Note that
        generators should return byte strings.
    dtype : data-type, optional
        Data-type of the resulting array; default: float.  If this is a
        structured data-type, the resulting array will be 1-dimensional, and
        each row will be interpreted as an element of the array.  In this
        case, the number of columns used must match the number of fields in
        the data-type.
    comments : str or sequence of str, optional
        The characters or list of characters used to indicate the start of a
        comment. None implies no comments. For backwards compatibility, byte
        strings will be decoded as 'latin1'. The default is '#'.
    delimiter : str, optional
        The string used to separate values. For backwards compatibility, byte
        strings will be decoded as 'latin1'. The default is whitespace.
    converters : dict, optional
        A dictionary mapping column number to a function that will parse the
        column string into the desired value.  E.g., if column 0 is a date
        string: ``converters = {0: datestr2num}``.  Converters can also be
        used to provide a default value for missing data (but see also
        `genfromtxt`): ``converters = {3: lambda s: float(s.strip() or 0)}``.
        Default: None.
    skiprows : int, optional
        Skip the first `skiprows` lines, including comments; default: 0.
    usecols : int or sequence, optional
        Which columns to read, with 0 being the first. For example,
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns.
        The default, None, results in all columns being read.

        .. versionchanged:: 1.11.0
            When a single column has to be read it is possible to use
            an integer instead of a tuple. E.g ``usecols = 3`` reads the
            fourth column the same way as ``usecols = (3,)`` would.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``.  When used with a
        structured data-type, arrays are returned for each field.
        Default is False.
    ndmin : int, optional
        The returned array will have at least `ndmin` dimensions.
        Otherwise mono-dimensional axes will be squeezed.
        Legal values: 0 (default), 1 or 2.

        .. versionadded:: 1.6.0
    encoding : str, optional
        Encoding used to decode the inputfile. Does not apply to input streams.
        The special value 'bytes' enables backward compatibility workarounds
        that ensures you receive byte arrays as results if possible and passes
        'latin1' encoded strings to converters. Override this value to receive
        unicode arrays and pass strings as input to converters.  If set to None
        the system default is used. The default value is 'bytes'.

        .. versionadded:: 1.14.0
    max_rows : int, optional
        Read `max_rows` lines of content after `skiprows` lines. The default
        is to read all the lines.

        .. versionadded:: 1.16.0
    ${ARRAY_FUNCTION_LIKE}

        .. versionadded:: 1.20.0

    Returns
    -------
    out : ndarray
        Data read from the text file.

    See Also
    --------
    load, fromstring, fromregex
    genfromtxt : Load data with missing values handled as specified.
    scipy.io.loadmat : reads MATLAB data files

    Notes
    -----
    This function aims to be a fast reader for simply formatted files.  The
    `genfromtxt` function provides more sophisticated handling of, e.g.,
    lines with missing values.
    '''
    return  np.loadtxt(filename, dtype=dtype, comments=comments, delimiter=delimiter,
                       converters=converters, skiprows=skiprows, usecols=usecols, unpack=unpack,
                       ndmin=ndmin, encoding=encoding, max_rows=max_rows)

def writeascii(filename, data,  fmt='%.18e', delimiter=' ', newline='\n', header='',
               footer='', comments='# ', encoding=None):
    '''
    Save an array to a text file.

    Parameters
    ----------
    fname : filename or file handle
        If the filename ends in ``.gz``, the file is automatically saved in
        compressed gzip format.  `loadtxt` understands gzipped files
        transparently.
    X : 1D or 2D array_like
        Data to be saved to a text file.
    fmt : str or sequence of strs, optional
        A single format (%10.5f), a sequence of formats, or a
        multi-format string, e.g. 'Iteration %d -- %10.5f', in which
        case `delimiter` is ignored. For complex `X`, the legal options
        for `fmt` are:

        * a single specifier, `fmt='%.4e'`, resulting in numbers formatted
          like `' (%s+%sj)' % (fmt, fmt)`
        * a full string specifying every real and imaginary part, e.g.
          `' %.4e %+.4ej %.4e %+.4ej %.4e %+.4ej'` for 3 columns
        * a list of specifiers, one per column - in this case, the real
          and imaginary part must have separate specifiers,
          e.g. `['%.3e + %.3ej', '(%.15e%+.15ej)']` for 2 columns
    delimiter : str, optional
        String or character separating columns.
    newline : str, optional
        String or character separating lines.

        .. versionadded:: 1.5.0
    header : str, optional
        String that will be written at the beginning of the file.

        .. versionadded:: 1.7.0
    footer : str, optional
        String that will be written at the end of the file.

        .. versionadded:: 1.7.0
    comments : str, optional
        String that will be prepended to the ``header`` and ``footer`` strings,
        to mark them as comments. Default: '# ',  as expected by e.g.
        ``numpy.loadtxt``.

        .. versionadded:: 1.7.0
    encoding : {None, str}, optional
        Encoding used to encode the outputfile. Does not apply to output
        streams. If the encoding is something other than 'bytes' or 'latin1'
        you will not be able to load the file in NumPy versions < 1.14. Default
        is 'latin1'.

        .. versionadded:: 1.14.0


    See Also
    --------
    save : Save an array to a binary file in NumPy ``.npy`` format
    savez : Save several arrays into an uncompressed ``.npz`` archive
    savez_compressed : Save several arrays into a compressed ``.npz`` archive

    Notes
    -----
    Further explanation of the `fmt` parameter
    (``%[flag]width[.precision]specifier``):

    flags:
        ``-`` : left justify

        ``+`` : Forces to precede result with + or -.

        ``0`` : Left pad the number with zeros instead of space (see width).

    width:
        Minimum number of characters to be printed. The value is not truncated
        if it has more characters.

    precision:
        - For integer specifiers (eg. ``d,i,o,x``), the minimum number of
          digits.
        - For ``e, E`` and ``f`` specifiers, the number of digits to print
          after the decimal point.
        - For ``g`` and ``G``, the maximum number of significant digits.
        - For ``s``, the maximum number of characters.

    specifiers:
        ``c`` : character

        ``d`` or ``i`` : signed decimal integer

        ``e`` or ``E`` : scientific notation with ``e`` or ``E``.

        ``f`` : decimal floating point

        ``g,G`` : use the shorter of ``e,E`` or ``f``

        ``o`` : signed octal

        ``s`` : string of characters

        ``u`` : unsigned decimal integer

        ``x,X`` : unsigned hexadecimal integer

    This explanation of ``fmt`` is not complete, for an exhaustive
    specification see [1]_.

    References
    ----------
    .. [1] `Format Specification Mini-Language
           <https://docs.python.org/library/string.html#format-specification-mini-language>`_,
           Python Documentation.
    '''
    data = np.array(data)

    np.savetxt(filename, data, fmt, delimiter, newline, header,
               footer, comments, encoding)

def loadarray(file, mmap_mode=None, allow_pickle=False, fix_imports=True,
              encoding='ASCII'):
    ''' 读取numpy保存输出的.npy文件 '''
    return np.load(file, mmap_mode, allow_pickle, fix_imports, encoding)


def savearray(filename, data, allow_pickle=True, fix_imports=True):
    ''' 输出.npy文件 '''
    return  np.save(filename, data, allow_pickle, fix_imports)


