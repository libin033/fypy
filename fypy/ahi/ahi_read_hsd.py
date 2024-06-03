# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : ahi_read_hsd.py

@Modify Time :  2022/11/10 14:45

@Author      : fypy Team

@Version     : 1.0

@Description :

'''
import os
import time
import numpy as np
from datetime import datetime, timedelta
import tempfile
import bz2
import dask.array as da
import xarray as xr
from contextlib import closing
import shutil
from subprocess import Popen, PIPE
from io import BytesIO
from .geometry import AreaDefinition

try:
    from shutil import which
except ImportError:
    # python 2 - won't be used, but needed for mocking in tests
    which = None

AHI_CHANNEL_NAMES = ("1", "2", "3", "4", "5",
                     "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16")

# logger = logging.getLogger('ahi_hsd')

# Basic information block:
_BASIC_INFO_TYPE = np.dtype([("hblock_number", "u1"),
                             ("blocklength", "<u2"),
                             ("total_number_of_hblocks", "<u2"),
                             ("byte_order", "u1"),
                             ("satellite", "S16"),
                             ("proc_center_name", "S16"),
                             ("observation_area", "S4"),
                             ("other_observation_info", "S2"),
                             ("observation_timeline", "<u2"),
                             ("observation_start_time", "f8"),
                             ("observation_end_time", "f8"),
                             ("file_creation_time", "f8"),
                             ("total_header_length", "<u4"),
                             ("total_data_length", "<u4"),
                             ("quality_flag1", "u1"),
                             ("quality_flag2", "u1"),
                             ("quality_flag3", "u1"),
                             ("quality_flag4", "u1"),
                             ("file_format_version", "S32"),
                             ("file_name", "S128"),
                             ("spare", "S40"),
                             ])

# Data information block
_DATA_INFO_TYPE = np.dtype([("hblock_number", "u1"),
                            ("blocklength", "<u2"),
                            ("number_of_bits_per_pixel", "<u2"),
                            ("number_of_columns", "<u2"),
                            ("number_of_lines", "<u2"),
                            ("compression_flag_for_data", "u1"),
                            ("spare", "S40"),
                            ])

# Projection information block
# See footnote 2; LRIT/HRIT Global Specification Section 4.4, CGMS, 1999)
_PROJ_INFO_TYPE = np.dtype([("hblock_number", "u1"),
                            ("blocklength", "<u2"),
                            ("sub_lon", "f8"),
                            ("CFAC", "<u4"),
                            ("LFAC", "<u4"),
                            ("COFF", "f4"),
                            ("LOFF", "f4"),
                            ("distance_from_earth_center", "f8"),
                            ("earth_equatorial_radius", "f8"),
                            ("earth_polar_radius", "f8"),
                            ("req2_rpol2_req2", "f8"),
                            ("rpol2_req2", "f8"),
                            ("req2_rpol2", "f8"),
                            ("coeff_for_sd", "f8"),
                            # Note: processing center use only:
                            ("resampling_types", "<i2"),
                            # Note: processing center use only:
                            ("resampling_size", "<i2"),
                            ("spare", "S40"),
                            ])

# Navigation information block
_NAV_INFO_TYPE = np.dtype([("hblock_number", "u1"),
                           ("blocklength", "<u2"),
                           ("navigation_info_time", "f8"),
                           ("SSP_longitude", "f8"),
                           ("SSP_latitude", "f8"),
                           ("distance_earth_center_to_satellite", "f8"),
                           ("nadir_longitude", "f8"),
                           ("nadir_latitude", "f8"),
                           ("sun_position", "f8", (3,)),
                           ("moon_position", "f8", (3,)),
                           ("spare", "S40"),
                           ])

# Calibration information block
_CAL_INFO_TYPE = np.dtype([("hblock_number", "u1"),
                           ("blocklength", "<u2"),
                           ("band_number", "<u2"),
                           ("central_wave_length", "f8"),
                           ("valid_number_of_bits_per_pixel", "<u2"),
                           ("count_value_error_pixels", "<u2"),
                           ("count_value_outside_scan_pixels", "<u2"),
                           ("gain_count2rad_conversion", "f8"),
                           ("offset_count2rad_conversion", "f8"),
                           ])

# Infrared band (Band No. 7 – 16)
# (Band No. 2 – 5: backup operation (See Table 4 bb))
_IRCAL_INFO_TYPE = np.dtype([("c0_rad2tb_conversion", "f8"),
                             ("c1_rad2tb_conversion", "f8"),
                             ("c2_rad2tb_conversion", "f8"),
                             ("c0_tb2rad_conversion", "f8"),
                             ("c1_tb2rad_conversion", "f8"),
                             ("c2_tb2rad_conversion", "f8"),
                             ("speed_of_light", "f8"),
                             ("planck_constant", "f8"),
                             ("boltzmann_constant", "f8"),
                             ("spare", "S40"),
                             ])

# Visible, near-infrared band (Band No. 1 – 6)
# (Band No. 1: backup operation (See Table 4 bb))
_VISCAL_INFO_TYPE = np.dtype([("coeff_rad2albedo_conversion", "f8"),
                              ("coeff_update_time", "f8"),
                              ("cali_gain_count2rad_conversion", "f8"),
                              ("cali_offset_count2rad_conversion", "f8"),
                              ("spare", "S80"),
                              ])

# 6 Inter-calibration information block
_INTER_CALIBRATION_INFO_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u2"),
    ("gsics_calibration_intercept", "f8"),
    ("gsics_calibration_slope", "f8"),
    ("gsics_calibration_coeff_quadratic_term", "f8"),
    ("gsics_std_scn_radiance_bias", "f8"),
    ("gsics_std_scn_radiance_bias_uncertainty", "f8"),
    ("gsics_std_scn_radiance", "f8"),
    ("gsics_correction_starttime", "f8"),
    ("gsics_correction_endtime", "f8"),
    ("gsics_radiance_validity_upper_lim", "f4"),
    ("gsics_radiance_validity_lower_lim", "f4"),
    ("gsics_filename", "S128"),
    ("spare", "S56"),
])

# 7 Segment information block
_SEGMENT_INFO_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u2"),
    ("total_number_of_segments", "u1"),
    ("segment_sequence_number", "u1"),
    ("first_line_number_of_image_segment", "u2"),
    ("spare", "S40"),
])

# 8 Navigation correction information block
_NAVIGATION_CORRECTION_INFO_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u2"),
    ("center_column_of_rotation", "f4"),
    ("center_line_of_rotation", "f4"),
    ("amount_of_rotational_correction", "f8"),
    ("numof_correction_info_data", "<u2"),
])

# 9 Observation time information block
_OBS_TIME_INFO_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u2"),
    ("number_of_observation_times", "<u2"),
])

# 10 Error information block
_ERROR_INFO_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u4"),
    ("number_of_error_info_data", "<u2"),
])

# 11 Spare block
_SPARE_TYPE = np.dtype([
    ("hblock_number", "u1"),
    ("blocklength", "<u2"),
    ("spare", "S256")
])



class ahi_read_hsd :

    def readhsd(self, hsdname, segnum):
        self.segment_number = segnum
        # self.segtotal = segtotal
        # self.data = None
        self.is_zipped = False

        self._unzipped = self.unzip_file(hsdname)
        if self._unzipped:
            # But if it is, set the filename to point to unzipped temp file
            self.is_zipped = True
            self.filename = self._unzipped
        else:
            self.filename = hsdname

        try:
            with open(self.filename) as fd:
                self.basic_info = np.fromfile(fd,
                                              dtype=_BASIC_INFO_TYPE,
                                              count=1)
                self.data_info = np.fromfile(fd,
                                             dtype=_DATA_INFO_TYPE,
                                             count=1)
                self.proj_info = np.fromfile(fd,
                                             dtype=_PROJ_INFO_TYPE,
                                             count=1)[0]
                self.nav_info = np.fromfile(fd,
                                            dtype=_NAV_INFO_TYPE,
                                            count=1)[0]
                fd.close()
        except BaseException as e :
            print(e)
            # os.remove(self.filename)
            return None

        self.platform_name = np2str(self.basic_info['satellite'])
        self.sensor = 'ahi'
        # self.segment_number = filename_info['segment']
        # self.total_segments = filename_info['total_segments']
        self.observation_area = np2str(self.basic_info['observation_area'])

        user_calibration=None
        calib_mode = 'NOMINAL'
        calib_mode_choices = ('NOMINAL', 'UPDATE')
        if calib_mode.upper() not in calib_mode_choices:
            raise ValueError('Invalid calibration mode: {}. Choose one of {}'.format(
                calib_mode, calib_mode_choices))

        self.calib_mode = calib_mode.upper()
        self.user_calibration = user_calibration

        with open(self.filename, "rb") as fp_:
            self._header = self._read_header(fp_)
            res = self._read_data(fp_, self._header)
            fp_.close()
        data = self.read_band(res)

        return data.values

    def read_band(self, res):

        res = self._mask_invalid(data=res, header=self._header)
        # Calibrate
        bandnum = self._header["block5"]['band_number'][0]
        if bandnum < 7:
            res = self.calibrate(res, 'reflectance', )
        else:
            res = self.calibrate(res, 'brightness_temperature')

        # Update metadata
        new_info = dict(
            # units=info['units'],
            # standard_name=info['standard_name'],
            # wavelength=info['wavelength'],
            resolution='resolution',
            # id=key,
            # name=key['name'],
            # scheduled_time=self.scheduled_time,
            platform_name=self.platform_name,
            sensor=self.sensor,
            satellite_longitude=float(self.nav_info['SSP_longitude']),
            satellite_latitude=float(self.nav_info['SSP_latitude']),
            satellite_altitude=float(self.nav_info['distance_earth_center_to_satellite'] -
                                     self.proj_info['earth_equatorial_radius']) * 1000,
            orbital_parameters={
                'projection_longitude': float(self.proj_info['sub_lon']),
                'projection_latitude': 0.,
                'projection_altitude': float(self.proj_info['distance_from_earth_center'] -
                                             self.proj_info['earth_equatorial_radius']) * 1000,
                # 'satellite_actual_longitude': actual_lon,
                # 'satellite_actual_latitude': actual_lat,
                # 'satellite_actual_altitude': actual_alt,
                'nadir_longitude': float(self.nav_info['nadir_longitude']),
                'nadir_latitude': float(self.nav_info['nadir_latitude'])}
        )
        res = xr.DataArray(res, attrs=new_info, dims=['y', 'x'])

        self.get_area_def()

        self.mask_space = True
        # Mask space pixels
        if self.mask_space:
            res = self._mask_space(res)

        return res

    def calibrate(self, data, calibration):
        """Calibrate the data."""
        tic = datetime.now()

        if calibration == 'counts':
            return data

        if calibration in ['radiance', 'reflectance', 'brightness_temperature']:
            data = self.convert_to_radiance(data)

        if calibration == 'reflectance':
            data = self._vis_calibrate(data)
        elif calibration == 'brightness_temperature':
            data = self._ir_calibrate(data)

        # print("Calibration time " + str(datetime.now() - tic))
        return data

    def convert_to_radiance(self, data):
        """Calibrate to radiance."""
        bnum = self._header["block5"]['band_number'][0]
        # Check calibration mode and select corresponding coefficients
        if self.calib_mode == "UPDATE" and bnum < 7:
            dn_gain = self._header['calibration']["cali_gain_count2rad_conversion"][0]
            dn_offset = self._header['calibration']["cali_offset_count2rad_conversion"][0]
            if dn_gain == 0 and dn_offset == 0:
                print("No valid updated coefficients, fall back to default values.")
                dn_gain = self._header["block5"]["gain_count2rad_conversion"][0]
                dn_offset = self._header["block5"]["offset_count2rad_conversion"][0]
        else:
            dn_gain = self._header["block5"]["gain_count2rad_conversion"][0]
            dn_offset = self._header["block5"]["offset_count2rad_conversion"][0]

        # Assume no user correction
        correction_type = None
        if isinstance(self.user_calibration, dict):
            # Check if we have DN correction coeffs
            if 'type' in self.user_calibration:
                correction_type = self.user_calibration['type']
            else:
                # If not, assume radiance correction
                correction_type = 'RAD'
            if correction_type == 'DN':
                # Replace file calibration with user calibration
                dn_gain, dn_offset = get_user_calibration_factors(self.band_name,
                                                                  self.user_calibration)
            elif correction_type == 'RAD':
                user_slope, user_offset = get_user_calibration_factors(self.band_name,
                                                                       self.user_calibration)

        data = (data * dn_gain + dn_offset).clip(0)
        # If using radiance correction factors from GSICS or similar, apply here
        if correction_type == 'RAD':
            data = apply_rad_correction(data, user_slope, user_offset)
        return data

    def _vis_calibrate(self, data):
        """Visible channel calibration only."""
        coeff = self._header["calibration"]["coeff_rad2albedo_conversion"]
        return (data * coeff * 100).clip(0)

    def _ir_calibrate(self, data):
        """IR calibration."""
        # No radiance -> no temperature
        data = da.where(data == 0, np.float32(np.nan), data)

        cwl = self._header['block5']["central_wave_length"][0] * 1e-6
        c__ = self._header['calibration']["speed_of_light"][0]
        h__ = self._header['calibration']["planck_constant"][0]
        k__ = self._header['calibration']["boltzmann_constant"][0]
        a__ = (h__ * c__) / (k__ * cwl)

        b__ = ((2 * h__ * c__ ** 2) / (data * 1.0e6 * cwl ** 5)) + 1

        Te_ = a__ / da.log(b__)

        c0_ = self._header['calibration']["c0_rad2tb_conversion"][0]
        c1_ = self._header['calibration']["c1_rad2tb_conversion"][0]
        c2_ = self._header['calibration']["c2_rad2tb_conversion"][0]

        return (c0_ + c1_ * Te_ + c2_ * Te_ ** 2).clip(0)

    def _read_header(self, fp_):
        """Read header."""
        header = {}

        fpos = 0
        header['block1'] = np.fromfile(
            fp_, dtype=_BASIC_INFO_TYPE, count=1)
        fpos = fpos + int(header['block1']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block1')
        fp_.seek(fpos, 0)
        header["block2"] = np.fromfile(fp_, dtype=_DATA_INFO_TYPE, count=1)
        fpos = fpos + int(header['block2']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block2')
        fp_.seek(fpos, 0)
        header["block3"] = np.fromfile(fp_, dtype=_PROJ_INFO_TYPE, count=1)
        fpos = fpos + int(header['block3']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block3')
        fp_.seek(fpos, 0)
        header["block4"] = np.fromfile(fp_, dtype=_NAV_INFO_TYPE, count=1)
        fpos = fpos + int(header['block4']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block4')
        fp_.seek(fpos, 0)
        header["block5"] = np.fromfile(fp_, dtype=_CAL_INFO_TYPE, count=1)
        # print("Band number = " + str(header["block5"]['band_number'][0]))
        # print('Time_interval: %s - %s', str(self.start_time), str(self.end_time))
        band_number = header["block5"]['band_number'][0]
        if band_number < 7:
            cal = np.fromfile(fp_, dtype=_VISCAL_INFO_TYPE, count=1)
        else:
            cal = np.fromfile(fp_, dtype=_IRCAL_INFO_TYPE, count=1)
        fpos = fpos + int(header['block5']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block5')
        fp_.seek(fpos, 0)

        header['calibration'] = cal

        header["block6"] = np.fromfile(
            fp_, dtype=_INTER_CALIBRATION_INFO_TYPE, count=1)
        fpos = fpos + int(header['block6']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block6')
        fp_.seek(fpos, 0)
        header["block7"] = np.fromfile(
            fp_, dtype=_SEGMENT_INFO_TYPE, count=1)
        fpos = fpos + int(header['block7']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block7')
        fp_.seek(fpos, 0)
        header["block8"] = np.fromfile(
            fp_, dtype=_NAVIGATION_CORRECTION_INFO_TYPE, count=1)
        # 8 The navigation corrections:
        ncorrs = header["block8"]['numof_correction_info_data'][0]
        dtype = np.dtype([
            ("line_number_after_rotation", "<u2"),
            ("shift_amount_for_column_direction", "f4"),
            ("shift_amount_for_line_direction", "f4"),
        ])
        corrections = []
        for _i in range(ncorrs):
            corrections.append(np.fromfile(fp_, dtype=dtype, count=1))
        fpos = fpos + int(header['block8']['blocklength'])
        self._check_fpos(fp_, fpos, 40, 'block8')
        fp_.seek(fpos, 0)
        header['navigation_corrections'] = corrections
        header["block9"] = np.fromfile(fp_,
                                       dtype=_OBS_TIME_INFO_TYPE,
                                       count=1)
        numobstimes = header["block9"]['number_of_observation_times'][0]

        dtype = np.dtype([
            ("line_number", "<u2"),
            ("observation_time", "f8"),
        ])
        lines_and_times = []
        for _i in range(numobstimes):
            lines_and_times.append(np.fromfile(fp_,
                                               dtype=dtype,
                                               count=1))
        header['observation_time_information'] = lines_and_times
        fpos = fpos + int(header['block9']['blocklength'])
        self._check_fpos(fp_, fpos, 40, 'block9')
        fp_.seek(fpos, 0)

        header["block10"] = np.fromfile(fp_,
                                        dtype=_ERROR_INFO_TYPE,
                                        count=1)
        dtype = np.dtype([
            ("line_number", "<u2"),
            ("numof_error_pixels_per_line", "<u2"),
        ])
        num_err_info_data = header["block10"][
            'number_of_error_info_data'][0]
        err_info_data = []
        for _i in range(num_err_info_data):
            err_info_data.append(np.fromfile(fp_, dtype=dtype, count=1))
        header['error_information_data'] = err_info_data
        fpos = fpos + int(header['block10']['blocklength'])
        self._check_fpos(fp_, fpos, 40, 'block10')
        fp_.seek(fpos, 0)

        header["block11"] = np.fromfile(fp_, dtype=_SPARE_TYPE, count=1)
        fpos = fpos + int(header['block11']['blocklength'])
        self._check_fpos(fp_, fpos, 0, 'block11')
        fp_.seek(fpos, 0)

        return header

    def _read_data(self, fp_, header):
        """Read data block."""

        nlines = int(header["block2"]['number_of_lines'][0])
        ncols = int(header["block2"]['number_of_columns'][0])
        return da.from_array(np.memmap(self.filename, offset=fp_.tell(),
                                       dtype='<u2', shape=(nlines, ncols), mode='r'),
                             chunks=4096)

    def unzip_file(self, filename, tmppath=None):
        from subprocess import Popen, PIPE
        from io import BytesIO
        from contextlib import closing
        import shutil
        import bz2

        try:
            from shutil import which
        except ImportError:
            # python 2 - won't be used, but needed for mocking in tests
            which = None

        """Unzip the file if file is bzipped = ending with 'bz2'."""
        try:
            if filename.endswith('bz2'):
                # fdn, tmpfilepath = tempfile.mkstemp()
                if tmppath is None :
                    tmpfilepath = filename.replace('.bz2','')
                else:
                    tmpfilepath = os.path.join(tmppath, os.path.basename(filename).replace('.bz2',''))

                if os.path.isfile(tmpfilepath) :
                    return tmpfilepath
                # print("解压bz2文件【%s】" %(tmpfilepath))
                # try pbzip2
                pbzip = which('pbzip2')
                # Run external pbzip2
                if pbzip is not None:
                    n_thr = os.environ.get('OMP_NUM_THREADS')
                    if n_thr:
                        runner = [pbzip,
                                  '-dc',
                                  '-p'+str(n_thr),
                                  filename]
                    else:
                        runner = [pbzip,
                                  '-dc',
                                  filename]
                    p = Popen(runner, stdout=PIPE, stderr=PIPE)
                    stdout = BytesIO(p.communicate()[0])
                    status = p.returncode
                    if status != 0:
                        raise IOError("pbzip2 error '%s', failed, status=%d"
                                      % (filename, status))
                    with closing(open(tmpfilepath, 'wb')) as ofpt:
                        try:
                            stdout.seek(0)
                            shutil.copyfileobj(stdout, ofpt)
                        except IOError:
                            import traceback
                            traceback.print_exc()
                            print("Failed to read bzipped file %s",
                                  str(filename))
                            os.remove(tmpfilepath)

                    return tmpfilepath

                # Otherwise, fall back to the original method
                bz2file = bz2.BZ2File(filename)
                with closing(open(tmpfilepath, 'wb')) as ofpt:
                    try:
                        ofpt.write(bz2file.read())
                    except IOError:
                        import traceback
                        traceback.print_exc()
                        print("Failed to read bzipped file %s", str(filename))
                        # os.remove(tmpfilepath)
                        return None
                return tmpfilepath
            else:
                return None
        except BaseException :
            return None

    def _check_fpos(self, fp_, fpos, offset, block):
        """Check file position matches blocksize."""
        if (fp_.tell() + offset != fpos):
            print("Actual "+block+" header size does not match expected")
        return

    def _mask_invalid(self, data, header):
        """Mask invalid data."""

        invalid = da.logical_or(data == header['block5']["count_value_outside_scan_pixels"][0],
                                data == header['block5']["count_value_error_pixels"][0])
        return da.where(invalid, np.float32(np.nan), data)

    def _mask_space(self, data):
        """Mask space pixels."""
        return data.where(get_geostationary_mask(self.area))
        # pass

    @property
    def start_time(self):
        """Get the start time."""
        return datetime(1858, 11, 17) + timedelta(days=float(self.basic_info['observation_start_time']))

    @property
    def end_time(self):
        """Get the end time."""
        return datetime(1858, 11, 17) + timedelta(days=float(self.basic_info['observation_end_time']))

    def get_area_def(self):
        """Get the area definition."""

        pdict = {}
        pdict['cfac'] = np.uint32(self.proj_info['CFAC'])
        pdict['lfac'] = np.uint32(self.proj_info['LFAC'])
        pdict['coff'] = np.float32(self.proj_info['COFF'])
        pdict['loff'] = -np.float32(self.proj_info['LOFF']) + 1
        pdict['a'] = float(self.proj_info['earth_equatorial_radius'] * 1000)
        pdict['h'] = float(self.proj_info['distance_from_earth_center'] * 1000 - pdict['a'])
        pdict['b'] = float(self.proj_info['earth_polar_radius'] * 1000)
        pdict['ssp_lon'] = float(self.proj_info['sub_lon'])
        pdict['nlines'] = int(self.data_info['number_of_lines'])
        pdict['ncols'] = int(self.data_info['number_of_columns'])
        pdict['scandir'] = 'N2S'

        pdict['loff'] = pdict['loff'] + (self.segment_number * pdict['nlines'])

        aex = get_area_extent(pdict)

        pdict['a_name'] = self.observation_area
        pdict['a_desc'] = "AHI {} area".format(self.observation_area)
        pdict['p_id'] = 'geosh8'

        area = get_area_definition(pdict, aex)

        self.area = area
        return area

    def __del__(self):

        # time.sleep(2)
        # if (self.is_zipped and os.path.exists(self.filename)):
        #     os.remove(self.filename)
        pass




def np2str(value):
    """Convert an `numpy.string_` to str.

    Args:
        value (ndarray): scalar or 1-element numpy array to convert

    Raises:
        ValueError: if value is array larger than 1-element or it is not of
                    type `numpy.string_` or it is not a numpy array

    """
    if hasattr(value, 'dtype') and \
            issubclass(value.dtype.type, (np.str_, np.string_, np.object_)) \
            and value.size == 1:
        value = value.item()
        if not isinstance(value, str):
            # python 3 - was scalar numpy array of bytes
            # otherwise python 2 - scalar numpy array of 'str'
            value = value.decode()
        return value
    else:
        raise ValueError("Array is not a string type or is larger than 1")

def get_user_calibration_factors(band_name, correction_dict):
    """Retrieve radiance correction factors from user-supplied dict."""
    if band_name in correction_dict:
        try:
            slope = correction_dict[band_name]['slope']
            offset = correction_dict[band_name]['offset']
        except KeyError:
            raise KeyError("Incorrect correction factor dictionary. You must "
                           "supply 'slope' and 'offset' keys.")
    else:
        # If coefficients not present, warn user and use slope=1, offset=0
        print("WARNING: You have selected radiance correction but "
                      " have not supplied coefficients for channel " +
                      band_name)
        return 1., 0.

    return slope, offset

def apply_rad_correction(data, slope, offset):
    """Apply GSICS-like correction factors to radiance data."""
    data = (data - offset) / slope
    return data

def get_area_extent(pdict):
    """Get the area extent seen by a geostationary satellite.

    Args:
        pdict: A dictionary containing common parameters:
            nlines: Number of lines in image
            ncols: Number of columns in image
            cfac: Column scaling factor
            lfac: Line scaling factor
            coff: Column offset factor
            loff: Line offset factor
            scandir: 'N2S' for standard (N->S), 'S2N' for inverse (S->N)
    Returns:
        aex: An area extent for the scene

    """
    # count starts at 1
    cols = 1 - 0.5

    if pdict['scandir'] == 'S2N':
        lines = 0.5 - 1
        scanmult = -1
    else:
        lines = 1 - 0.5
        scanmult = 1
    # Lower left x, y scanning angles in degrees
    ll_x, ll_y = get_xy_from_linecol(lines * scanmult,
                                     cols,
                                     (pdict['loff'], pdict['coff']),
                                     (pdict['lfac'], pdict['cfac']))

    cols += pdict['ncols']
    lines += pdict['nlines']
    # Upper right x, y scanning angles in degrees
    ur_x, ur_y = get_xy_from_linecol(lines * scanmult,
                                     cols,
                                     (pdict['loff'], pdict['coff']),
                                     (pdict['lfac'], pdict['cfac']))
    if pdict['scandir'] == 'S2N':
        ll_y *= -1
        ur_y *= -1

    # Convert degrees to radians and create area extent
    aex = make_ext(ll_x=ll_x, ur_x=ur_x, ll_y=ll_y, ur_y=ur_y, h=pdict['h'])

    return aex


def get_area_definition(pdict, a_ext):
    """Get the area definition for a geo-sat.

    Args:
        pdict: A dictionary containing common parameters:
            nlines: Number of lines in image
            ncols: Number of columns in image
            ssp_lon: Subsatellite point longitude (deg)
            a: Earth equatorial radius (m)
            b: Earth polar radius (m)
            h: Platform height (m)
            a_name: Area name
            a_desc: Area description
            p_id: Projection id
        a_ext: A four element tuple containing the area extent (scan angle)
               for the scene in radians
    Returns:
        a_def: An area definition for the scene

    .. note::

        The AreaDefinition `proj_id` attribute is being deprecated.

    """
    proj_dict = {'a': float(pdict['a']),
                 'b': float(pdict['b']),
                 'lon_0': float(pdict['ssp_lon']),
                 'h': float(pdict['h']),
                 'proj': 'geos',
                 'units': 'm'}

    a_def = AreaDefinition(
        pdict['a_name'],
        pdict['a_desc'],
        pdict['p_id'],
        proj_dict,
        int(pdict['ncols']),
        int(pdict['nlines']),
        a_ext)

    return a_def

def get_xy_from_linecol(line, col, offsets, factors):
    """Get the intermediate coordinates from line & col.

    Intermediate coordinates are actually the instruments scanning angles.
    """
    loff, coff = offsets
    lfac, cfac = factors
    x__ = float(col - coff) / (float(cfac) / 2 ** 16)
    y__ = float(line - loff) / (float(lfac) / 2 ** 16)

    return x__, y__


def make_ext(ll_x, ur_x, ll_y, ur_y, h):
    """Create the area extent from computed ll and ur.

    Args:
        ll_x: The lower left x coordinate (m)
        ur_x: The upper right x coordinate (m)
        ll_y: The lower left y coordinate (m)
        ur_y: The upper right y coordinate (m)
        h: The satellite altitude above the Earth's surface
    Returns:
        aex: An area extent for the scene

    """
    aex = (np.deg2rad(ll_x) * h, np.deg2rad(ll_y) * h,
           np.deg2rad(ur_x) * h, np.deg2rad(ur_y) * h)

    return aex


def get_geostationary_mask(area):
    """Compute a mask of the earth's shape as seen by a geostationary satellite.

    Args:
        area (pyresample.geometry.AreaDefinition) : Corresponding area
                                                    definition

    Returns:
        Boolean mask, True inside the earth's shape, False outside.

    """
    # Compute projection coordinates at the earth's limb
    h = area.proj_dict['h']
    xmax, ymax = get_geostationary_angle_extent(area)
    xmax *= h
    ymax *= h

    # Compute projection coordinates at the centre of each pixel
    x, y = area.get_proj_coords(chunks=4096)

    # Compute mask of the earth's elliptical shape
    return ((x / xmax) ** 2 + (y / ymax) ** 2) <= 1

def get_geostationary_angle_extent(geos_area):
    """Get the max earth (vs space) viewing angles in x and y."""
    # TODO: take into account sweep_axis_angle parameter

    # get some projection parameters
    try:
        crs = geos_area.crs
        a = crs.ellipsoid.semi_major_metre
        b = crs.ellipsoid.semi_minor_metre
        if np.isnan(b):
            # see https://github.com/pyproj4/pyproj/issues/457
            raise AttributeError("'semi_minor_metre' attribute is not valid "
                                 "in older versions of pyproj.")
    except AttributeError:
        # older versions of pyproj don't have CRS objects
        from pyresample.utils import proj4_radius_parameters
        a, b = proj4_radius_parameters(geos_area.proj_dict)

    req = float(a) / 1000
    rp = float(b) / 1000
    h = float(geos_area.proj_dict['h']) / 1000 + req

    # compute some constants
    aeq = 1 - req**2 / (h ** 2)
    ap_ = 1 - rp**2 / (h ** 2)

    # generate points around the north hemisphere in satellite projection
    # make it a bit smaller so that we stay inside the valid area
    xmax = np.arccos(np.sqrt(aeq))
    ymax = np.arccos(np.sqrt(ap_))
    return xmax, ymax
