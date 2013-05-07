"""
Routines for reading IQ data from IRIS TS files (RAW IRIS files).

"""

import os

import numpy as np

import decode_highsnr


class IrisTSFile:
    """
    A class for working with IRIS TS files which hold IQ data.

    This class reads the global and pulse specific metadata from a file
    formatted in the IRIS TS record data format.  Pulse data can be retrieved
    using the `read_pulse_data` and `read_all_pulse_data` methods.

    Parameters
    ----------
    filename : str
        Filename of IRIS TS file to read data from.

    Attributes
    ----------
    filename : str
        Filename of IRIS TS file being read.
    npulses : int
        Number of pulses in the file.
    pulse_info : dic
        Dictionary of metadata contained in the TS file pulseinfo section.
    pulse_headers : list of dics
        Lists of length npulses which contain metadata dictionaries for
        a given pulse.
    npoints : array
        Array of length npulses, which specify the number of points in a
        given pulse.
    max_points : int
        Maximum number of points in any pulse.

    """
    # private attributes
    # _f : file
    #   Open file from which data is read from.
    # _data_locations: array
    #   Seek locations in _f from which pulse data can be read.

    def __init__(self, filename):
        """ Initalize. """

        # open the file for reading
        self.filename = filename
        f = open(filename, 'r')
        self._f = open(filename, 'r')

        # read in the pulse info, store as pulse_info attribute
        self.pulse_info = self._get_rvptspulseinfo_dic(f)

        # read through the file, for each pulse record the header, the number
        # of points in that pulse and the location of the pulse data.
        self.pulse_headers = []
        data_locations = []
        npoints = []

        while self._f.tell() != os.stat(filename).st_size:

            # read and store the pulse header
            current_pulsehdr = self._get_rvptspulsehdr_dic(self._f)
            self.pulse_headers.append(current_pulsehdr)

            # record the location of the pulse data
            data_locations.append(self._f.tell())

            # record the number of points in the pulse data
            numvecs = int(current_pulsehdr['iNumVecs'])
            viqperbin = int(current_pulsehdr['iVIQPerBin'])
            npoints.append(numvecs * viqperbin * 2)

            # skip over the pulse data
            bytes_in_current_pulse = 2 * 2 * numvecs * viqperbin
            self._f.seek(bytes_in_current_pulse, 1)

        # set pulse attributes
        self._data_locations = np.array(data_locations, dtype='int64')
        self.npoints = np.array(npoints, dtype='int32')
        self.npulses = len(self.npoints)
        self.max_points = np.max(self.npoints)
        return

    def _get_rvptspulseinfo_dic(self, f):
        """ Read an ASCII rvptsPulseInfo structure, returns a dict. """
        dic = {}
        while True:
            line = f.readline()
            if '=' in line:
                key, value = line.split('=')
                dic[key] = value.strip('\n')
            if line.startswith('rvptsPulseInfo end'):
                break
        return dic

    def _get_rvptspulsehdr_dic(self, f):
        """ Read an ASCII rvptsPulseHdr structure, returns a dict. """
        dic = {}
        while True:
            line = f.readline()
            if '=' in line:
                key, value = line.split('=')
                dic[key] = value.strip('\n')
            if line.startswith('rvptsPulseHdr end'):
                break
        return dic

    def read_pulse_data(self, pulse_id):
        """
        Read and decode pulse data.

        Parameters
        ----------
        pulse_id : int
            Pulse number (0 based) to retrieve data from.

        Returns
        -------
        data : array
            Array of raw pulse data decoded to float32 format.  This is the
            raw data, IQ is still interleaved and if there are multiple
            channels they are still appended together.

        """
        if pulse_id < 0 or pulse_id >= self.npulses:
            raise ValueError('invalid pulse_id')

        # determine number of points in the pulse
        n_points = self.npoints[pulse_id]

        # read pulse data from file
        data_location = self._data_locations[pulse_id]
        self._f.seek(data_location)
        buf = self._f.read(n_points * 2)
        highsnr_data = np.fromstring(buf, np.uint16)

        # return decoded data
        return decode_highsnr.decode_highsnr(highsnr_data, n_points)
        # if the decode_highsnr module is not available, the python
        # implementation can be used, but with significantly slower
        # performance
        # return self._decode_highsnr_python(highsnr_data, n_points)

    def read_all_pulse_data(self):
        """
        Read and decode all pulse data.

        Returns
        -------
        data : array
            2D array of raw pulse data decoded to float32 format.  This is the
            raw data, IQ is still interleaved and if there are multiple
            channels they are still appended together.

        """
        data = np.empty((self.npulses, self.max_points), dtype=np.float32)
        for i, npoints in enumerate(self.npoints):
            data[i, :npoints] = self.read_pulse_data(i)
        return data

    def _decode_highsnr_python(self, highsnr_data, npoints):
        """ Decode a vector encoded in High SNR format. """
        # a much faster implementation of this function is available in the
        # Cython module decode_highsnr.  This function is kept for reference.
        decoded_data = np.empty((npoints), dtype=np.float32)
        for i, p in enumerate(highsnr_data):
            decoded_data[i] = self._decode_highsnr_point_python(p)
        return decoded_data

    def _decode_highsnr_point_python(self, p, debug=False):
        """ Decode a point encoded in High SNR format. """
        # A much faster implementation of this function is available in the
        # Cython module decode_highsnr. This function is kept for reference.

        # find the sign_bit
        #sign_bit = (p & (2 ** 11)) / (2 ** 11)
        #sign_bit = (p & (2048)) / (2048)
        sign_bit = p >> 11 & 1

        # find the exponent
        #exponent = (p & (2 ** 12 + 2 ** 13 + 2 ** 14 + 2 ** 15)) / (2 ** 12)
        #exponent = (p & (61440)) / (4096)
        exponent = p >> 12

        # find the mantissa
        #mantissa = p & (2*11-1) # and with 0000 0 11111111111
        mantissa = p & 2047

        if debug:
            print "exponent:", exponent
            print "sign_bit:", sign_bit
            print "mantissa:", mantissa

        if exponent != 0:
            if sign_bit:  # sign == 1
                # 13-bit signed int with leading 10
                # calculate two's compliment
                #val = -1 * (2 ** 12 - mantissa)
                val = -4096 + mantissa
            else:  # sign == 0
                # 13-bit signed int with leading 01
                #val = 2 ** 11 + mantissa
                val = 2048 + mantissa
            return val * 2. ** (exponent - 25)
        else:  # exponent == 0
            if sign_bit:
                # two compliment
                val = -2048 + mantissa
            else:
                val = mantissa
            return val * 2. ** -24


def write_iristsfile_to_netcdf(filename, iristsfile, format='NETCDF4'):
    """
    Write the IQ data from a IrisTSFile object to a netCDF file.

    Parameters
    ----------
    filename : str
        Name of file to write
    iristsfile : IrisTSFile
        IrisTSFile object to write to netCDF file.
    format : str
        Format of the netCDF file to create.

    """
    import netCDF4

    # create the dataset
    dataset = netCDF4.Dataset(filename, 'w', format=format)

    npulses = iristsfile.npulses
    nsamples = iristsfile.max_points / 2
    if iristsfile.pulse_info['iPolarization'] == '3':
        dual_polarized = True
        nsamples /= 2
    else:
        dual_polarized = False

    # create the pulse and sample dimensions
    dataset.createDimension('pulse', npulses)
    dataset.createDimension('sample', nsamples)

    # create the I and Q variables
    comment = 'First data point in each pulse is from the burst pulse'
    i0 = dataset.createVariable('I_0', 'f4', ('pulse', 'sample'))
    i0.setncattr('long_name', 'I data for first channel')
    i0.setncattr('units', 'voltage')
    i0.setncattr('comment', comment)
    q0 = dataset.createVariable('Q_0', 'f4', ('pulse', 'sample'))
    q0.setncattr('long_name', 'Q data for first channel')
    q0.setncattr('units', 'voltage')
    q0.setncattr('comment', comment)
    if dual_polarized:
        i1 = dataset.createVariable('I_1', 'f4', ('pulse', 'sample'))
        i1.setncattr('long_name', 'I data for second channel')
        i1.setncattr('units', 'voltage')
        i1.setncattr('comment', comment)
        q1 = dataset.createVariable('Q_1', 'f4', ('pulse', 'sample'))
        q1.setncattr('long_name', 'Q data for second channel')
        q1.setncattr('units', 'voltage')
        q1.setncattr('comment', comment)

    # populate I and Q variables
    # TODO differing number of samples per pulse
    for pulse_id in xrange(npulses):
        pulse_data = iristsfile.read_pulse_data(pulse_id)
        i_data = pulse_data[::2]    # even samples
        q_data = pulse_data[1::2]   # odd samples
        if dual_polarized:
            half = len(i_data) / 2
            i0[pulse_id] = i_data[:half]    # first half
            q0[pulse_id] = q_data[:half]    # second half
            i1[pulse_id] = i_data[half:]    # first half
            q1[pulse_id] = q_data[half:]    # second half
        else:
            i0[pulse_id] = i_data[:]
            q0[pulse_id] = q_data[:]

    metadata = _create_metadata(iristsfile.pulse_info)
    for k, v in metadata.iteritems():
        if type(v) == dict:
            _create_ncvar(v, dataset, k, ())
        else:
            dataset.setncattr(k, v)

    pulse_metadata = _create_pulse_metadata(iristsfile.pulse_headers)
    for k, v in pulse_metadata.iteritems():
        _create_ncvar(v, dataset, k, ('pulse', ))

    dataset.close()


def _create_metadata(pulse_info):
    """ Create a dictionary of metadata from a pulse_info dictionary. """

    def float_array(key):
        """ Create array of floats. """
        return np.array(float(pulse_info[key]), dtype='float32')

    d = {}
    info_keys = pulse_info.keys()
    if 'iPolarization' in info_keys:
        modes = {'0': 'Horizontal',
                 '1': 'Vertical',
                 '2': 'Alternating_horizontal_and_vertical',
                 '3': 'Horizontal_and_vertical'}
        d['polarization'] = modes[pulse_info['iPolarization']]

    if 'sSiteName' in info_keys:
        d['site_name'] = pulse_info['sSiteName']

    if 'fPWidthUSec' in info_keys:
        d['pulse_width'] = {
            'data': float_array('fPWidthUSec'),
            'long_name': 'Pulse width',
            'units': 'microseconds',
        }

    if 'fDBzCalib' in info_keys:
        d['calibration_reflectivity'] = {
            'data': float_array('fDBzCalib'),
            'long_name': 'Calibration reflectivity',
            'units': 'dBz',
            'comment': 'calibration at 1 km'
        }

    if 'iSampleSize' in info_keys:
        d['pulses_per_ray'] = int(pulse_info['iSampleSize'])

    if 'fGdrOffset' in info_keys:
        d['total_gain_ratio'] = {
            'data': float_array('fGdrOffset'),
            'long_name': 'Total gain ratio of the reciever channels',
            'units': 'dB',
            'comment': 'Valid in co-reciever systems',
        }

    if 'fXdrOffset' in info_keys:
        d['reciever_gain_ratio'] = {
            'data': float_array('fXdrOffset'),
            'long_name': 'Ratio of reciever gain',
            'units': 'dB',
            'comment': 'Valid in dual reciever systems'
        }

    if 'fSyClkMHz' in info_keys:
        d['system_clock_frequency_mhz'] = float(pulse_info['fSyClkMHz'])

    if 'fWavelengthCM' in info_keys:
        d['wavelength_cm'] = float(pulse_info['fWavelengthCM'])

    if 'fSaturationDBM' in info_keys:
        d['saturation_power'] = {
            'data': float_array('fSaturationDBM'),
            'long_name': 'Saturation power',
            'units': 'dBm',
            'comment': ('Power in dBm corresponding to a power '
                        'magnitude of 1.0'),
        }

    # TODO Ranges
    if 'iRangeMask' in info_keys:
        r_str = r_str = pulse_info['iRangeMask'].split()
        u16r = np.array([int(i) for i in r_str], dtype='uint16')
        u8_0, u8_1 = divmod(u16r, 256)  # 2**8
        u8 = np.empty((len(r_str) * 2,),  dtype='uint8')
        u8[::2] = u8_0      # XXX these might be switched
        u8[1::2] = u8_1
        # need to determine array and uint8 bit order.

    if 'fNoiseDBm' in info_keys:
        noise = pulse_info['fNoiseDBm'].split()
        d['noise_level_0'] = {
            'data': np.array([float(noise[0])], dtype='float32'),
            'long_name': 'Noise level of first channel',
            'units': 'dBm'
        }

        if len(noise) > 1:
            d['noise_level_1'] = {
                'data': np.array([float(noise[1])], dtype='float32'),
                'long_name': 'Noise level for second channel',
                'units': 'dBm'
            }

    if 'fNoiseStdvDB' in info_keys:
        noise = pulse_info['fNoiseDBm'].split()
        d['noise_standard_deviation_0'] = {
            'data': np.array([float(noise[0])], dtype='float32'),
            'long_name': ('Standard deviation of noise measuremnt for '
                          'first channel'),
            'units': 'dB',
        }
        if len(noise) > 1:
            d['noise_standard_deviation_1'] = {
                'data': np.array([float(noise[1])], dtype='float32'),
                'long_name': ('Standard deviation of noise measurement '
                              'for second channel'),
                'units': 'dB',
            }

    if 'fNoiseRangeKM' in info_keys:
        d['noise_range'] = {
            'data': float_array('fNoiseRangeKM'),
            'long_name': 'Range where noise sample was taken',
            'units': 'km',
        }

    if 'fNoisePRFHz' in info_keys:
        d['noise_prf'] = {
            'data': float_array('fNoisePRFHz'),
            'long_name': 'PRF rate when noise sample taken',
            'units': 'hz',
        }

    if 'sVersionString' in info_keys:
        d['iris_ts_file_version'] = pulse_info['sVersionString']

    return d


def _create_pulse_metadata(pulse_headers):
    """ Create a dictionary of pulse metadata from pulse headers. """

    def int_array(key):
        """ Create array of ints. """
        return np.array([int(d[key]) for d in pulse_headers], dtype='int32')

    def degree_array(key):
        """ Create array of degrees measurements. """
        d = np.array([float(d[key]) for d in pulse_headers], dtype='float32')
        d *= 180. / 32768.
        return d

    def volt_array(key):
        """ Create array of volt measurements. """
        d = np.array([float(d[key]) for d in pulse_headers], dtype='float32')
        d *= 0.446  # 0.446 volts = 1
        return d

    hdr_keys = pulse_headers[0].keys()
    d = {}
    if 'iSeqNum' in hdr_keys:
        d['sequence_number'] = {
            'data': int_array('iSeqNum'),
            'long_name': 'Sequence Number',
            'comment': ('Pulse sequence number, non-sequential values '
                        'indicate pulses are missing'),
        }

    if 'iFlags' in hdr_keys:
        d['flags'] = {
            'data': int_array('iFlags'),
            'long_name': 'Pulse status flags',
            'comment': ('bit 1 set if data is valid\n'
                        'bit 2 set if prior pulses are missing\n'
                        'bit 3 set if first pulse in a trigger bank\n'
                        'bit 4 set if last pulse in a trigger bank\n'
                        'bit 5 set if trigger bank is beginning\n'
                        'bit 6 set if triggers were blanked'),
        }

    if 'iTimeUTC' in hdr_keys:
        d['time'] = {
            'data': int_array('iTimeUTC'),
            'long_name': 'Time in seconds since the Epoch',
            'units': 'seconds since 1970-01-01T00:00:00Z',
            'comment': ('Seconds since the UNIX epoch.\n'
                        'Must be combined with time_nanoseconds to '
                        'determine the time of the pulse'),
        }

    if 'iNanoUTC' in hdr_keys:
        d['time_nanoseconds'] = {
            'data': int_array('iNanoUTC'),
            'long_name': 'Nanoseconds since time',
            'units': 'nanoseconds since time',
            'ancillary_variables': 'time',
        }

    if 'iTxPhase' in hdr_keys:
        d['transmit_phase'] = {
            'data': degree_array('iTxPhase'),
            'long_name': 'Transmit phase',
            'units': 'degrees',
        }

    if 'iPedAz' in hdr_keys:
        d['pedestal_azimuth'] = {
            'data': degree_array('iPedAz'),
            'long_name': 'Pedestal azimuth angle',
            'units': 'degrees',
        }

    if 'iPedEl' in hdr_keys:
        d['pedestal_elevation'] = {
            'data': degree_array('iPedEl'),
            'long_name': 'Pedestal elevation angle',
            'units': 'degrees',
        }

    if 'iAzV' in hdr_keys:
        d['azimuth_velocity'] = {
            'data': degree_array('iAzV'),
            'long_name': 'Azimuth velocity',
            'units': 'degrees_per_second',
        }

    if 'iElV' in hdr_keys:
        d['elevation_velocity'] = {
            'data': degree_array('iElV'),
            'long_name': 'Elevation velocity',
            'units': 'degrees_per_second',
        }

    if 'iAz' in hdr_keys:
        d['azimuth'] = {
            'data': degree_array('iAz'),
            'long_name': 'Azimuth angle',
            'units': 'degrees',
        }

    if 'iEl' in hdr_keys:
        d['elevation'] = {
            'data': degree_array('iEl'),
            'long_name': 'Elevation angle',
            'units': 'degrees'
        }

    if 'iNumVecs' in hdr_keys:
        d['number_of_samples'] = {
            'data': int_array('iNumVecs'),
            'long_name': 'Samples in pulse',
        }

    if 'RX[0].fBurstMag' in hdr_keys:
        d['pulse_burst_magnitude_0'] = {
            'data': volt_array('RX[0].fBurstMag'),
            'long_name': 'Magnitude of first channel pulse burst',
            'units': 'volt',
            'comment': 'Amplitude not power',
        }

    if 'RX[0].iBurstArg' in hdr_keys:
        d['pulse_burst_phase_change_0'] = {
            'data': degree_array('RX[0].iBurstArg'),
            'long_name': 'Phase change from previous pulse',
            'units': 'degrees',
        }

    if 'RX[1].fBurstMag' in hdr_keys:
        d['pulse_burst_magnitude_1'] = {
            'data': volt_array('RX[1].fBurstMag'),
            'long_name': 'Magnitude of second channel pulse burst',
            'units': 'volt',
            'comment': 'Ampltide not power',
        }

    if 'RX[1].iBurstArg' in hdr_keys:
        d['pulse_burst_phase_change_1'] = {
            'data': degree_array('RX[1].iBurstArg'),
            'long_name': 'Phase change from previous pulse',
            'units': 'degrees'
        }

    return d


# taken from pyart/pyart/io/netcdf.py
def _create_ncvar(dic, dataset, name, dimensions):
    """
    Create and fill a Variable in a netCDF Dataset object.

    Parameters
    ----------
    dic : dict
        Radar dictionary to containing variable data and meta-data
    dataset : Dataset
        NetCDF dataset to create variable in.
    name : str
        Name of variable to create.
    dimension : tuple of str
        Dimension of variable.

    """
    # create array from list, etc.
    data = dic['data']
    if isinstance(data, np.ndarray) is not True:
        print "Warning, converting non-array to array:", name
        data = np.array(data)

    # convert string array to character arrays
    if data.dtype.char is 'S' and data.dtype != 'S1':
        import netCDF4
        data = netCDF4.stringtochar(data)

    # create the dataset variable
    if 'least_significant_digit' in dic:
        lsd = dic['least_significant_digit']
    else:
        lsd = None
    if "_FillValue" in dic:
        fill_value = dic['_FillValue']
    else:
        fill_value = None

    ncvar = dataset.createVariable(name, data.dtype, dimensions,
                                   zlib=True, least_significant_digit=lsd,
                                   fill_value=fill_value)

    # set all attributes
    for key, value in dic.iteritems():
        if key not in ['data', '_FillValue']:
            ncvar.setncattr(key, value)

    # set the data
    if data.shape == ():
        data.shape = (1,)
    ncvar[:] = data[:]
    #if type(data) == np.ma.MaskedArray:
    #    ncvar[:] = data.data
    #else:
    #    ncvar[:] = data
