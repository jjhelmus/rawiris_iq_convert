"""
Routines for reading IQ data from IRIS TS files (RAW IRIS files).

"""

import os

import numpy as np

import decode_highsnr


class IrisTSFile:
    """
    A class for working with IRIS TS files which hold IQ data

    Parameters
    ----------
    filename : str
        Filename of IRIS TS file to read data from.

    Attributes
    ----------
    filename : str
        Filename of IRIS TS file being read.
    npulse : int
        Number of pulse in the file.
    pulse_info : dic
        Dictionary of

    Private attributes
    ------------------
    _f : filehandler

    _data_locations : array

    _npoints : array

    _maxpoints : int

    _hdr_dics : list

    """

    def __init__(self, filename):
        """ initalize. """

        length_of_file = os.stat(filename).st_size

        f = open(filename, 'r')

        self._hdr_dics = []
        data_locations = []
        npoints = []

        # read in the pulse info, store as pulse_info attribute
        self.pulse_info = self._get_rvptspulseinfo_dic(f)

        # read through the file recording the data locations and pulse
        # header dictionaries
        while f.tell() != length_of_file:

            # read and store the pulse header
            current_pulsehdr = self._get_rvptspulsehdr_dic(f)
            self._hdr_dics.append(current_pulsehdr)

            # record the location of the pulse data
            data_locations.append(f.tell())

            # record the number of points in the pulse data
            numvecs = int(current_pulsehdr['iNumVecs'])
            viqperbin = int(current_pulsehdr['iVIQPerBin'])
            npoints.append(numvecs * viqperbin * 2)

            # skip over the pulse data
            bytes_in_current_pulse = 2 * 2 * numvecs * viqperbin
            f.seek(bytes_in_current_pulse, 1)

        # set attributes
        self._f = f
        self.filename = filename
        self._data_locations = np.array(data_locations, dtype='int64')
        self._npoints = np.array(npoints, dtype='int32')
        self.npulses = len(self._npoints)
        self._max_points = np.max(self._npoints)

        return

    def _get_rvptspulseinfo_dic(self, f):
        """ Read an ASCII rvptspulseinfo header. """
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
        """ Read an ASCII rvptspulsehdr header. """
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
        Read and decode pulse data
        """
        if pulse_id < 0 or pulse_id >= self.npulses:
            raise ValueError('invalid pulse_id')

        # determine number of points in the pulse
        n_points = self._npoints[pulse_id]

        # read pulse data from file
        data_location = self._data_locations[pulse_id]
        self._f.seek(data_location)
        buf = self._f.read(n_points * 2)
        highsnr_data = np.fromstring(buf, np.uint16)

        # return decoded data
        return decode_highsnr.decode_highsnr(highsnr_data, n_points)

    def read_all_pulse_data(self):
        """ Read and return all data from all pulses. """
        data = np.empty((self.npulses, self._max_points), dtype=np.float32)
        for i, npoints in enumerate(self._npoints):
            data[i, :npoints] = self.read_pulse_data(i)
        return data

    def _decode_highsnr_python(self, highsnr_data, npoints):
        """ Decode a pulsedata vector encoded in High SNR format. """
        # a much faster implementatio of this function is available in
        # Cython module decode_highsnr.pyx.  This is kept for
        # reference
        decoded_data = np.empty((npoints), dtype=np.float32)
        for i, p in enumerate(highsnr_data):
            decoded_data[i] = self._decode_highsnr_point_python(p)
        return decoded_data

    def _decode_highsnr_point_python(self, p, debug=False):
        """ Decode a point encoded in High SNR format. """
        # a much faster implementatio of this function is available in
        # Cython module decode_highsnr.pyx.  This is kept for
        # reference

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
    Write a IrisTSFile object to netCDF file.

    Parameters
    ----------
    filename : str
        Name of file to write
    iristsfile : IrisTSFile
        IrisTSFile object to write to netCDF file.
    format : str


    """
    import netCDF4

    # create the dataset
    dataset = netCDF4.Dataset(filename, 'w', format=format)

    npulses = iristsfile.npulses
    nsamples = iristsfile._max_points / 2
    pulse_info = iristsfile.pulse_info
    if pulse_info['iPolarization'] == '3':
        dual_polarized = True
        nsamples /= 2
    else:
        dual_polarized = False

    # create the pulse and sample dimensions
    dataset.createDimension('pulse', npulses)
    dataset.createDimension('sample', nsamples)

    # create the I and Q variables
    i0 = dataset.createVariable('I_0', 'f4', ('pulse', 'sample'))
    q0 = dataset.createVariable('Q_0', 'f4', ('pulse', 'sample'))
    if dual_polarized:
        i1 = dataset.createVariable('I_1', 'f4', ('pulse', 'sample'))
        q1 = dataset.createVariable('Q_1', 'f4', ('pulse', 'sample'))

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

    # Find and add the azimuth and elevation data
    if 'iEl' in iristsfile._hdr_dics[0]:
        iel = np.array([int(d['iEl']) for d in iristsfile._hdr_dics],
                       dtype='float32')
        elev = dataset.createVariable('elevation', 'f4', ('pulse', ))
        elev[:] = iel / 32768. * 180.
        elev.units = 'degrees'
        elev.comment = 'Elevation'

    if 'iAz' in iristsfile._hdr_dics[0]:
        iaz = np.array([int(d['iAz']) for d in iristsfile._hdr_dics],
                       dtype='float32')
        azim = dataset.createVariable('azimuth', 'f4', ('pulse', ))
        azim[:] = iaz / 32768. * 180.
        azim.units = 'degrees'
        azim.comment = 'Azimuth'

    dataset.close()
