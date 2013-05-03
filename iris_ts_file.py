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
    filename

    n_pulse : int
        Number of pulse in the file.
    pulse_info : dic
        Dictionary of
    _hdr_locations : list




    """

    def __init__(self, filename):
        """ initalize. """
        self.filename = filename

        length_of_file = os.stat(filename).st_size

        f = open(filename, 'r')

        self._hdr_locations = []
        self._hdr_dics = []
        self._data_locations = []
        self._numvecs = []
        self._viqperbin = []
        self._npoints = []

        # read in the pulse info
        self.pulse_info = self._get_rvptspulseinfo_dic(f)

        # read through the file recording the pulse headers
        # and data locations
        while f.tell() != length_of_file:

            # record the location of the pulse header
            self._hdr_locations.append(f.tell())

            # read the pulse header
            current_pulsehdr = self._get_rvptspulsehdr_dic(f)
            self._hdr_dics.append(current_pulsehdr)

            # record the location of the pulse data
            self._data_locations.append(f.tell())

            # record size of pulse data
            numvecs = int(current_pulsehdr['iNumVecs'])
            viqperbin = int(current_pulsehdr['iVIQPerBin'])
            self._numvecs.append(numvecs)
            self._viqperbin.append(viqperbin)
            self._npoints.append(numvecs * viqperbin * 2)

            # skip over the pulse data
            bytes_in_current_pulse = 2 * 2 * numvecs * viqperbin
            f.seek(bytes_in_current_pulse, 1)

        self.f = f
        self.n_pulses = len(self._numvecs)
        self._max_points = max(self._npoints)
        return

    def read_pulse_data(self, pulse_id):
        """
        Read and decode pulse data
        """
        if pulse_id < 0 or pulse_id >= self.n_pulses:
            raise ValueError('invalid pulse_id')

        # determine number of points in the pulse
        numvecs = self._numvecs[pulse_id]
        viqperbin = self._viqperbin[pulse_id]
        n_points = numvecs * viqperbin * 2

        # read pulse data from file
        data_location = self._data_locations[pulse_id]
        self.f.seek(data_location)
        buf = self.f.read(n_points * 2)
        highsnr_data = np.fromstring(buf, np.uint16)

        # return decoded data
        return decode_highsnr.decode_highsnr(highsnr_data)

    def read_all_pulse_data(self):
        """ Read and return all data from all pulses. """
        data = np.empty((self.n_pulses, self._max_points), dtype=np.float32)
        for i, npoints in enumerate(self._npoints):
            data[i, :npoints] = self.read_pulse_data(i)
        return data

    def _decode_highsnr_python(self, highsnr_data):
        """ Decode a pulsedata vector encoded in High SNR format. """
        decoded_data = np.empty(highsnr_data.shape, dtype=np.float32)
        for i, p in enumerate(highsnr_data):
            decoded_data[i] = self._decode_highsnr_point_python(p)
        return decoded_data

    def _decode_highsnr_point_python(self, p, debug=False):
        """ Decode a point encoded in High SNR format. """

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
