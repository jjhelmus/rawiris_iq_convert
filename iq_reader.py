import os

import numpy as np


def read_iq_ascii(filename):

    length_of_file = os.stat(filename).st_size

    f = open(filename, 'r')

    pulsehdr_locations = []
    pulsehdr_dics = []
    pulsedata_locations = []
    n_points = []
    f.seek(0)
    pulse_info = read_ascii_rvptsPulseInfo(f)

    while f.tell() != length_of_file:

        # record the location of the pulse header
        pulsehdr_locations.append(f.tell())

        # read the pulse header
        current_pulsehdr = read_ascii_rvptsPulseHdr(f)
        pulsehdr_dics.append(current_pulsehdr)

        # record the location of the pulse data
        pulsedata_locations.append(f.tell())

        # skip over the pulse data
        iNumVecs = int(current_pulsehdr['iNumVecs'])
        iVIQPerBin = int(current_pulsehdr['iVIQPerBin'])
        n_points.append(iNumVecs * iVIQPerBin)
        bytes_in_current_pulse = 2 * 2 * iNumVecs * iVIQPerBin
        f.seek(bytes_in_current_pulse, 1)

    n_pulses = len(n_points)
    max_points = max(n_points)

    data = np.empty((n_pulses, max_points * 2), dtype=np.float32)
    for i, (npoints) in enumerate(n_points):
        #npoints = n_points[i]
        pulse_start = pulsedata_locations[i]
        f.seek(pulse_start)
        data[i] = read_pulsedata(f, npoints)

    return pulse_info, pulsehdr_dics, data


def read_ascii_rvptsPulseInfo(f):
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


def read_ascii_rvptsPulseHdr(f):
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


def read_pulsedata(f, npoints):
    """ npoints is number I|Q data (numvecs * VIQperbin) """

    buf = f.read(npoints * 4)
    high_snr_data = np.fromstring(buf, np.uint16)

    decoded_data = np.empty(high_snr_data.shape, dtype=np.float32)
    for i, p in enumerate(high_snr_data):
        decoded_data[i] = decode_high_snr_point(p)
    return decoded_data


def decode_high_snr_point(p, debug=False):

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
