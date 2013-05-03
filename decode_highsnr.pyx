# cython: profile=True
import numpy as np


def decode_highsnr(unsigned short [:] highsnr_data, int npoints):
    """ Decode a pulsedata vector encoded in High SNR format. """
    cdef float [:] decoded_data = np.empty((npoints, ), dtype=np.float32) 
    cdef unsigned short p
    for i in range(npoints):
        p = highsnr_data[i]
        decoded_data[i] = decode_highsnr_point(p)
    return decoded_data


cpdef decode_highsnr_point(unsigned short p, debug=False):
    """ Decode a point encoded in High SNR format. """

    # find the sign_bit, exponent, and mantissa
    cdef unsigned short sign_bit = p >> 11 & 1
    cdef unsigned short exponent = p >> 12
    cdef unsigned short mantissa = p & 2047
    cdef short val
    cdef float r

    if exponent != 0:
        if sign_bit:  # sign == 1
            val = -4096 + mantissa
        else:  # sign == 0
            val = 2048 + mantissa
        r = val * 2.**(exponent -25)
        return r
    else:  # exponent == 0
        if sign_bit:
            # two compliment
            val = -2048 + mantissa
        else:
            val = mantissa
        r = val * 2. ** -24
        return r
