import numpy as np
#cimport numpy as np

def f(x):
    return x**2-x

def integrate_f(a, b, N):
    s = 0
    dx = (b-a)/N
    for i in range(N):
        s += f(a+i*dx)
    return s * dx


def decode_highsnr(highsnr_data):
    """ Decode a pulsedata vector encoded in High SNR format. """
    decoded_data = np.empty(highsnr_data.shape, dtype=np.float32)
    for i, p in enumerate(highsnr_data):
        decoded_data[i] = decode_highsnr_point(p)
    return decoded_data

def decode_highsnr_point(unsigned short p, debug=False):
    """ Decode a point encoded in High SNR format. """

    # find the sign_bit, exponent, and mantissa
    sign_bit = p >> 11 & 1
    exponent = p >> 12
    mantissa = p & 2047

    #print "exponent:", exponent
    #print "sign_bit:", sign_bit
    #print "mantissa:", mantissa

    if exponent != 0:
        if sign_bit:  # sign == 1
            val = -4096 + mantissa
        else:  # sign == 0
            val = 2048 + mantissa
        return val * 2. ** (exponent - 25)
    else:  # exponent == 0
        if sign_bit:
            # two compliment
            val = -2048 + mantissa
        else:
            val = mantissa
        return val * 2. ** -24
