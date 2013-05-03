import numpy as np

import iq_reader

t = iq_reader.read_iq_ascii('sample_data/example_1')

data = t[-1]
np.save('data.npy', data)
