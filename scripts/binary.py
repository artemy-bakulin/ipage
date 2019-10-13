import numpy as np
import hashlib
def write_np_array(inp_array, out_filename):
    shape_n_dimensions = len(inp_array.shape)
    n_dimentions_bitstring = np.uint8(shape_n_dimensions).tobytes()
    shape_array_bitstring = np.array(inp_array.shape, dtype=np.uint32).tobytes()

    array_bitstring = inp_array.tobytes()
    md5 = hashlib.md5()
    md5.update(array_bitstring)
    md5_checksum = md5.digest()
    full_bytestring = n_dimentions_bitstring + shape_array_bitstring + array_bitstring + md5_checksum

    with open(out_filename, 'wb') as wf:
        wf.write(full_bytestring)


def read_np_array(inp_filename, dtype):
    with open(inp_filename, 'rb') as rf:
        bitstring = rf.read()

    n_dimentions_bitstring = bitstring[0: 1]
    n_dimentions = np.frombuffer(n_dimentions_bitstring, dtype=np.uint8)[0]
    shape_array_bitstring = bitstring[1: 1 + 4 * n_dimentions]
    shape_array = np.frombuffer(shape_array_bitstring, dtype=np.uint32)
    flatten_length = int(np.prod(shape_array))

    output_bitstring = bitstring[1 + 4 * n_dimentions :
                                1 + 4 * n_dimentions + dtype.itemsize * flatten_length]
    md5_checksum = bitstring[1 + 4 * n_dimentions + dtype.itemsize * flatten_length : ]
    output_array = np.frombuffer(output_bitstring, dtype=dtype)
    reshaped_array = np.reshape(output_array, shape_array, order='C')

    array_bitstring = reshaped_array.tobytes()
    md5 = hashlib.md5()
    md5.update(array_bitstring)
    md5_read = md5.digest()
    assert(md5_checksum == md5_read)
    return reshaped_array