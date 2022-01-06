import os

import numpy as np

"Functions for decoding and reading CAEN .dat binary files"

# eventSize = (header[0] - 24) //2
# Total data size per wave = 4120
# Header size = 6 * uint32 = 24
# Event size = 1024 * float32(??) = 4096
# CAEN say samples are two bytes in length but apparently they are 4 bytes???


# def def_dtype(hsize=6, rsize=1024, hdtype=np.int32, wdtype=np.float32):
#     dt = np.dtype([('header', np.int32, 6),
#                    ('data', np.float32, 1024)])
#     return dt

def def_dtype(type, size):
    dt = np.dtype((type, (size,)))
    return dt


def read_one(f, hsize=6, rsize=1024, hdtype=np.int32, wdtype=np.float32):
    # Reads one waveform from the last position the file was read.
    # Should be used inside with open(...) statement.
    # For V1742, use wdtype = np.float32
    # For V1730B, use wdtype = np.int16 (not tested)
    # If digitiser has different header or record length, change hsize or rsize
    # Size is the number of that datatype to read, not bytes
    if f.tell() >= os.stat(f.name).st_size:
        return
    header = np.fromfile(f, dtype=hdtype, count=hsize)
    wave = np.fromfile(f, dtype=wdtype, count=rsize)
    if header is None or wave is None:
        return
    return header, wave


def read_one_new(f, hdtype, rdtype):
    # Define header, wave combined dtype using def_dtype() and pass into here
    header = np.fromfile(f, dtype=hdtype, count=1)
    wave = np.fromfile(f, dtype=rdtype, count=1)
    return header, wave


def read_dat(fname, hsize=6, rsize=1024, hdtype=np.int32, wdtype=np.float32):
    # Generator to read entire file sequentially. Use as iterable in for loop.
    # i.e. for header, wave in read_dat(...): do_stuff
    # Same notes about data types from read_one apply here.
    read_func = read_one
    with open(fname, "rb") as f:
        while f.tell() < os.stat(f.name).st_size:
            header, wave = read_func(f, hsize, rsize, hdtype, wdtype)
            yield header, wave


def find_wave(f, bytes_to_seek, hsize=6, rsize=1024, hdtype=np.int32, wdtype=np.float32):
    # header_bytes = hdtype(0).nbytes
    # wave_bytes = wdtype(0).nbytes
    # bytes_to_seek = n_ev * (header_bytes + wave_bytes)
    f.seek(bytes_to_seek)
    header, wave = read_one(f, hsize, rsize, hdtype, wdtype)
    return header, wave


def find_dats(dir):
    datfiles = []
    for file in os.listdir(dir):
        if file.endswith(".dat"):
            datfiles.append(os.path.join(dir, file))
    return datfiles


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i: i + n]


class CAENReader():

    def __init__(self, filename, header_bytes, header_dtype, record_bytes, record_dtype) -> None:
        pass
