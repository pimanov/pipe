import struct as _st


_dtsize = {'i': 4, 'd': 8}
def _read(f, dtype): return _st.unpack(dtype, f.read(_dtsize[dtype]))[0]
def _write(f, dtype, value): f.write(_st.pack(dtype, value))


def read_list(f, fmt):
    _read(f, 'i')

    new_list = []
    for dtype in fmt: 
        new_list.append(_read(f, dtype))

    _read(f, 'i')
    return new_list


def read_array(f, dtype, array):
    _read(f, 'i')

    for i in range(len(array)): 
        array[i] = _read(f, dtype)

    _read(f, 'i')
    return array


def write_list(f, fmt, date):
    mem_size = 0
    for dtype in fmt: 
        mem_size += _dtsize[dtype]

    _write(f, 'i', mem_size)

    for dtype, value in zip(fmt, date): 
        _write(f, dtype, value)

    _write(f, 'i', mem_size)


def write_array(f, dtype, data):
    n = len(data)
    mem_size = n * _dtsize[dtype]

    _write(f, 'i', mem_size)

    for value in data: 
        _write(f, dtype, value)

    _write(f, 'i', mem_size)

