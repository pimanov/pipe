import numpy as np
import struct


def fget_d(f):
    d = f.read(8)
    res, = struct.unpack('d',d)
    return res



def fget_i(f):
    i = f.read(4)
    res, = struct.unpack('i',i)
    return res



def fget_dline(f,n):
    
    dline = np.zeros(n)
    
    fget_i(f)
    for i in range(n): dline[i] = fget_d(f)
    fget_i(f)

    return dline



def nans(a):
    res = np.zeros(a)
    res[:] = float('nan')
    return res



def fget_vfield(f,Im,Jm,Km):

    u = nans((Km+2, Jm+2, Im+2))
    v = nans((Km+2, Jm+2, Im+2))
    w = nans((Km+2, Jm+2, Im+2))

    for u1,v1,w1 in zip(u,v,w)[1:-1]:
        for u2,v2,w2 in zip(u1,v1,w1)[1:-1]:
            u2[1:-1] = fget_dline(f,Im)
            v2[1:-1] = fget_dline(f,Im)
            w2[1:-1] = fget_dline(f,Im)

    return u,v,w



def fget_2Dvfield(f,Im,Jm):

    u = nans((Jm+2, Im+2))
    v = nans((Jm+2, Im+2))
    w = nans((Jm+2, Im+2))

    for u1,v1,w1 in zip(u,v,w)[1:-1]:
        u1[1:-1] = fget_dline(f,Im)
        v1[1:-1] = fget_dline(f,Im)
        w1[1:-1] = fget_dline(f,Im)

    return u,v,w



def fget_2Dfield(f,Im,Jm):

    u = nans((Jm+2, Im+2))

    for u1 in u[1:-1]:
        u1[1:-1] = fget_dline(f,Im)

    return u
