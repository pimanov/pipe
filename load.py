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

    u = nans((Im+2, Jm+2, Km+2))
    v = nans((Im+2, Jm+2, Km+2))
    w = nans((Im+2, Jm+2, Km+2))

    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            u[1:-1,j,k] = fget_dline(f,Im)
            v[1:-1,j,k] = fget_dline(f,Im)
            w[1:-1,j,k] = fget_dline(f,Im)

    return u,v,w



def fget_2Dvfield(f,Im,Jm):

    u = nans((Im+2, Jm+2))
    v = nans((Im+2, Jm+2))
    w = nans((Im+2, Jm+2))

    for j in range(1,Jm+1):
        u[1:-1,j] = fget_dline(f,Im)
        v[1:-1,j] = fget_dline(f,Im)
        w[1:-1,j] = fget_dline(f,Im)

    return u,v,w
