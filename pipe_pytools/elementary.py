import numpy as _np
import grid as _grid

def get_ucl(u):
    return u[1:-1,1,:].mean(0)


def get_uPois((X,R,Th)):
    u = np.zeros((Th.m+2, R.m+2, X.m+2), order='C')
    for k in xrange(1, Th.m+1):
        for j in xrange(1, R.m+1):
            for i in xrange(1, X.m+1):
                u[k,j,i] = 1.0 - R.f[j]**2
    return u


def get_vPois((X,R,Th)):
    vel = np.zeros((3, Th.m+2, R.m+2, X.m+2), order='C')
    for k in xrange(1, Th.m+1):
        for j in xrange(1, R.m+1):
            for i in xrange(1, X.m+1):
                vel[0,k,j,i] = 1.0 - R.f[j]**2
    return vel


def _get_uPois(lx, Xmax, Jm, epsr, lt, nsym):
    X = _grid.X(2**lx, Xmax)
    R = _grid.R(Jm, epsr)
    Th = _grid.Th(2**lt, nsym)
    return get_uPois((X,R,Th))
    

def _get_vPois(lx, Xmax, Jm, epsr, lt, nsym):
    X = _grid.X(2**lx, Xmax)
    R = _grid.R(Jm, epsr)
    Th = _grid.Th(2**lt, nsym)
    return get_vPois((X,R,Th))


def cs_mean(u,(R,Th)):
    res = 0.0
    ss = 0.0
    for k in xrange(1, Th.m+1):
        for j in xrange(1, R.m+1):
            res += R.f1[j] * R.f[j] * Th.h * u[k,j]
            ss  += R.f1[j] * R.f[j] * Th.h
    return res / ss


def pipe_mean(u, (X,R,Th)):
    res = 0.0
    ll = 0.0
    for i in xrange(1, X.m+1):
        res += cs_mean_(u[:,:,i], (R,Th)) * X.h
        ll  += X.h
    return res / ll


def _cs_mean(u, Jm, epsr, lt, nsym):
    R = _grid.R(Jm, epsr)
    Th = _grid.Th(2**lt, nsym)
    return cs_mean(u, (R,Th))


def _pipe_mean(u, lx, Xmax, Jm, epsr, lt, nsym):
    X = _grid.X(2**lx, Xmax)
    R = _grid.R(Jm, epsr)
    Th = _grid.Th(2**lt, nsym)
    return pipe_mean(u, (X,R,Th))
    

