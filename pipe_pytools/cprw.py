import numpy as _np
import forio as _io


def _read_vfield(f, Im, Jm, Km):
    vel = _np.zeros((3, Km+2, Jm+2, Im+2), order='C')
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            for d in range(3):
                _io.read_array(f, 'd', vel[d,k,j,1:-1])

    return vel


def _write_vfield(f, vel):
    Km = vel.shape[1] - 2
    Jm = vel.shape[2] - 2
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            for d in range(3):
                _io.write_array(f, 'd', vel[d,k,j,1:-1])
            
    return 


def get_scp(fname):
    f = open(fn, "rb")
    t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym = _io.read_list(f, "ddddddiiii")
    _read_vfield(f, 2**lx, Jm, 2**lt)
    f.close()
    return (t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym), vel

def get_cp(fname):
    f = open(fn, "rb")
    t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt = _io.read_list(f, "ddddddiii")
    _read_vfield(f, 2**lx, Jm, 2**lt)
    f.close()
    return (t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt), vel


def put_scp(fn, t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym, vel):
    f = open(fn,"wb")
    _io.write_list(f, "ddddddiiii", (t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym))
    _write_vfield(f, vel)
    f.close()    

def put_cp(fn, t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt, vel):
    f = open(fn,"wb")
    _io.write_list(f, "ddddddiii", (t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt))
    _write_vfield(f, vel)
    f.close()    


def put_car(cf, dtmax, tmax, kwrt=10000, kprt=10, tol=1.e-2, cpfn="tmp.cp", datfn="a0.dat"):
    f = open("pipe.car","w")
    f.write("%22.15f    -tol\n" % tol)
    f.write("%i         -kprt\n" % kprt)
    f.write("%i         -kwrt\n" % kwrt)
    f.write("%22.15f    -tmax\n" % tmax)
    f.write("%18.15f    -dtmax\n" % dt)
    f.write("%18.15f    -cf\n" % cf)
    f.write("%s\n" % cpfn)
    f.write("%s\n" % datfn)
    f.close()
    return

