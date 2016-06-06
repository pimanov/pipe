import struct
import os
import numpy as np


_get_i = lambda(f): struct.unpack('i', f.read(4))[0]
_get_d = lambda(f): struct.unpack('d', f.read(8))[0]

def get_line(f,line):
    n = line.shape[0]
    
    n1 = _get_i(f)
    if (n1 != n*8): print "tools.py: get_line(): warning! fortran binary file is incorrect"
    
    for i in range(n): 
        line[i] = _get_d(f)
        
    n1 = _get_i(f)
    if (n1 != n*8): print "tools.py: get_line(): warning! fortran binary file is incorrect"


def get_vfield(f,Im,Jm,Km):
    vel = np.zeros((3, Km+2, Jm+2, Im+2), order='C')
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            get_line(f, vel[0,k,j,1:-1])
            get_line(f, vel[1,k,j,1:-1])
            get_line(f, vel[2,k,j,1:-1])

    return vel

def get_scp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,_ = i(f),d(f),d(f),d(f),d(f),d(f),d(f),i(f),i(f),i(f),i(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel



def get_dcp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,_ = i(f),d(f),d(f),d(f),d(f),d(f),d(f),i(f),i(f),i(f),d(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel


def get_cp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,_ = i(f),d(f),d(f),d(f),d(f),d(f),d(f),i(f),i(f),i(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,vel


def get_bcp(fn):
    i = _get_i
    d = _get_d

    with open(fn, "rb") as f:
        _,cf,Xmax,epsr,lx,Jm,lt,nsym,_ = i(f),d(f),d(f),d(f),i(f),i(f),i(f),i(f),i(f)
        vel = get_vfield(f, 2**lx, Jm, 2**lt)

    return cf,Xmax,epsr,lx,Jm,lt,nsym,vel


def put_line(f,line):
    n = line.shape[0]
    f.write(struct.pack("i",n*8))
    for l in line: 
        f.write(struct.pack("d",l))
    f.write(struct.pack("i",n*8))
    return 


def put_vfield(f,vel):
    Km = vel.shape[1] - 2
    Jm = vel.shape[2] - 2
    
    for k in range(1,Km+1):
        for j in range(1,Jm+1):
            put_line(f,vel[0,k,j,1:-1])
            put_line(f,vel[1,k,j,1:-1])
            put_line(f,vel[2,k,j,1:-1])
            
    return 

def put_scp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym):
    m = 4*4+6*8
    f.write(struct.pack("i",m))
    f.write(struct.pack("d",t))
    f.write(struct.pack("d",dt))
    f.write(struct.pack("d",Dp))
    f.write(struct.pack("d",Re))
    f.write(struct.pack("d",Xmax))
    f.write(struct.pack("d",epsr))
    f.write(struct.pack("i",lx))
    f.write(struct.pack("i",Jm))
    f.write(struct.pack("i",lt))
    f.write(struct.pack("i",nsym))
    f.write(struct.pack("i",m))
    return

def put_scp(fn,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel,comm="comment"):
    
    with open(fn,"wb") as f:
        put_scp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym)
        put_vfield(f,vel)
        f.write(comm)
        f.close()
    
    return


def put_bcp_header(f,cf,Xmax,epsr,lx,Jm,lt,nsym):
    m = 4*4+3*8
    f.write(struct.pack("i",m))
    f.write(struct.pack("d",cf))
    f.write(struct.pack("d",Xmax))
    f.write(struct.pack("d",epsr))
    f.write(struct.pack("i",lx))
    f.write(struct.pack("i",Jm))
    f.write(struct.pack("i",lt))
    f.write(struct.pack("i",nsym))
    f.write(struct.pack("i",m))
    return

def put_bcp(fn,cf,Xmax,epsr,lx,Jm,lt,nsym,vel,comm="comment"):
    
    with open(fn,"wb") as f:
        put_bcp_header(f,cf,Xmax,epsr,lx,Jm,lt,nsym)
        put_vfield(f,vel)
        f.write(comm)
        f.close()
    
    return


def put_dcp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym):
    m = 3*4+7*8
    f.write(struct.pack("i",m))
    f.write(struct.pack("d",t))
    f.write(struct.pack("d",dt))
    f.write(struct.pack("d",Dp))
    f.write(struct.pack("d",Re))
    f.write(struct.pack("d",Xmax))
    f.write(struct.pack("d",epsr))
    f.write(struct.pack("i",lx))
    f.write(struct.pack("i",Jm))
    f.write(struct.pack("i",lt))
    f.write(struct.pack("d",nsym))
    f.write(struct.pack("i",m))
    return

def put_dcp(fn,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel):
    
    with open(fn,"wb") as f:
        put_dcp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym)
        put_vfield(f,vel)
    
    return


def put_cp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt):
    m = 3*4+6*8
    f.write(struct.pack("i",m))
    f.write(struct.pack("d",t))
    f.write(struct.pack("d",dt))
    f.write(struct.pack("d",Dp))
    f.write(struct.pack("d",Re))
    f.write(struct.pack("d",Xmax))
    f.write(struct.pack("d",epsr))
    f.write(struct.pack("i",lx))
    f.write(struct.pack("i",Jm))
    f.write(struct.pack("i",lt))
    f.write(struct.pack("i",m))
    return

def put_cp(fn,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,vel):
    
    with open(fn,"wb") as f:
        put_cp_header(f,t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt)
        put_vfield(f,vel)
    
    return


def put_car(tmax, dtmax, cf, cpfn, tol=0.01, kprt=10, kwrt=100000, prtfn="a0.dat"):
    f = open("pipe.car","w")
    f.write("%f - tol\n" % tol)
    f.write("%d - kprt\n" % kprt)
    f.write("%d - kwrt\n" % kwrt)
    f.write("%22.15f - tmax\n" % tmax)
    f.write("%18.15f - dtmax\n" % dtmax)
    f.write("%18.15f - cf\n" % cf)
    f.write("%s\n" % cpfn)
    f.write("%s" % prtfn)
    f.close()
    return




def put_bcar(tmax, dtmax, bcpfn, cpfn, tol=0.01, kprt=10, kwrt=100000, prtfn="a0.dat"):
    f = open("pipe.car","w")
    f.write("%f - tol\n" % tol)
    f.write("%d - kprt\n" % kprt)
    f.write("%d - kwrt\n" % kwrt)
    f.write("%22.15f - tmax\n" % tmax)
    f.write("%18.15f - dtmax\n" % dtmax)
    f.write("%s\n" % bcpfn)
    f.write("%s\n" % cpfn)
    f.write("%s" % prtfn)
    f.close()
    return


