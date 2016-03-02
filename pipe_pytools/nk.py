import elementary as _el
import subprocess as _sp
import cprw as _rw
import numpy as np

pos = {'dt':0, 'cf':1, 'Re':2, 'Xmax':3, 'epsr':4, 'lx':5, 'Jm':6, 'lt':7, 'nsym':8}
DT = (0,0,0,0)
CF = (0,0,0,1)
RE = (0,0,0,2)
XMAX = (0,0,0,3)
EPSR = (0,0,0,4)
LX = (0,0,0,5)
JM = (0,0,0,6)
LT = (0,0,0,7)
NSYM = (0,0,0,8)


def put_scp(fname, velb, time=0.0):
    dt,cf,Re,Xmax,epsr,lx,Jm,lt,nsym = velb[0,0,0,:9]
    _rw.put_scp(fname, time, dt, cf, Re, Xmax, epsr, lx, Jm, lt, nsym, velb)

def get_scp(fname, cf):
    (t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym), vel = get_scp(fname)
    vel[:9] = dt,cf,Re,Xmax,epsr,lx,Jm,lt,nsym
    return vel
	
def dot(v1, v2):
    if (v1[0,0,0,3:9] != v2[0,0,0,3:9]).any(): return 0.0/0.0
    Xmax,epsr,lx,Jm,lt,nsym = v1[0,0,0,3:9]
    return _el._pipe_mean((v1*v2).sum(0), Xmax,epsr,lx,Jm,lt,nsym)

def norm(v1): 
	return dot(v1, v1)**0.5

def calc(vel0, nstep):
    dt, cf, Re = vel0[0,0,0,:3]
    print "calc: dt=%20.15f, cf=%20.15f, Re=%20.15f,  Do ...     " % (dt,cf,Re), 
    
    _sp.call(['rm','-f','a0.dat'])
    _rw.put_car(cf, dt, dt*nstep, cpfn="tmp.cp", datfn="a0.dat")
    put_scp("tmp.scp", vel0)
    
    _sp.call(['mpirun','-np4','./pipe.out'])
    
    velT = get_scp('tmp.scp', cf)
    print "Done!"
    return velT

def F(vel): 
	return calc(vel, 200) - vel

def calc_Jdx(x0, res0, dx):
    eps = 1e-7
    return (res0 - F(x0 + eps*dx)) / eps


def _modify_KS(KS, x0, res0, q, res_lose, dx_accum):
    p = calc_Jdx(x0, res0, q)

    for q_old, p_old in KS:
        k = dot(p, p_old)
        p -= p_old * k
        q -= q_old * k

    p_norm = norm(p)
    q /= p_norm
    p /= p_norm

    k = dot(res_lose, p)
    res_lose -= k * p
    dx_accum += k * q
    KS.append((q, p))

    res_norm = norm(res_lose)
    print "Done, p_norm = %20.10e, effect = %20.10e, res_norm=%20.10e" % (p_norm, k, res_norm)

    return p.copy()



def solve_Jdx_is_r(x0, res0, add_param, tol=1e-7):
    print "GMRES start: res_norm = %20.10e" % norm(res0)
    
    KS = []
    dx_accum = np.zeros_like(x0)
    res_lose = res0.copy()
    
    
    for param in add_param:
        q = np.zeros_like(x0)
        q[(0,0,0,pos[param])] = 1.0
        print "#Extra %s " % param,
        _modify_KS(KS, x0, res0, q, res_lose, dx_accum)
    
    q = res_lose / norm(res_lose)
    while norm(res_lose) > tol:
        print "#Regular ", 
        q = _modify_KS(KS, x0, res0, q, res_lose, dx_accum)
    
    return dx_accum


def NewtonKrylow(x0, res0, add, tol):
    dx = solve_Jdx_is_r(x0, res0, add, tol)
    print "dx: norm=%20.15f, dt=%20.15f, cf=%20.15f, Re=%20.15f" % (norm(dx), dx[0,0,0,0], dx[0,0,0,1], dx[0,0,0,2])
    return x0 + dx


def extrapol2(vel0, vel1, vel2, param, p_new):
    p0 = vel0[(0,0,0,pos[param])]
    p1 = vel1[(0,0,0,pos[param])]
    p2 = vel2[(0,0,0,pos[param])]
    print [p0, p1, p2]
    
    X = np.matrix([[p0**2, p0, 1], [p1**2, p1, 1], [p2**2, p2, 1]])
    
    x4 = np.matrix([[p_new**2, p_new, 1]])
    
    c = x4 * X ** (-1)
    print "c=", c[0]
    
    vel_new = vel0*c[0,0] + vel1*c[0,1] + vel2*c[0,2]
    return vel_new


def extrapol1(vel0, vel1, param, p_new):
    p0 = vel0[(0,0,0,pos[param])]
    p1 = vel1[(0,0,0,pos[param])]
    print [p0, p1]
    
    X = np.matrix([[p0, 1], [p1, 1]])
    
    x3 = np.matrix([[p_new, 1]])
    
    c = x3 * X ** (-1)
    print "c=", c[0]
    
    vel_new = vel0*c[0,0] + vel1*c[0,1]
    return vel_new


