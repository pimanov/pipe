import numpy as np
import tools
import os

    
def chec_lam(vel, turb_vel):
    amp0 = np.abs(turb_vel[2,1:-1,1:-1,1:-1]).max()
    amp = np.abs(vel[2,1:-1,1:-1,1:-1]).max()
    if amp > 0.9*amp0: return -1
    if amp < 0.1*amp0: return 1
    return 0


def check_scp_lam_tend(sname, vel, al):
    
    al_vel = vel*al
    time = 0.0
    tmax = 0.0
    tend = 0
    while tend == 0:
        tmax += 200
        tools.put_scp("tmp.scp",time,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym, al_vel)
        tools.put_car(tmax, dt, cf, "tmp.scp")
        
        _ = os.execv('rm', ['-f', 'a0.dat'])
        _ = os.execv('mpirun', ['-np', '4', './pipe.out'])
        
        time,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel_al = tools.get_scp("tmp.scp")
        tend = chec_lam(vel_al, vel)
        
    a_name = sname + "_prt/a%s.dat" % str(n).zfill(3)
    os.execv('mv', ['a0.dat', a_name])
    
    return 0 if tend == -1 else 1
    
    
def init(sname, n=0, al=1.0, dal=1.0, lam_tend=0):
    f = open(sname,"a")
    f.write("%3d %25.22f %25.22f %d\n" % (n,al,dal,lam_tend))
    f.close()
    
    
def step(sname, vel, check_lam_tend):
        
    n, al, dal, lam_tend = np.loadtxt(sname+".txt", usecols=range(4), unpuck=False)[-1]
    print sname+".txt[-1]:", n, al, dal, lam_tend
    
    n += 1
    dal /= 2
    if lam_tend == 1: 
        al += dal
    else: 
        al -= dal
        
    lam_tend = check_lam_tend(sname, vel, al)
    
    f = open(sname,"a")
    f.write("%3d %25.22f %25.22f %d\n" % (n,al,dal,lam_tend))
    f.close()
