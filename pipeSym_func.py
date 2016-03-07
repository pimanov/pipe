
# coding: utf-8

# In[2]:

import pipeSym_wrapper as wr
import numpy as _np
import math as _math


def look():
    print "Xmax=%f, epsr=%f, nsym=%d" % (wr.dim.Xmax, wr.dim.epsr, wr.dim.nsym)
    print "Im=%d, Jm=%d, Km=%d" % (wr.dimx.Im, wr.dimr.Jm, wr.dimt.Km)
    print "hx=%f, ht=%f, real_nsym=%f" % (wr.dimx.hx, wr.dimt.ht, wr.dimt.Km * wr.dimt.ht)
    print "Re=%f" % (wr.Re.Re)

# In[13]:

def new_z_vfield(Im = None, Jm = None, Km = None):
    if(Im == None): Im = wr.dimx.Im
    if(Jm == None): Jm = wr.dimr.Jm
    if(Km == None): Km = wr.dimt.Km
    return _np.zeros((3,Km+2, Jm+2, Im+2), order='C')


# In[12]:

def new_z_pfield(Im = None, Jm = None, Km = None):
    if(Im == None): Im = wr.dimx.Im
    if(Jm == None): Jm = wr.dimr.Jm
    if(Km == None): Km = wr.dimt.Km
    return _np.zeros((Km+2, Jm+2, Im+2), order='C')


# In[4]:

def update_ht(ht):
    wr.dimt.ht = ht
    Km = wr.dimt.Km
        
    for k in range(0,Km):
        wr.rlt.rlt[k] = (2.0 / ht * _math.sin(0.5 * k * _math.pi / Km))**2


# In[6]:

def update_hx(hx):
    wr.dimx.hx = hx
    Im = wr.dimx.Im
    wr.dim.Xmax = wr.dimx.hx * Im

    for i in xrange(0, Im/2 +1):
        wr.rlx.rlx[i] = (2.0 / hx * _math.sin(i * _math.pi / Im))**2

    for i in range(1, Im/2):
        i1 = Im - i
        wr.rlx.rlx[i1] = wr.rlx.rlx[i]
    


# In[7]:

def divmax(vel):
    return wr.divmax(vel[0].T, vel[1].T, vel[2].T)


# In[8]:

def ffmean(u):
    return wr.mean(u.T)


# In[18]:

def init(Xmax=None, epsr=None, nsym=None, lx=None, Jm=None, lt=None, Im=None, Km=None, 
         vel=None, hx=None, ht=None, Re=None):
    
    if(Re != None): wr.Re.Re = Re
    if(Xmax != None): wr.dim.Xmax = Xmax
    if(epsr != None): wr.dim.epsr = epsr
        
    if(nsym != None):
        if(nsym % 1.0 == 0.0):
            wr.dim.nsym = nsym
            calc_ht = False
        else:
            wr.dim.nsym = 2.0
            calc_ht = True
        
    if(vel != None): 
        Im = vel.shape[1] - 2
        Jm = vel.shape[2] - 2
        Km = vel.shape[3] - 2
        
    if(Im != None): 
        lx = _math.log(Im,2)
        if(2**lx != Im):
            print "###Error! Im != 2**n"
            return
            
    if(Km != None): 
        lt = _math.log(Km,2)
        if(2**lt != Km):
            print "###Error! Km != 2**n"
            return
        
    if(lx != None): wr.dimx.lx = lx
    if(Jm != None): wr.dimr.Jm = Jm
    if(lt != None): wr.dimt.lt = lt
        
    wr.com()
    
    if(calc_ht == True): ht = _math.pi / (nsym * wr.dimt.Km)
    
    if(ht != None): update_ht(ht)
    if(hx != None): update_hx(hx)
        
    return


# In[9]:

def step(t, dt, vel, velt, om, p, cf, tol=0.01):
    t, dt = wr.step(t, dt, tol, vel[0].T, vel[1].T, vel[2].T, velt[0].T, velt[1].T, velt[2].T, 
                    om[0].T, om[1].T, om[2].T, p.T, cf)
    
    return t, dt, vel, velt, om, p


# In[10]:

def rp(t, vel, velt, om):
    wr.rp(t, vel[0].T, vel[1].T, vel[2].T, velt[0].T, velt[1].T, velt[2].T, om[0].T, om[1].T, om[2].T)
    return velt, om


# In[11]:

def pres(vel, p, ub):
    wr.pres(vel[0].T, vel[1].T, vel[2].T, p.T, ub)
    return vel, p


# In[11]:

def calc(vel, dtmax, cf=0.0, Re=None, velt=None, om=None, p=None, maxnstep=None, tmax=None, ht=None, hx=None, 
         tol=0.01, time=0.0, dt=None, prt=None, ampmax=None, ampmin=None, const_dt_chec=False):    
        
    if(hx != None): update_hx(hx)
    if(ht != None): update_ht(ht)
    if(Re != None): wr.Re.Re = Re
    if(dt == None): dt = dtmax

    if(om == None):
        om = new_z_vfield()
        
    if(velt == None): 
        velt = new_z_vfield()
        velt,om = rp(time, vel, velt, om)
        
    if(p == None):
        p = new_z_pfield()
        velt,p = pres(velt, p, 0.0)
        
    nstep = 0
    break_cond_chec = False
    while(True):
        if(tmax != None):
            break_cond_chec = True
            if(time > tmax + 0.1*dt): 
                print "Planned break! tmax reached"
                break
    
        time, dt1, vel, velt, om, p = step(time, dt, vel, velt, om, p, cf, tol)
        
        if(const_dt_chec == True):
            if(dt1 != dt): 
                print "###Forced break! dt has changed"
                break
        else: 
            dt = min(dt1, dtmax)
        
        if(prt != None): 
            amp = prt(time, dt, vel, velt, om, p, cf)
            
            if(ampmax != None): 
                break_cond_chec = True
                if(amp > ampmax): 
                    print "Planned break! ampmax reached"
                    break
            
            if(ampmin != None): 
                break_cond_chec = True
                if(amp < ampmin): 
                    print "Planned break! ampmin reached"
                    break
                
        nstep += 1
        if(maxnstep != None):
            break_cond_chec = True
            if (nstep >= maxnstep): 
                print "Planned break! maxnstep reached"
                break
                
                
        if(break_cond_chec == False):
            print "No break condition set. Forced break!"
            break
    
    return time, dt, vel, velt, om, p


# In[14]:

def read_vel(fname):
    vel = new_z_vfield()
    t, dt, Dp = wr.load(fname, vel[0].T, vel[1].T, vel[2].T)
    Re = wr.Re.Re
    return t, dt, Dp, Re, vel


# In[60]:

def read_cp(fname):
    wr.init_like(fname)
    return read_vel(fname)


# In[15]:

def fullread(fname):
    t, dt, Dp, Re, vel = read_cp(fname)

    velt = new_z_vfield()
    om = new_z_vfield()
    p = new_z_pfield()

    rp(time, vel, velt, om)
    pres(velt, p, 0.0)    
    
    return t, dt, Dp, Re, vel, velt, om, p


# In[16]:

def write_cp(fname, vel, t=0.0, dt=0.1, Dp=0.0, Re=None):
    if(Re != None): wr.Re.Re = Re
    wr.down(fname, vel[0].T, vel[1].T, vel[2].T, t, dt, Dp)


# In[ ]:



