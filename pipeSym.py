
# coding: utf-8

# In[2]:

import pipeSym_wrapper as wr


# In[35]:

def update_ht(ht):
    wr.dimt.ht = ht
    Km = wr.dimt.Km
        
    for k in range(0,Km):
        wr.rlt.rlt[k+1] = (2.0 / ht * _math.sin(0.5 * k * _math.pi / Km))**2


# In[36]:

def update_hx(hx):
    wr.dimx.hx = hx
    Im = wr.dimx.Im
    wr.dim.Xmax = wr.dimx.hx * Im

    for i in xrange(0, Im/2 +1):
        wr.rlx.rlx[i+1] = (2.0 / hx * _math.sin(i * _math.pi / Im))**2

    for i in range(1, Im/2):
        i1 = Im + 1 - i
        wr.rlx.rlx[i1] = wr.rlx.rlx[i+1]
    


# In[38]:

def divmax(vel):
    return wr.divmax(vel[0], vel[1], vel[2])


# In[41]:

def ffmean(u):
    return wr.mean(u)


# In[18]:

def init(Xmax=None, epsr=None, nsym=None, lx=None, lt=None, Im=None, Jm=None, Km=None, 
         vel=None, hx=None, ht=None, Re=None):
    
    if(Re != None): wr.Re.Re = Re
    if(Xmax != None): wr.dim.Xmax = Xmax
    if(epsr != None): wr.dim.epsr = epsr
        
    if(nsym != None):
        if(nsym % 1.0 == 0.0):
            wr.dim.nsym = nsym
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


# In[13]:

def step(t, dt, vel, velt, om, p, cf, tol=0.01):
    t, dt = wr.step(t, dt, tol, vel[0], vel[1], vel[2], velt[0], velt[1], velt[2], om[0], om[1], om[2], p, cf)
    return t, dt, vel, velt, om, p


# In[15]:

def rp(t, vel, velt, om):
    wr.rp(t, vel[0], vel[1], vel[2], velt[0], velt[1], velt[2], om[0], om[1], om[2])
    return velt, om


# In[17]:

def pres(vel, p, ub):
    wr.pres(vel[0], vel[1], vel[2], p, ub)
    return vel, p


# In[11]:

import numpy as _np
import math as _math

def calc(vel, dtmax, cf=0.0, Re=None, velt=None, om=None, p=None, maxnstep=None, tmax=None, ht=None, hx=None, 
         tol=0.01, time=0.0, dt=None, prt=None, ampmax=None, ampmin=None, const_dt_chec=False):    
        
    if(hx != None): update_hx(hx)
    if(ht != None): update_ht(ht)
    if(Re != None): wr.Re.Re = Re
    if(dt == None): dt = dtmax

    if(om == None):
        om = _np.zeros_like(vel)
        
    if(velt == None): 
        velt = _np.zeros_like(vel)
        wr.rp(time, vel[0], vel[1], vel[2], velt[0], velt[1], velt[2], om[0], om[1], om[2])
        
    if(p == None):
        p = _np.zeros_like(vel[0])
        wr.pres(vel[0], vel[1], vel[2], p, 0.5)
        
    nstep = 0
    break_cond_chec = False
    while(true):
        if(tmax != None):
            break_cond_chec = True
            if(time > tmax + 0.1*dt): 
                print "Planned break! tmax reached"
                break
    
        time, dt1 = step(time, dt, tol, vel[0], vel[1], vel[2], velt[0], velt[1], velt[2], om[0], om[1], om[2], p, cf)
        
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
        if(nstep != None):
            break_cond_chec = True
            if (nstep >= maxnstep): 
                print "Planned break! maxnstep reached"
                break
                
                
        if(break_cond_chec == False):
            print "No break condition set. Forced break!"
            break
    
    return time, dt, vel, vel1, om, p


# In[46]:

def new_z_vfield(Im = None, Jm = None, Km = None):
    if(Im == None): Im = wr.dimx.Im
    if(Jm == None): Jm = wr.dimr.Jm
    if(Km == None): Km = wr.dimt.Km
    return _np.zeros((3,Km+2, Jm+2, Im+2), order='C')


# In[47]:

def new_z_pfield(Im = None, Jm = None, Km = None):
    if(Im == None): Im = wr.dimx.Im
    if(Jm == None): Jm = wr.dimr.Jm
    if(Km == None): Km = wr.dimt.Km
    return _np.zeros((Km+2, Jm+2, Im+2), order='C')


# In[59]:

def read_vel(fname):
    vel = new_z_vfield()
    t, dt, Dp = wr.load(fname, vel[0], vel[1], vel[2])
    Re = wr.Re.Re
    return t, dt, Dp, Re, vel


# In[60]:

def read_cp(fname):
    wr.init_like(fname)
    return read_vel(fname)


# In[61]:

def fullread(fname):
    t, dt, Dp, Re, vel = read_cp(fname)

    velt = new_z_vfield()
    om = new_z_vfield()
    p = new_z_pfield()

    wr.rp(time, vel[0], vel[1], vel[2], velt[0], velt[1], velt[2], om[0], om[1], om[2])
    wr.pres(vel[0], vel[1], vel[2], p, 0.5)    
    
    return t, dt, Dp, Re, vel, velt, om, p


# In[64]:

def write_cp(fname, vel, t=0.0, dt=0.1, Dp=0.0, Re=None):
    if(Re != None): wr.Re.Re = Re
    wr.down(fname, vel[0], vel[1], vel[2], t, dt, Dp)


# In[ ]:



