
# coding: utf-8

# In[2]:

import ductSym_wrap as wr
import numpy as np
import math


def look():
    print "Xmax=%f, epsr=%f, dsym=%f" % (wr.dim.Xmax, wr.dim.epsr, wr.dim.dsym)
    print "Im=%d, Jm=%d, Km=%d" % (wr.dimx.Im, wr.dimr.Jm, wr.dimt.Km)
    print "hx=%f, ht=%f" % (wr.dimx.hx, wr.dimt.ht)
    print "Re=%f" % (wr.Re.Re), "duct"

# In[13]:

def new_z_vfield(Im=None, Jm=None, Km=None):
    if Im == None: Im = wr.dimx.Im
    if Jm == None: Im = wr.dimr.Jm
    if Km == None: Im = wr.dimt.Km
    return np.zeros((3, wr.dimt.Km+2, wr.dimr.Jm+2, wr.dimx.Im+2), order='C')


# In[12]:

def new_z_pfield(Im=None, Jm=None, Km=None):
    if Im == None: Im = wr.dimx.Im
    if Jm == None: Im = wr.dimr.Jm
    if Km == None: Im = wr.dimt.Km
    return np.zeros((wr.dimt.Km+2, wr.dimr.Jm+2, wr.dimx.Im+2), order='C')


# In[7]:

def divmax(vel):
    return wr.divmax(vel[0].T, vel[1].T, vel[2].T)


# In[8]:

def init(Xmax=None, epsr=None, dsym=None, lx=None, Jm=None, lt=None, Re=None, Im=None, Km=None, vel=None):
    
    if(Re != None): wr.Re.Re = Re
    if(Xmax != None): wr.dim.Xmax = Xmax
    if(epsr != None): wr.dim.epsr = epsr
    if(dsym != None): wr.dim.dsym = dsym
        
    if(vel != None): 
        Im = vel.shape[1] - 2
        Jm = vel.shape[2] - 2
        Km = vel.shape[3] - 2
        
    if(Im != None): 
        lx = int(round(math.log(Im,2)))
        if(2**lx != Im):
            print "###Error! Im != 2**n"
            return
            
    if(Km != None): 
        lt = int(round(math.log(Km,2)))
        if(2**lt != Km):
            print "###Error! Km != 2**n"
            return
        
    if(lx != None): wr.dimx.lx = lx
    if(Jm != None): wr.dimr.Jm = Jm
    if(lt != None): wr.dimt.lt = lt
        
    wr.com()
    return


# In[9]:

def step(t, dt, vel, velt, om, p, cf, tol=0.01):
    t, dt = wr.step(t, dt, tol, vel[0].T, vel[1].T, vel[2].T, velt[0].T, velt[1].T, velt[2].T, om[0].T, om[1].T, om[2].T, p.T, cf)
    return t, dt, vel, velt, om, p


# In[10]:
def bc_om(vel, om=None):
    if type(om) == type(None): om = new_z_vfield()
    wr.bc_om(vel[0].T, vel[1].T, vel[2].T, om[0].T, om[1].T, om[2].T)
    return om

def add_nl(vel, om, velt=None):
    if type(velt) == type(None): velt = new_z_vfield()
    wr.add_nl(vel[0].T,vel[1].T,vel[2].T,om[0].T,om[1].T,om[2].T,velt[0].T,velt[1].T,velt[2].T)
    return velt 

def visc(om, velt=None):
    if type(velt) == type(None): velt = new_z_vfield()
    wr.visc(om[0].T,om[1].T,om[2].T,velt[0].T,velt[1].T,velt[2].T)
    return velt

def rp(vel, t=0.0, velt=None, om=None):
    if type(velt) == type(None): velt = new_z_vfield()
    if type(om) == type(None): om = new_z_vfield()
    wr.rp(t, vel[0].T, vel[1].T, vel[2].T, velt[0].T, velt[1].T, velt[2].T, om[0].T, om[1].T, om[2].T)
    return velt, om


# In[11]:

def pres(vel, ub, inplace=False, p=None):
    if not inplace:
        vel = vel.copy()

    if type(p) == type(None): p = new_z_pfield()

    wr.pres(vel[0].T, vel[1].T, vel[2].T, p.T, ub)
    return vel, p

def grad(vel, p, inplace=False):
    if not inplace:
        vel = vel.copy()

    wr.gradp(vel[0].T, vel[1].T, vel[2].T, p.T)
    return vel

# In[11]:

def calc(vel, dtmax, cf=0.5, Re=None, velt=None, om=None, p=None, maxnstep=None, tmax=None, 
         tol=0.01, time=0.0, dt=None, prt=None, ampmax=None, ampmin=None, const_dt_chec=False, inplace=False, ready=False):    
    """
    amp = prt(nstep, time, dt, vel, velt, om, p, cf)
    return time, dt, vel, velt, om, p
    """
        
    if not inplace:
        vel = vel.copy()
        
    if(Re != None): wr.Re.Re = Re
    if(dt == None): dt = dtmax

    if not ready: 
        velt, om = rp(vel, time, velt, om)
        velt, p = pres(velt, 0.0, True, p)
        
    nstep = 0
    break_cond_chec = False
    while True:
        if(tmax != None):
            break_cond_chec = True
            if(time > tmax + 0.1*dt): 
                print "Planned break! tmax reached"
                break
    
        time1, dt1, vel, velt, om, p = step(time, dt, vel, velt, om, p, cf, tol)
        
        if(const_dt_chec == True):
            if(dt != dtmax):
                print "###Forced break! dt != dtmax"
                break

            if(time1 != time + dt): 
                print "###Forced break! dt has changed in step()"
                break

            time = time1
        else:
            dt = min(dt1, dtmax)
            time = time1

        if(prt != None): 
            amp = prt(nstep, time, dt, vel, velt, om, p, cf)
            
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

def read_dcp(fname):
    wr.init_like_dcp(fname)
    vel = new_z_vfield()
    t,dt,Dp = wr.read_dcp(fname, vel[0].T, vel[1].T, vel[2].T)
    return t,dt,Dp,vel

# In[15]:

def fullread_dcp(fname):
    t, dt, Dp, vel = read_dcp(fname)

    velt, om = rp(vel, t)
    velt, p = pres(velt, 0.0, True)    
    
    return t, dt, Dp, vel, velt, om, p

# In[16]:


def write_dcp(fname, vel, t=0.0, dt=0.25, Dp=0.0, Re=None):
    if(Re != None): wr.Re.Re = Re
    wr.write_dcp(fname, vel[0].T, vel[1].T, vel[2].T, t, dt, Dp)

# In[ ]:



