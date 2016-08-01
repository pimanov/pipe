
# coding: utf-8

# In[114]:

import numpy as np
import math


# In[115]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# # Grid

# In[116]:

class Rrt:
    def __init__(self, y0, y01):
        yp = lambda a,b: 0.01 * math.sqrt(0.3*(b+1)) * (a*(58.*b-21.-21.*b*b)+100.)
        y1p = lambda a,b: 1. + a * (b - 0.45 * (b+1.)**2)
        dyda = lambda a,b: 0.01 * math.sqrt(0.3*(b+1)) * (58.*b - 21. - 21.*b*b)
        dydb = lambda a,b: 0.01 * math.sqrt(0.3*(b+1)) * (0.5 * (a*(58.*b - 21. - 21.*b*b) + 100.)/(b+1) + a*(58.-42*b))
        dy1da = lambda a,b: b - 0.45 * (b+1.)**2
        dy1db = lambda a,b: a * (1. - 0.9*(b+1.))

        a = 0.0
        b = 10. / 3.*y0*y0 - 1.

        while True:
            c1 = dyda(a,b)
            c2 = dydb(a,b)
            c3 = y0 - yp(a,b)
            c4 = dy1da(a,b)
            c5 = dy1db(a,b)
            c6 = y01 - y1p(a,b)
            d  = c1*c5 - c2*c4
            da = (c3*c5 - c2*c6) / d
            db = (c1*c6 - c3*c4) / d
            a = a + da
            b = b + db
            if abs(da) < 1.e-5 and abs(db) < 1.e-5: 
                break

        self.aset = a
        self.bset = b

    i0 = lambda self, x: x * (1. + self.aset * (1. - x**2) * (self.bset - x**2))
    i1 = lambda self, x: 1. + self.aset * (self.bset - 3.*(1.+self.bset)*x**2 + 5.*x**4)


# In[117]:

Xmax,Zmax,epsr,nsym = 4*[float("nan")]
xf,xn,hx,lx,Im = None, None, float('nan'), 0, 1
rn,yn,rt1,rf,yf,yt1,Jm = 6*[None] + [1]
thn,thf,ht,lt,Km = None, None, float('nan'), 0, 1
Re = float('nan')
cf = float('nan')


# In[118]:

def init_r():
    global rn, yn, rt1, rf, yf, yt1
    
    hr = 1.0 / Jm
    rrt = Rrt(1.0, epsr)

    ros = np.linspace(0, 1+hr, Jm+2)
    rn = np.array([rrt.i0(ro) for ro in ros])
    yn = 1.0 - rn
    rt1 = np.array([rrt.i1(ro)*hr for ro in ros])

    ros -= hr/2
    rf = np.array([rrt.i0(ro) for ro in ros])
    yf = 1.0 - rf
    yt1 = np.array([rrt.i1(ro)*hr for ro in ros])
    
    

def init_th():
    global thn, thf, ht, Km, Zmax
    
    Km = 2**lt
    Zmax = 2.4 / nsym
    ht = Zmax / Km
    thn = np.linspace(0, Zmax+ht, Km+2)
    thf = thn - 0.5*ht

def init_x():
    global xn, xf, hx, Im
    
    Im = 2**lx
    hx = Xmax / Im
    xn = np.linspace(0, Xmax+hx, Im+2)
    xf = xn - 0.5*hx

def __str__():
    ret = "## duct sym geometry module\n"
#    ret += "## characteristic units is R and Umax=2*Ubulk \n"
    ret += "Xmax * Ymax * Zmax = %f x %f x %f\n" % (Xmax, 1.0, Zmax)
    ret += "Im * Jm * Km = %d x %d x %d\n" % (Im,Jm,Km)
    ret += "epsr=%f, nsym=%f, hx=%f, ht=%f\n" % (epsr, nsym, hx, ht)
    ret += "Re=%f, cf=%f" % (Re, cf)
    return ret


# In[119]:

print __str__()


# # Calc

# In[120]:

def mod_bc(vel):
        
    vel[:,:,:,0] = vel[:,:,:,-2]
    vel[:,:,:,-1] = vel[:,:,:,1]
    
    u,v,w = vel
    
    u[0,:,:] = u[1,:,:]
    u[-1,:,:] = u[-2,:,:]
    
    u[:,0,:] = u[:,1,:]
    u[:,-1,:] = - u[:,-2,:]

    v[:,0,:] = 0.0
    v[:,-2,:] = 0.0
    v[:,-1,:] = float('nan')
            
    v[0,:,:] = v[1,:,:]
    v[-1,:,:] = v[-2,:,:]

    w[0,:,:] = 0.0
    w[-2,:,:] = 0.0
    w[-1,:,:] = float('nan')

    w[:,0,:] = w[:,1,:]
    w[:,-1,:] = - w[:,-2,:]


# In[121]:

def get_om(vel):
    mod_bc(vel)
    u,v,w = vel
    
    om = np.zeros((3,Km+2,Jm+2,Im+2))
    ox,on,ot = om
    
    for i in range(1,Im+1):
        for j in range(0,Jm+1):
            for k in range(0,Km+1):
                w0 = w.T[i,j,k]
                w1 = w.T[i,j+1,k]
                v0 = v.T[i,j,k]
                v1 = v.T[i,j,k+1]
                ox.T[i,j,k] = (w1 - w0) / rt1[j] - (v1 - v0) / ht
                
    for k in range(0,Km+1):
        for j in range(1,Jm+1):
            for i in range(0,Im+1):
                u0 = u.T[i,j,k]
                u1 = u.T[i,j,k+1]
                w0 = w.T[i,j,k]
                w1 = w.T[i+1,j,k]
                on.T[i,j,k] = (u1 - u0) / ht - (w1 - w0) / hx
                
    for k in range(1,Km+1):
        for j in range(0,Jm+1):
            for i in range(0,Im+1):
                u0 = u.T[i,j,k]
                u1 = u.T[i,j+1,k]
                v0 = v.T[i,j,k]
                v1 = v.T[i+1,j,k]
                ot.T[i,j,k] = (v1 - v0) / hx - (u1 - u0) / rt1[j]
    
    return om


# In[122]:

def get_ox(v,w):
    ox = np.zeros_like(v)

    for i in range(1,Im+1):
        for j in range(0,Jm+1):
            for k in range(0,Km+1):
                w0 = w[k,j]
                w1 = w[k,j+1]
                v0 = v[k,j]
                v1 = v[k+1,j]
                ox[k,j] = (w1 - w0) / rt1[j] - (v1 - v0) / ht
                
    return ox


# # Read/Write

# In[123]:

import pipe_pytools.tools as tools


# In[124]:

def init():
    init_x()
    init_r()
    init_th()


# In[125]:

def read_dcp(fname, is_cf_in=False):
    
    global cf,Re,Xmax,epsr,lx,Jm,lt,nsym
    t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym,vel = tools.get_dcp(fname)
    
    if is_cf_in: 
        cf = Dp
        
    init()
    mod_bc(vel)
    
    return t,dt,vel


# In[126]:

def write_dcp(fname, t, dt, vel, is_cf_in=True):
    Dp = cf if is_cf_in else 0.0
    tools.put_dcp(fname, t, dt, Dp, Re, Xmax, epsr, lx, Jm, lt, nsym, vel)


# # Calc

# # Elementary

# In[127]:

def thmean(vel):
    VEL = np.empty_like(vel)
    VEL[:] = vel[1:-1].mean(0)
    return VEL


# In[128]:

def xmean(vel):
    VEL = np.empty_like(vel)
    VEL.T[:] = vel.T[1:-1].mean(0)
    return VEL


# In[129]:

def cs_mean(u):
    return ((u[1:-1].mean(0).T * rt1)[1:-1].sum(-1))


# # Plots

# In[130]:

def cs_plot(*args, **kwargs):
    plt.contourf(*args, **kwargs)
    plt.xlim(0,Zmax)
    plt.ylim(0,1)
    plt.gca().set_aspect('equal')


# In[131]:

def aplot(n):
    a = np.loadtxt("a0.dat", usecols=[0,n], unpack=True)
    plt.plot(a[0],a[n])


# In[132]:

def xmean_decomp(vel):
    VEL = xmean(vel)
    puls = (vel**2 - VEL**2).sum(0)[:,:,1:-1].mean(-1)**0.5
    return VEL[:,:,:,1], puls


# In[133]:

def xmeans_plot(vel):
    (U,V,W),puls = xmean_decomp(vel)
    OX = get_ox(V,W)
    
    plt.subplot(1,3,1)
    umax = U[1:-1,1:-1].max()
    cs_plot(thf, yf, U.T, np.linspace(0, 1, 21))
    plt.xlabel("U, max=%f" % umax, fontsize=14)
    plt.ylabel("x averaged", fontsize=14)
    
    plt.subplot(1,3,2)
    oxm = OX[1:-1,1:-1].max()
    cs_plot(thn, yn, OX.T, np.linspace(-oxm, oxm, 21))
    plt.xlabel("OX, max=%f" % oxm, fontsize=14)
    
    plt.subplot(1,3,3)
    pm = puls[1:-1,1:-1].max()
    cs_plot(thn, yn, puls.T, np.linspace(-pm, pm, 21))
    plt.xlabel("puls, max=%f" % pm, fontsize=14)


# In[ ]:



