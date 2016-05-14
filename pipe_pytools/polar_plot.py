from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import com


def polar_contourf(r, th, u, amp, cmap_param="jet"):
    (umin, umax) =amp
    
    radius = r
    angles = np.repeat( th[ ..., np.newaxis ], r.shape[0], axis = 1 )

    x = (radius * np.cos( angles )).flatten()
    y = (radius * np.sin( angles )).flatten()
    triang = tri.Triangulation(x, y)

    z = u.flatten()
    
    levels = MaxNLocator(nbins=16).tick_values(umin, umax)
    scmap = plt.get_cmap(cmap_param)

    plt.gca().set_aspect('equal')
    plt.tricontourf( triang, z, levels = levels, cmap = scmap)
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)

    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    return
    
    
def polar_contour(r, th, u, amp, ampmax=None, nlines=1, line_param='-k'):

    if type(ampmax) == type(None): 
        ampmax = amp
    
    radius = r
    angles = np.repeat( th[ ..., np.newaxis ], r.shape[0], axis = 1 )

    x = (radius * np.cos( angles )).flatten()
    y = (radius * np.sin( angles )).flatten()
    triang = tri.Triangulation(x, y)

    z = u.flatten()
    
    level0 = MaxNLocator(nbins=nlines).tick_values(amp, ampmax)

    plt.gca().set_aspect('equal')
    plt.tricontour( triang, z, colors=line_param, levels=level0)
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)

    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    return
    

def polar_arrow(v, w, rs, ths, alpha, color='k'):
    
    def x(r,th):
        return r * math.cos(th)
    
    def y(r,th):
        return r * math.sin(th)
    
    def ux(v,w,th):
        return v*math.cos(th) - w*math.sin(th)
    
    def uy(v,w,th):
        return v*math.sin(th) + w*math.cos(th)
    
    for th, v1, w1 in zip(ths, v, w):
        for r, v2, w2, in zip(rs, v1, w1):
            plt.arrow(x(r,th), y(r,th), alpha*ux(v2,w2,th), alpha*uy(v2,w2,th), head_width=0.01, fc=color, ec=color)
         
    plt.gca().set_aspect('equal')
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)

    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    return



def sym_polar_plot( grid, u, amp):
    
    (R, O) = grid
    (umin, umax) =amp

    radius = R.face[1:-1]
    angles = np.linspace( O.h/2.0, 2.0*math.pi + O.h/2.0, 4.0*O.m, endpoint = False )
    angles = np.repeat( angles[ ..., np.newaxis ], R.m, axis = 1 )

    x = (radius * np.cos( angles )).flatten()
    y = (radius * np.sin( angles )).flatten()
    triang = tri.Triangulation(x, y)

    z = angles.copy()
    z[:O.m,:] = u[1:-1, 1:-1]
            
    for k in range(O.m):
        z[2*O.m-k-1] = z[k]
        
    for k in range(2*O.m):
        z[4*O.m-k-1] = z[k]
        
    z = z.flatten()

    
    level0 = MaxNLocator(nbins=1).tick_values(umin, umax)
    levels = MaxNLocator(nbins=16).tick_values(umin, umax)
    scmap = plt.get_cmap()

    plt.gca().set_aspect('equal')
    plt.tricontourf( triang, z, levels = levels, cmap = scmap)
    plt.tricontour ( triang, z, colors = 'k', levels = level0)
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)

    return

    
    
    
def quart_polar_plot( grid, u, amp):
    
    (R, O) = grid
    (umin, umax) =amp

    radius = R.face
    radius[0] = 0.0
    radius[-1] = 1.0
    
    
    angles = O.face
    angles[0] = angles[1] - O.h/2
    angles[-1] = angles[-2] + O.h/2
    angles = np.repeat( angles[ ..., np.newaxis ], R.m+2, axis = 1 )

    x = (radius * np.cos(angles) ).flatten()
    y = (radius * np.sin(angles) ).flatten()
    triang = tri.Triangulation(x, y)

    ucl = u[1:-1,1].mean()
    
    z = u.copy()
    
    for z1 in z:
        z1[0] = ucl
        z1[-1] = 0.0
        
    z[0]  = z[1]
    z[-1] = z[-2]
    
    z = z.flatten()
    
    level0 = MaxNLocator(nbins=1).tick_values(umin, umax)
    levels = MaxNLocator(nbins=16).tick_values(umin, umax)
    scmap = plt.get_cmap()

    plt.gca().set_aspect('equal')
    plt.tricontourf( triang, z, levels = levels, cmap = scmap)
    plt.tricontour ( triang, z, colors = 'k', levels = level0)
    plt.xlim(0,1)
    plt.ylim(0,1)

    return    



def polar_plot( grid, u, amp):
    
    (R, O) = grid
    (umin, umax) =amp

    radius = R.face
    radius[0]  = 0.0
    radius[-1] = 1.0
    
    angles = O.face[1:-1]
    angles = np.repeat( angles[ ..., np.newaxis ], R.m+2, axis = 1 )

    x = (radius * np.cos( angles )).flatten()
    y = (radius * np.sin( angles )).flatten()
    triang = tri.Triangulation(x, y)

    ucl = u[1:-1,1].mean()
    
    z = u[1:-1].copy()
    for z1 in z:
        z1[0] = ucl
        z1[-1] = 0.0
        
    
    
    level0 = MaxNLocator(nbins=1).tick_values(umin, umax)
    levels = MaxNLocator(nbins=16).tick_values(umin, umax)
    scmap = plt.get_cmap()

    plt.gca().set_aspect('equal')
    plt.tricontourf( triang, z, levels = levels, cmap = scmap)
    plt.tricontour ( triang, z, colors = 'k', levels = level0)
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)

    return

    
