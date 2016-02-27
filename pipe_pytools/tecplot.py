
import math
import numpy as np

def tec_write(plot_name, grid, u, step = (1,1,1), lim = (1,-1), ishift = 0):

    plot_file = open(plot_name,'w')

    x,r,o = grid
    xstep, rstep, ostep = step
    xmin, xmax = lim
    if (xmax < 0): xmax += x.m

    ucl = u[1:-2, 1, :].mean(0)
    

    plot_file.write('ZONE I=%d, J=%d, K=%d F=POINT \n'%((xmax-xmin)/xstep +1, r.m/rstep +2, 4*o.m/ostep +1))
    
    for kk in range(0,o.m*4+1,ostep):
        if (kk == 0):
            k = 1
        elif (kk <= o.m):
            k = kk
        elif (kk <= 2*o.m):
            k = 2*o.m - kk + 1
        elif (kk <= 3*o.m):
            k = kk - 2*o.m
        else:
            k = 4*o.m - kk + 1
        
        for j in [0]: 
            for i in range(xmin,xmax+1,xstep):
                ii = (i+ishift)%x.m
                plot_file.write( '%e %e %e %e %e \n' % ( x.face[i], 0.0, 0.0, ucl[ii] - 1.0, 0.0 ) )
        
        for j in range(rstep,r.m+1,rstep):
            for i in range(xmin,xmax+1,xstep):
                ii = (i+ishift)%x.m
                y = r.face[j] * math.cos( o.h*(kk-0.5) )
                z = r.face[j] * math.sin( o.h*(kk-0.5) )
                plot_file.write( '%e %e %e %e %e \n' % 
                    ( x.face[i], y, z, u[k,j,ii] - (1.0-r.face[j]**2), r.face[j]))
               
        for j in [r.m+1]:
            for i in range(xmin,xmax+1,xstep):
                ii = (i+ishift)%x.m
                y = 1.0 * math.cos( o.h*(kk-0.5) )
                z = 1.0 * math.sin( o.h*(kk-0.5) )
                plot_file.write( '%e %e %e %e %e \n' % ( x.face[i], y, z, 0.0, 1.0 ) )
            
    plot_file.close()
    
