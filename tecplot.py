
import math
import numpy as np

def plot_write(plot_name,u):
    plot_file = open(plot_name,'w')

    zshift = 7
    ishift = 300
    xmin = 200
    xmax = 800
    xstep = 2
    rstep = 2
    ostep = 2
    
    plot_file.write('ZONE I=%d, J=%d, K=%d F=POINT \n'%((xmax-xmin)/xstep +1, r.m/rstep +1, 4*o.m/ostep +1))
    
    for kk in range(0,o.m*4+1,ostep):
        if (kk < o.m):
            k = kk
        elif (kk < 2*o.m):
            k = 2*o.m - kk
        elif (kk < 3*o.m):
            k = kk - 2*o.m
        else:
            k = 4*o.m - kk
        
        # j = 0
        for i in range(xmin,xmax+1,xstep):
            ii = (i+ishift)%x.m
            plot_file.write( '%e %e %e %e %e \n' % ( x.n[i], 0.0, 0.0 + zshift, u[ii,0,k] - 1.0, 0.0 ) )
        
        for j in range(rstep,r.m,rstep):
            for i in range(xmin,xmax+1,xstep):
                ii = (i+ishift)%x.m
                y = r.n[j] * math.cos( o.h*(kk-0.5) )
                z = r.n[j] * math.sin( o.h*(kk-0.5) ) + zshift
                plot_file.write( '%e %e %e %e %e \n' % 
                    ( x.n[i], y, z, (u[ii,j,k] + u[ii,j+1,k])*0.5 - (1.0-r.n[j]**2), r.n[j]))
               
        # j = Jm
        for i in range(xmin,xmax+1,xstep):
            ii = (i+ishift)%x.m
            y = 1.0 * math.cos( o.h*(kk-0.5) )
            z = 1.0 * math.sin( o.h*(kk-0.5) ) + zshift
            plot_file.write( '%e %e %e %e %e \n' % ( x.n[i], y, z, 0.0, 1.0 ) )
            
    plot_file.close()
    
    print ("u speed have been read in %s \n" % (plot_name))
