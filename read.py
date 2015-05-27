import com
import load
import numpy as np

d = load.fget_d
i = load.fget_i



def cp_read(cp_name):
    
    with open(cp_name, 'rb') as f:
        
        _, t, dt, Dp, Re, Xmax, epsr, lx, Jm, lt, _ = i(f), d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f)

        Im = 2**lx
        Km = 2**lt
        u, v, w = load.fget_vfield(f, Im, Jm, Km)
    
    x = com.xgrid( Im, Xmax )
    r = com.rgrid( Jm, epsr )
    th = com.thgrid( Km )
    
    return Re, t, (x, r, th), (u, v, w)



def scp_read(file_name):
    
    with open( file_name, "rb") as f:
    
        _, t, dt, Dp, Re, Xmax, epsr, lx, Jm, lt, nsym, _ = i(f), d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f), i(f)
    
        Im = 2**lx
        Km = 2**lt
        u, v, w = load.fget_vfield( f, Im, Jm, Km )
    
    x = com.xgrid( Im, Xmax )
    r = com.rgrid( Jm, epsr )
    th = com.symthgrid( Km, nsym )
    
    return Re, t, (x, r, th), (u, v, w)



def bss_read( fn ):
    
    with open( fn, "rb" ) as f:
    
        _, epsr, Jm ,lt, nsym, _ = i(f), d(f), i(f), i(f), i(f), i(f)
    
        Km = 2**lt
        u,v,w = load.fget_2Dvfield( f, Jm, Km )
    
    r = com.rgrid( Jm, epsr )
    th = com.symthgrid( Km, nsym )
    
    return (r, th), (u, v, w)


scs_read = bss_read


    