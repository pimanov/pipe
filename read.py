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


def sas_read( fn ):
    
    with open( fn, "rb" ) as f:
    
        _, Re, Xmax, epsr, lx, Jm, kas, _ = i(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f)
    
        Im = 2**lx
        
        res = []
        while True:
            s = f.read(4)
            if not s: break
            t, _ = d(f), i(f)
            u = load.fget_2Dfield( f, Im, Jm )
            res.append((t,u))
            
    
        x = com.xgrid( Im, Xmax )
        r = com.rgrid( Jm, epsr )
    
        return (x, r), res


def sap_read( fn ):

    with open( fn, "rb" ) as f:

        _, t, dt, Re, epsr, Jm, lt, nsym, alpha, _ = i(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), d(f), i(f)

        Km = 2**lt

        (u1, v1, w1), (u2, v2, w2) = load.fget_2Dcvfield( f, Jm, Km )

        r = com.rgrid( Jm, epsr )
        th = com.symthgrid( Km, nsym )

        return (r, th), (u1, v1, w1), (u2, v2, w2)


def bsp_read( fn ):
    
    with open( fn, "rb" ) as f:
 
        _, cf, Xmax, epsr, lx, Jm, lt, nsym, _ = i(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f), i(f)

        Im = 2**lx
        Km = 2**lt

        (u, v, w) = load.fget_vfield( f, Im, Jm, Km )

        x = com.xgrid( Im, Xmax )
        r = com.rgrid( Jm, epsr )
        th = com.symthgrid( Km, nsym )

        return cf, (x, r, th), (u, v, w)


def msp_read( fn ):

    with open( fn, "rb" ) as f:

        _, t, dt, Dp, Re, Xmax, epsr, lx, Jm, lt, nsym, _ = i(f), d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f), i(f) 

        Im = 2**lx
        Km = 2**lt

        (u, v, w) = load.fget_vfield( f, Im, Jm, Km )
        (uu, vv, ww) = load.fget_vfield( f, Im, Jm, Km )
        (uv, vw, wu) = load.fget_vfield( f, Im, Jm, Km )

        x = com.xgrid( Im, Xmax )
        r = com.rgrid( Jm, epsr )
        th = com.symthgrid( Km, nsym )

        return Re, (x,r,th), (u,v,w), (uu,vv,ww), (uv,vw,wu)
