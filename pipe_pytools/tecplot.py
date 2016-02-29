import math as _m


def plot(fname, grid, u, uname='U', step=(1,1,1), ishift=0, time=0.0, nframe=0):
    X,R,Th = grid
    istep,jstep,kstep = step

    if (n == 0): mode = 'w'
    else: mode = 'a'

    with open(fname, mode) as f:
        if (nframe == 0): f.write('VARIABLES="X" "Y" "Z" "R" "%s"\n' % uname)
        f.write('ZONE I=%d J=%d K=%d\n' % (X.m/istep, R.m/jstep+2, 2*Th.nsym*Th.m/kstep+1))
        f.write('SOLUTIONTIME=%f\n' % time)
        if (nframe > 0): f.write('VARSHARELIST=([1-4])\n')

        ucl = u[1:-1,1,:].mean(0)
        
        for kk in range(1, 2*Th.nsym*Th.m+2, kstep):
            if (kk != 2*Th.nsym*Th.m+1):
                k = kk
                while k > Th.m:
                    k -= 2*Th.m
                if (k <= 0): 
                    k = 1 - k
            else:
                k = 1
                
            for j in range(0, Jm+2, jstep):
                for i in range(1, Im+1, istep):
                    x = X.n[i]
                    th = kk*Th.h
                    ii = (i + ishift - 1) % Im + 1
                    
                    if (j==0): 
                        r=0.0
                        uu = ucl[ii] - 1.0
                    elif (j==Jm+1): 
                        r=1.0
                        uu = 0.0
                    else: 
                        r = R.f[j]
                        uu = u[k,j,ii] - (1.0 - r**2)
                    
                    y = r * math.cos(th)
                    z = r * math.sin(th)
                    
                    if (n == 0): f.write("%6.3f %6.3f %6.3f %6.3f %6.3f\n" % (x,y,z,r,u))
                    else: f.write("%6.3f\n" % u)
                    
        f.close()
    return

"""
def plot2(fname, grid, uS, unameS = [], step=(1,1,1), ishift=0, time=0.0, nframe=0):
    X,R,Th = grid
    istep,jstep,kstep = step

    for i in range(len(unameS), len(uS)):
        unameS.append('U%d' % i)

    if (n == 0): mode = 'w'
    else: mode = 'a'

    with open(fname, mode) as f:
        if (nframe == 0):
            var_header = 'VARIABLES="X" "Y" "Z" "R" '
            for uname in unameS: var_header += '"%s" ' % uname
            var_header += '\n'
            f.write(var_header)

        f.write('ZONE I=%d J=%d K=%d\n' % (X.m, R.m+2, 2*Th.nsym*Th.m+1))

        f.write('SOLUTIONTIME=%f\n' % time)

        if (nframe > 0): f.write('VARSHARELIST=([1-4])\n')

        uclS = []
        for u in uS
            uclS.append(u[1:-1,1,:].mean(0))
        
        for kk in range(1, 2*Th.nsym*Th.m+2, kstep):
            if (kk != 2*Th.nsym*Th.m+1): 
                k = kk
                while k > Th.m:
                    k -= 2*Th.m
                if k <= 0: k = 1 - k 
            else: 
                k = 1
                
            for j in range(0,Jm+2, jstep):
                for i in range(1,Im+1, istep):
                    x = X.n[i]
                    th = kk*Th.h
                    ii = (i + ishift - 1) % Im + 1
                    
                    if (j==0): 
                        r=0.0
                        uu = []
                        for u,ucl in zip(Us, uclS):
                            uu.append(ucl[ii] - 1.0)
                    elif (j==Jm+1): 
                        r=1.0
                        u = 0.0
                    else: 
                        r = R.f[j]
                        u = vel[k,j,ii] - (1.0 - r**2)
                    
                    y = r * math.cos(th)
                    z = r * math.sin(th)
                    
                    if (n == 0): f.write("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n" % (x,y,z,r,u,c))
                    else: f.write("%6.3f\n" % u)
                    
        f.close()
    return
"""
