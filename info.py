#!/usr/bin/python
import struct

def d(f):
    d = f.read(8)
    res, = struct.unpack('d',d)
    return res

def i(f):
    i = f.read(4)
    res, = struct.unpack('i',i)
    return res
    
def fget_param(f):
 
    t = i(f)
    nsym={60: 0, 64: 1}.get(t,-1)

    if(nsym > 0):
        res = ( d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f) )
    else:
        res = ( d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), nsym )
    fget_i(f)

    return res
  

import sys

def main():

    for cp_name in sys.argv[1:]:
    
        with open(cp_name, 'rb') as cp_file:
            t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym = fget_param(cp_file)
            Im = 2**lx
            Km = 2**lt

        if(nsym < 0): 
            print cp_name + ': type not be recognized'
        elif (nsym > 0):
            print cp_name + ': Re=%f, time=%f, dt=%f'%(Re,t,dt) +', Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km) + ', Xmax=%f, epsr=%f, nsym=%d'%(Xmax,epsr,nsym)
        else: # nsym == 0
            print cp_name + ': Re=%f, time=%f, dt=%f'%(Re,t,dt) +', Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km) + ', Xmax=%f, epsr=%f, NO SYM'%(Xmax,epsr)


if __name__ == '__main__': main()

