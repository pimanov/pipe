#!/usr/bin/python
import struct

def fget_x(f):
  f.read(4)

def fget_d(f):
  d = f.read(8)
  res, = struct.unpack('d',d)
  return res

def fget_i(f):
  i = f.read(4)
  res, = struct.unpack('i',i)
  return res

def fget_param(f):
  d = fget_d
  i = fget_i
  t=i(f)
  nsym={60: 0, 64: 1}.get(t,-1)
  if(nsym > 0):
    res = ( d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f) )
  else:
    res = ( d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), nsym )
    
  fget_x(f)
  return res
  

def fget_dline(f,n):
  fget_x(f)
  dline = []
  for i in range(n):
    dline.append(fget_d(f))
  fget_x(f)
  return dline

def flose_vfield(f,Im,Jm,Km):
  for k in range(Km):
    for j in range(Jm):
      for n in range(3):
        fget_x(f)
        for i in range(Im):
          fget_d(f)
        fget_x(f)

def fget_vfield(f,Im,Jm,Km):
  u = []
  v = []
  w = []
  for k in range(Km):
    uq = []
    vq = []
    wq = []
    for j in range(Jm):
      uq.append( fget_dline(f,Im))
      vq.append( fget_dline(f,Im))
      wq.append( fget_dline(f,Im))
    u.append(uq)
    v.append(vq)
    w.append(wq)
  return u,v,w

import sys

def main():

  for cp_name in sys.argv[1:]:
    
    with open(cp_name, 'rb') as cp_file:
      t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym = fget_param(cp_file)

      if(nsym < 0): 
        print 'type not be recognized'
        exit()

      Im = 2**lx
      Km = 2**lt
      vec = flose_vfield(cp_file,Im,Jm,Km)
      lines = cp_file.readlines()

    if (nsym > 0):
      print cp_name + ': Re=%f, time=%f, dt=%f'%(Re,t,dt) +', Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km) + ', Xmax=%f, epsr=%f, nsym=%d'%(Xmax,epsr,nsym) + '; '.join(lines)
    else:
      print cp_name + ': Re=%f, time=%f, dt=%f'%(Re,t,dt) +', Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km) + ', Xmax=%f, epsr=%f, NO SYM'%(Xmax,epsr) + '; '.join(lines)
  
if __name__ == '__main__': main()

