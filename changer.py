#!/usr/bin/python
import math
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
  
def fput_x(f):
  res = struct.pack('i',-1)
  f.write(res)

def fput_d(f,p):
  res = struct.pack('d',p)
  f.write(res)

def fput_i(f,p):
  res = struct.pack('i',p)
  f.write(res)

def fget_param(f):
  d = fget_d
  i = fget_i
  fget_x(f)
  res = ( d(f), d(f), d(f), d(f), d(f), d(f), i(f), i(f), i(f), i(f) )
  fget_x(f)
  return res
  
def fput_param(f, d1, d2, d3, d4, d5, d6, i7, i8, i9, i0):
  fput_x(f)
  fput_d(f,d1)
  fput_d(f,d2)
  fput_d(f,d3)
  fput_d(f,d4)
  fput_d(f,d5)
  fput_d(f,d6)
  fput_i(f,i7)
  fput_i(f,i8)
  fput_i(f,i9)
  fput_i(f,i0)
  fput_x(f)
  return

def fget_dline(f,n):
  fget_x(f)
  dline = []
  for i in range(n):
    dline.append(fget_d(f))
  fget_x(f)
  return dline

def fput_dline(f,l):
  fput_x(f)
  for p in l:
    fput_d(f,p)
  fput_x(f)

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

def fput_vfield(f,vec):
  for q in zip(vec[0],vec[1],vec[2]):
    for l in zip(q[0],q[1],q[2]):
      fput_dline(f,l[0])
      fput_dline(f,l[1])
      fput_dline(f,l[2])

def alpha_reducting(alpha,vec):
  u = vec[0]
  v = vec[1]
  w = vec[2]

  Km = len(u)
  Jm = len(u[0])
  Im = len(u[0][0])

  for i in range(Im):
    for j in range(Jm):
      ua = 0.0
      va = 0.0
      wa = 0.0

      for k in range(Km):
        ua += u[k][j][i]
        va += v[k][j][i]
        wa += w[k][j][i]

      ua /= Km
      va /= Km
      wa /= Km

      u[k][j][i] = alpha * u[k][j][i] + (1.0 - alpha) * ua
      v[k][j][i] = alpha * v[k][j][i] + (1.0 - alpha) * va
      w[k][j][i] = alpha * w[k][j][i] + (1.0 - alpha) * wa

  return

import optparse

def parse_args():
  parser = optparse.OptionParser(usage="Usage: %prog IN <options>", version='%prog version 0.1', description="change some value of cp stored in IN file")
  
  parser.add_option('-o', '--out', help='set output to file OUT', default='ch.cp')
  parser.add_option('-r', '--Re', help='chage the Reynolds number to RE', type='float')
  parser.add_option('-t', '--time', help='change the time to TIME', type='float')
  parser.add_option('-d', '--dt', help='change the dt to DT', type='float')
  parser.add_option('-l', '--Xmax', help='change the Xmax to XMAX', type='float')
  parser.add_option('-e', '--epsr', help='change the epsr to EPSR', type='float')
  parser.add_option('-n', '--nsym', help='change the nsym to NSYM', type='int')
  parser.add_option('-a', '--alpha', help='set and rate of 3D structure', metavar='A', type='float')
  
  (opts,args) = parser.parse_args()
  
  res = {}
  if not opts.out is None: res['out'] = opts.out
  if not opts.Re is None: res['Re'] = opts.Re
  if not opts.time is None: res['time'] = opts.time
  if not opts.dt is None: res['dt'] = opts.dt
  if not opts.Xmax is None: res['Xmax'] = opts.Xmax
  if not opts.epsr is None: res['epsr'] = opts.epsr
  if not opts.nsym is None: res['nsym'] = opts.nsym
  if not opts.alpha is None: res['Xmax'] = opts.alpha
  if len(args) is 1: res['in'] = args[-1]
  else: 
    parser.print_help()
    exit(-1)

  return res


def main():
  opts = parse_args()
 
  with open(opts['in'], 'rb') as cp_file:
    t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt,nsym = fget_param(cp_file)
    Im = 2**lx
    Km = 2**lt
    vec = fget_vfield(cp_file,Im,Jm,Km)
    lines = cp_file.readlines()

  print 'input file %s was read: '%(opts['in'])
  print '    Re=%f, time=%f, dt=%f'%(Re,t,dt)
  print '    Xmax=%f, epsr=%f, nsym=%d'%(Xmax,epsr,nsym)
  print '    Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km)
  print '>>  ' + '    '.join(lines)

  logs = []
  param = {'Re':Re, 'time':t, 'dt':dt, 'Xmax':Xmax, 'epsr':epsr, 'nsym':nsym}
  for key in opts:
    if key is 'alpha':
      alpha = opts['alpha']
      logs.append('3D structure of flow have been reducted %f times'%(alpha))
      alpha_reducting(vec, alpha)
    elif not key is 'out' and not key is 'in':
      logs.append('%s have been changed from %f to %f'%(key,param[key],opts[key]))
      param[key] = opts[key]

  if len(logs) == 0: logs = ['No changes']
  log = 'changer.py: ' + ', '.join(logs) + '\n'
  print log

  with open( opts['out'], 'wb') as cp_file:
    fput_param( cp_file, param['time'], param['dt'], Dp, param['Re'], param['Xmax'], param['epsr'], lx, Jm, lt, param['nsym'])
    fput_vfield( cp_file, vec)
    cp_file.writelines( lines)
    cp_file.writelines( [log])

  print 'new control point have been saved in %s'%(opts['out'])
  print '    Re=%f, time=%f, dt=%f'%(param['Re'],param['time'],param['dt'])
  print '    Xmax=%f, epsr=%f, nsym=%d'%(param['Xmax'],param['epsr'],param['nsym'])
  print '    Im=%d, Jm=%d, Km=%d'%(Im,Jm,Km)
  print ''
  
if __name__ == '__main__': main()

