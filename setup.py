#!/usr/bin/env python2 

from numpy.distutils.core import Extension

ext1 = Extension(name = 'ductSym_wrap',
                 sources = ['ductSym_wrap.pyf', 'com.for', 'divmax.for', 'fft8.for', 'bc_om.for',  
                            'lin.for', 'read.for', 'pres.for', 'prog3.for', 'rp.for', 'visc.for',
                            'rrt8.for', 'write.for', 'step.for', 'init_like.for', 'add_nl.for']
                )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'ductSym_wrap',
          description       = "ductSym core fortran routines wrapper",
          author            = "Vladimir Pimanov",
          author_email      = "pimanov-vladimir@yandex.ru",
          ext_modules = [ext1]
         )
     

# End of setup.py
