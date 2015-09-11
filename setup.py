# python setup.py install

from numpy.distutils.core import Extension

ext1 = Extension(name = 'pipeSym',
                 sources = ['pipeSym.pyf', 
                            'com.for', 'divmax.for', 'fft8.for', 
                            'init_like.for', 'lin.for', 'load.for', 
                            'pres.for', 'prog3.for', 'rp.for', 
                            'rrt8.for', 'down.for', 'step.for'
                           ]
                )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'pipeSym',
          description       = "pipeSym core fortran routines wrapper",
          author            = "Vladimir Pimanov",
          author_email      = "pimanov-vladimir@yandex.ru",
          ext_modules = [ext1]
         )
# End of setup.py

