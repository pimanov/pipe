#!/usr/bin/env python
 
from distutils.core import setup
 
setup(name='pipe_pytools',
 version='0.1',
 description='python tools vor pipe programs',
 author='Vladimir Pimanov',
 author_email='pimanov-vladimir@yandex.ru',
 url='https://github.com/pimanov/pipe',
 install_requires=['numpy', 'scipy', 'math', 'os'],
 py_modules = ['com','forio','cprw','elementary']
)
