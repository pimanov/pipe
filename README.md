# pipe
Fortran routine to simulate Navier-Stokes equations in pipe geometry

How is this repository?
=======================
primary branches: 
+ pipe - simplest version of code to time intgration and support tools like init program. 
+ pipeSym - code to time integration with symmetry reduction (reflection and rotation symmetry).
+ pipeMPI - MPI version of code to time integration, without init or other support tools. 
+ pipeSymMPI - MPI version of code to time integration with simmetry reduction (reflection and rotation symmetry), without other tools. 
+ pytools - python scriptes to manipulation with control points files (change of greed, Re, sym, ... ). 
+ master - nothing special

Other branches are produced by primary branches with name like (primary_branch_name)-(feature) and contain some spetial cases of code.

  Expirements used code is supplied with tag


# 

branch pipeSymMPI-3Dlin
Расчет эволюции малых возмущений, наложенных на базовое течение. 

lin - не меняется, так как неявно решниется только вязкое слагаемое, а оно такое же, как и в общем случае. 
step - меняется только в том смысле, что в него передаются новые параметры для rp, и расход теперь 0. 
cf - определяется однозначно базовым течением - это скорость той системы отсчета, в которой оно покоется. 

