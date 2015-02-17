# pipe
Fortran routine to simulate Navier-Stokes equations in pipe geometry

How is this repository?
=======================
primary branches: 
+ pipe - simplest version of code to time intgration and support tools like init program. 
+ pipeSym - code to time integration with symmetry reduction (reflection and rotation symmetry).
+ pipeMPI - MPI version of code to time integration, without init or other support tools. 
+ pipeSymMPI - MPI version of code to time integration with simmetry reduction (reflection and rotation symmetry), without other tools. 
+ pytools - python scriptes to manipulation with control points (change of greed, Re, sym to common, ... ). 
+ master - nothing special

  Файлы в различных ветках первого эшелона не могут быть общими, так как изменения в них затрагивают практически каждый файл. Идея о том, что com подпрограмма, например, могла бы быть объщая для всех верчий кода, чтобы избежать личшних ошибок и сократить время на можификацию кода, не реализуема. 
  Для того, чтобы иметь возможность совмещать в исследовании расчеты в случае с симметрией и без нее, можно выгружать различные ветки в репозиторий и, скомпилировав файл, выкладывать ошник в папку bin. Собрав в bin все необходимые ошники, начинать эксперимент. Это возможно, так как расчет проводится толькопо коду, который загружен на гитхаб. 

Производные ветки отделяются от первичных и содержат различные модификации кода для решения специальных задач. Они именуются, как <primary_branch_name>-<feature>.  

  Когда с помощью какого-то кода проводится расчет, на него ставится tag 
