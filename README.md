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

Программа для расчета эволюции малых возмущений, развивающихся на потоке, однородном вроль трубы, постоянном по времени. Базовое течение необязательно должно удовлетворять условию несживаемости и быть стационарным, то есть может быть любым. В таких условиях решение можно представить действительной частью суммы элементарных возмущений вида u(r,\theta,t) exp(i \alpha x), где alpha - действительный параметр, u - комплексная амплитуда возмущения, i - комплексна единица. Alpha не может иметь комплексную составляющею, так как трубе имеет бесконечнуб протяженность и в этом случае амплитуда возмущения нарастала бы неограничено. Нас интересует вопрос устойчивости, если существует хотябы одно alpha, при котором возмущения наростают, профиль скорости неустойчив. Задача сводится к тому, чтобы найти те alpha, при которых возмущения нарастают, или убедиться, что таких нет. 

Задача требет модификаии: поле скорости, давление, завихренность - комплексные функции. Редукция к двумерному случаю, произвожная по x вырождается в умножение на i alpha. Многопроцессорность решино оставить, на каждом процессоре параллельно вычисляется одно сеение трубы но при различных значениях alpha. Введено базовое течение - действительное. Задача для давления решается отдельно для действ и мним частей, но требет модификации для 2D.  

.sap - Symmetry alpha point - контрольная точка комплексного поля скорости (2D) для некоторого значения альфа с условиями симметрии. 
.scs - Symmtry cross section - контрольная точка однородного вдоль трубы течения (2D) с условиями симметрии. 
.scp - Symmtry control point - контрольная точка трехмерного поля скорости (3D) с условиями симметрии. 
.cp - Control point - котрольная точка трехмерного поля скорости (3D) в общем виде. 

step: размерар шага контролируется теперь не из аобсолютной одибки, а из относительной. 

На течени Пуазйля протестировано. 
