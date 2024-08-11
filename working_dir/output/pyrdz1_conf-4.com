%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-4.chk
# B3LYP gen 6D pseudo=read symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrdz1.xyz

0 1
C    -3.70840    -0.00490    2.01170
C    -3.37890    0.21100    0.52300
C    -4.07970    1.46040    -0.04190
C    -1.87160    0.24280    0.27380
N    -1.36810    -0.74510    -0.51520
N    -0.05140    -0.79080    -0.77990
C    0.79210    0.14870    -0.26540
C    2.23430    0.02410    -0.62180
C    2.95000    1.16190    -1.02310
C    4.30020    1.07870    -1.36430
C    4.95800    -0.14690    -1.30640
I    6.98100    -0.27480    -1.81610
C    4.26500    -1.28660    -0.90730
C    2.91550    -1.20650    -0.56600
F    2.28710    -2.33560    -0.16290
C    0.31350    1.18050    0.55010
C    -1.05170    1.22930    0.82790
H    -4.78090    -0.12860    2.16550
H    -3.21900    -0.90140    2.39540
H    -3.38290    0.83590    2.62480
H    -3.78630    -0.64880    -0.01640
H    -3.77120    2.36770    0.47800
H    -5.16340    1.38100    0.05000
H    -3.85160    1.59290    -1.10050
H    2.45090    2.11800    -1.08680
H    4.83760    1.96200    -1.67790
H    4.76990    -2.23890    -0.85670
H    0.98310    1.91730    0.96760
H    -1.46140    2.00520    1.45680

C H N F 0
6-31G(d,p)
****
I 0
LANL2DZ
****

I 0
LANL2DZ

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-4.chk
# M062X gen pseudo=read int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrdz1.xyz

0 1

C H N F 0
def2tzvp
****
I 0
SDD
****

I 0
SDD

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-4.chk
# M062X gen pseudo=read int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrdz1.xyz

0 1

C H N F 0
def2tzvp
****
I 0
SDD
****

I 0
SDD

