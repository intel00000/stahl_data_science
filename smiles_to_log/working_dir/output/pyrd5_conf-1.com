%nprocs=16
%mem=16GB
%chk=pyrd5_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd5.xyz

0 1
Cl    0.21300    3.65240    -0.61010
C    0.06690    1.95760    -0.32560
N    -1.17390    1.40970    -0.28710
C    -1.28970    0.07090    -0.06480
C    -2.71810    -0.48270    -0.03620
F    -2.73300    -1.79940    0.18630
F    -3.42000    0.10660    0.93320
F    -3.32240    -0.25550    -1.20360
C    -0.16840    -0.74950    0.12530
C    1.10490    -0.17420    0.08870
C    2.34080    -1.03280    0.27550
C    1.22390    1.19850    -0.14440
H    -0.28500    -1.80990    0.29920
H    3.19880    -0.60840    -0.24730
H    2.58430    -1.12980    1.33380
H    2.18690    -2.02990    -0.13880
H    2.19110    1.67670    -0.18430

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd5_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd5.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd5_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd5.xyz

0 1

