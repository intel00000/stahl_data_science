%nprocs=16
%mem=16GB
%chk=pyrd10_conf-8.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.08400    -0.12870    0.46850
C    4.44920    -1.95240    -0.01960
C    2.91240    -2.01430    0.01370
C    2.31390    -1.80300    1.42260
C    0.79830    -1.95980    1.54250
C    0.21450    -1.98180    2.81490
N    -1.12030    -2.11450    3.02230
C    -1.92210    -2.22700    1.93450
C    -1.43070    -2.21050    0.62860
C    -0.05450    -2.07420    0.43390
H    4.83170    -2.17330    -1.01710
H    4.90130    -2.66430    0.67260
H    2.51550    -1.26830    -0.67670
H    2.59740    -2.98800    -0.36510
H    2.57220    -0.80720    1.78750
H    2.77650    -2.51050    2.11330
H    0.82070    -1.89040    3.70460
H    -2.98020    -2.32950    2.12940
H    -2.10030    -2.29890    -0.21430
H    0.32830    -2.05580    -0.57490

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-8.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-8.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

