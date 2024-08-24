%nprocs=16
%mem=16GB
%chk=pyrd10_conf-9.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.37650    -2.28450    0.29370
C    4.35150    -0.58580    0.46980
C    3.02450    -0.83130    1.20840
C    2.04780    -1.74960    0.43950
C    0.67460    -1.95950    1.07670
C    0.33630    -1.48660    2.35370
N    -0.88390    -1.67570    2.91770
C    -1.81380    -2.35930    2.20550
C    -1.56450    -2.86930    0.93090
C    -0.30420    -2.66450    0.36600
H    4.98180    0.11870    1.01430
H    4.19260    -0.18150    -0.53100
H    3.24120    -1.25530    2.19040
H    2.54480    0.13300    1.38420
H    2.50050    -2.73320    0.30190
H    1.89010    -1.34170    -0.56060
H    1.03990    -0.94280    2.96500
H    -2.77500    -2.49550    2.68060
H    -2.32810    -3.41170    0.39290
H    -0.09050    -3.05390    -0.61920

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-9.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-9.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

