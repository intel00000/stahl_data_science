%nprocs=16
%mem=16GB
%chk=pyrd10_conf-12.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    3.18880    -0.42090    2.20700
C    3.96880    -1.44620    0.68790
C    2.97900    -1.54310    -0.49050
C    1.87660    -2.61770    -0.34990
C    0.82980    -2.32240    0.71390
C    0.75750    -3.08240    1.88720
N    -0.16380    -2.85910    2.85870
C    -1.04810    -1.84690    2.67290
C    -1.04810    -1.04070    1.53290
C    -0.09560    -1.28680    0.54180
H    4.87570    -0.91830    0.38990
H    4.26730    -2.42810    1.05710
H    2.54090    -0.56310    -0.68670
H    3.55360    -1.78700    -1.38510
H    2.32640    -3.59530    -0.16800
H    1.35080    -2.72300    -1.30050
H    1.45170    -3.88820    2.07380
H    -1.76490    -1.68770    3.46570
H    -1.76550    -0.24130    1.42120
H    -0.07370    -0.67360    -0.34710

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-12.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-12.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

