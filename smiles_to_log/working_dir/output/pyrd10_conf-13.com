%nprocs=16
%mem=16GB
%chk=pyrd10_conf-13.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    3.51830    -3.02080    1.77770
C    3.87650    -1.87640    0.18750
C    3.07530    -0.56040    0.24920
C    1.59390    -0.65680    -0.18160
C    0.69510    -1.40120    0.79380
C    0.27480    -0.79780    1.98450
N    -0.53250    -1.41610    2.88370
C    -0.94170    -2.67920    2.60390
C    -0.57150    -3.35500    1.43960
C    0.25660    -2.70260    0.52370
H    4.94880    -1.67650    0.18440
H    3.64850    -2.45560    -0.70810
H    3.15990    -0.12150    1.24480
H    3.55900    0.15060    -0.42210
H    1.51440    -1.10420    -1.17390
H    1.18640    0.34920    -0.29710
H    0.58780    0.20360    2.24060
H    -1.57690    -3.14960    3.34070
H    -0.91190    -4.36280    1.25320
H    0.56250    -3.21010    -0.37920

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-13.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-13.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

