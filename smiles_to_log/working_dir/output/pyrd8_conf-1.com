%nprocs=16
%mem=16GB
%chk=pyrd8_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd8.xyz

0 1
Br    -1.12580    -4.05330    -0.88780
C    -0.91220    -2.21110    -0.46290
N    0.34270    -1.70900    -0.34520
C    0.49710    -0.39040    -0.04340
C    1.89900    0.12980    0.07630
C    2.36730    1.41590    0.37450
N    3.67220    1.59400    0.41480
C    4.30660    0.48340    0.15650
S    3.30550    -0.88930    -0.16130
C    -0.60280    0.45360    0.14820
C    -1.89300    -0.07100    0.03040
C    -3.10470    0.82090    0.21770
C    -2.04960    -1.42310    -0.28410
H    1.73800    2.26980    0.57340
H    5.38180    0.40380    0.13720
H    -0.45680    1.49520    0.38520
H    -3.37700    0.88000    1.27160
H    -3.95900    0.45000    -0.34970
H    -2.90220    1.82870    -0.14650
H    -3.02960    -1.86330    -0.38830

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd8_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd8.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd8_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd8.xyz

0 1

