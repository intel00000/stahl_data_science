%nprocs=16
%mem=16GB
%chk=pyrmd7_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd7.xyz

0 1
Cl    -2.40180    3.17190    -0.17150
C    -1.07220    2.08050    -0.05400
N    0.17350    2.61550    0.01490
C    1.20920    1.74380    0.10700
C    1.02270    0.35860    0.13120
C    2.19980    -0.59810    0.23500
C    -0.29740    -0.10330    0.05540
O    -0.53650    -1.46530    0.07550
C    -1.88080    -1.93730    -0.00170
N    -1.35380    0.75360    -0.03770
H    2.19600    2.18100    0.16130
H    3.15760    -0.07920    0.28420
H    2.21220    -1.26830    -0.62560
H    2.09650    -1.22100    1.12440
H    -1.88680    -3.02670    0.02670
H    -2.36020    -1.62640    -0.93150
H    -2.47830    -1.57960    0.83860

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd7_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd7.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd7_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd7.xyz

0 1

