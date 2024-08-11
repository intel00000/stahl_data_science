%nprocs=16
%mem=16GB
%chk=pyrd1_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd1.xyz

0 1
Br    1.84980    2.93050    -0.57160
C    1.69530    1.06420    -0.21530
C    0.44240    0.45320    -0.04530
C    0.42590    -0.93600    0.19450
O    -0.72720    -1.66310    0.42020
C    -1.91920    -1.09350    -0.12390
C    -2.06060    0.36610    0.34240
C    -0.87310    1.22110    -0.14340
C    1.62790    -1.63760    0.27540
N    2.83800    -1.04610    0.11400
C    2.85390    0.28870    -0.12860
H    -1.88680    -1.13490    -1.21480
H    -2.77510    -1.68670    0.20030
H    -3.00340    0.79140    -0.00280
H    -2.09790    0.37800    1.43300
H    -1.01960    1.50930    -1.18570
H    -0.81750    2.14890    0.42850
H    1.62240    -2.69720    0.48000
H    3.82510    0.74360    -0.25690

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd1.xyz

0 1

