%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd3.xyz

0 1
Cl    4.37380    0.96510    -0.69090
C    2.72840    0.59690    -0.32910
N    2.33510    0.69610    0.96730
C    1.03650    0.40270    1.23520
C    0.12760    0.01500    0.24390
C    -1.31970    -0.30090    0.61450
C    -2.23680    -0.71720    -0.54970
C    0.63380    -0.05660    -1.06210
N    1.92690    0.23030    -1.36270
H    0.73650    0.48540    2.27000
H    -1.31240    -1.09540    1.36300
H    -1.74380    0.57770    1.10440
H    -3.24650    -0.92220    -0.19220
H    -2.31260    0.07020    -1.30060
H    -1.87620    -1.62300    -1.03890
H    0.01730    -0.34520    -1.90020

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd3.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd3.xyz

0 1

