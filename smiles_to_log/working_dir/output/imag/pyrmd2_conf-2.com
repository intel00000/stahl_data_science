%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd2.xyz

0 1
Cl    3.55180    0.83780    -3.14560
C    2.99510    0.39970    -1.57330
C    3.26300    1.24630    -0.49570
C    2.79490    0.85500    0.76020
Cl    3.09950    1.86560    2.12380
N    2.10500    -0.30040    0.94220
C    1.89080    -1.06210    -0.16470
C    1.12050    -2.37360    -0.02200
C    0.63150    -2.74570    1.37640
C    -0.63690    -2.34180    1.80220
C    -1.09200    -2.70070    3.06840
C    -0.27980    -3.46710    3.90210
C    0.98340    -3.87800    3.48100
C    1.43890    -3.51920    2.21480
N    2.30750    -0.76240    -1.42250
H    3.80990    2.16650    -0.62770
H    1.75150    -3.17910    -0.39970
H    0.26350    -2.33750    -0.69590
H    -1.26590    -1.74320    1.15890
H    -2.06920    -2.38420    3.40360
H    -0.63210    -3.74340    4.88560
H    1.60780    -4.47000    4.13460
H    2.42070    -3.83440    1.89180

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd2.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd2.xyz

0 1

