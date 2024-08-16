%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd1.xyz

0 1
N    -4.00580    2.10370    0.07620
C    -2.82480    1.40160    0.04710
N    -2.86290    0.18000    -0.51590
C    -1.69280    -0.51140    -0.54270
C    -1.81470    -1.88120    -1.19300
C    -0.49130    0.01090    -0.01170
C    0.81950    -0.72100    -0.02410
O    1.82530    -0.01950    0.56120
C    3.14940    -0.53900    0.63530
C    3.91220    -0.28740    -0.67300
O    0.98120    -1.84250    -0.50650
C    -0.59330    1.29820    0.54640
N    -1.75010    2.00690    0.58440
H    -4.00240    3.01840    0.49220
H    -4.82170    1.67030    -0.31930
H    -1.53150    -2.65910    -0.48450
H    -2.82950    -2.09210    -1.53150
H    -1.15480    -1.94550    -2.05760
H    3.14060    -1.60260    0.88370
H    3.66280    -0.03550    1.45500
H    3.44350    -0.81050    -1.50720
H    3.93720    0.77510    -0.91460
H    4.94040    -0.63990    -0.59680
H    0.26110    1.79910    0.98160

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd1.xyz

0 1

