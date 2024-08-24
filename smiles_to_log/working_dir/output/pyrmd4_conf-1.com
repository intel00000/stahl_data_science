%nprocs=16
%mem=16GB
%chk=pyrmd4_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd4.xyz

0 1
Cl    0.61560    -1.70120    -3.04860
C    0.45410    -0.69120    -1.66030
N    1.06680    0.52260    -1.65930
C    0.90450    1.27410    -0.53930
C    1.55390    2.57600    -0.52220
N    2.07060    3.61260    -0.50800
N    0.19520    0.92930    0.56680
C    -0.40330    -0.29040    0.53050
C    -1.19420    -0.67470    1.76390
C    -0.30040    -1.14370    -0.57320
H    -0.54580    -0.68540    2.64080
H    -1.64520    -1.66220    1.66490
H    -1.99340    0.04590    1.94040
H    -0.77810    -2.11170    -0.59640

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd4_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd4.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd4_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd4.xyz

0 1

