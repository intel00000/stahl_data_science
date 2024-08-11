%nprocs=16
%mem=16GB
%chk=pyrd8_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd8.xyz

0 1
Br    -1.07420    -4.05240    -0.88730
C    -0.89160    -2.20700    -0.46190
N    0.35480    -1.68520    -0.33970
C    0.48820    -0.36530    -0.03760
C    1.90400    0.11160    0.07200
C    3.10800    -0.58450    -0.08970
N    4.22790    0.09120    0.06630
C    3.98390    1.34020    0.35550
S    2.31140    1.76930    0.45560
C    -0.62670    0.46230    0.15000
C    -1.90780    -0.08280    0.02740
C    -3.13380    0.79020    0.21040
C    -2.04190    -1.43720    -0.28740
H    3.16900    -1.63580    -0.33290
H    4.75920    2.07030    0.52420
H    -0.50190    1.50690    0.38760
H    -3.40250    0.85610    1.26480
H    -3.98410    0.40040    -0.35010
H    -2.94830    1.79760    -0.16400
H    -3.01450    -1.89260    -0.39510

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd8_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd8.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd8_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd8.xyz

0 1

