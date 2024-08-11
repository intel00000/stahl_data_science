%nprocs=16
%mem=16GB
%chk=pyrz3_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz3.xyz

0 1
C    -3.64660    -0.56360    0.83080
C    -2.25190    -0.23710    0.33820
N    -1.24710    -0.35900    1.24400
C    0.01130    -0.07000    0.81920
C    0.27460    0.33930    -0.50080
C    1.64960    0.66310    -0.99310
O    2.60070    0.51980    -0.03730
C    3.96410    0.78730    -0.33650
O    1.89850    1.02330    -2.14670
N    -0.74520    0.45680    -1.39820
C    -2.00380    0.16830    -0.97540
H    -3.91800    0.09120    1.65940
H    -4.39160    -0.44050    0.04440
H    -3.69200    -1.59520    1.18140
H    0.79730    -0.17270    1.55390
H    4.32570    0.13430    -1.13250
H    4.58020    0.61970    0.54690
H    4.09930    1.82250    -0.65440
H    -2.79360    0.26880    -1.70400

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz3_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz3.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz3_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz3.xyz

0 1

