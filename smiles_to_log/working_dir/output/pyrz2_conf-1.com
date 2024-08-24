%nprocs=16
%mem=16GB
%chk=pyrz2_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz2.xyz

0 1
Cl    -4.20320    0.97240    -0.32240
C    -2.49300    0.76560    -0.17900
C    -1.91370    -0.50880    -0.20840
C    -2.73380    -1.78120    -0.36600
N    -0.56930    -0.67010    -0.09560
C    0.21700    0.43580    0.04810
C    1.70050    0.21320    0.16970
O    2.06180    -1.09520    0.11810
C    3.42890    -1.46950    0.21920
O    2.52150    1.12600    0.30420
C    -0.38300    1.71140    0.07570
N    0.33890    2.88910    0.21820
N    -1.72320    1.87290    -0.03670
H    -2.11240    -2.67840    -0.36680
H    -3.28960    -1.75900    -1.30410
H    -3.45100    -1.87240    0.45070
H    3.52440    -2.55350    0.15770
H    3.85450    -1.14320    1.16950
H    4.01570    -1.02940    -0.58870
H    -0.13770    3.77560    0.23190
H    1.34670    2.79870    0.30470

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz2.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz2.xyz

0 1

