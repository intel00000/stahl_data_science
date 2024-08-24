%nprocs=16
%mem=16GB
%chk=pyrd3_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd3.xyz

0 1
Br    -2.53810    2.43790    -0.10310
C    -1.41700    0.90190    -0.05420
C    -0.02970    1.04640    -0.00880
C    0.75850    -0.10710    0.02620
C    2.27020    -0.03490    0.07630
N    0.23180    -1.35850    0.01790
C    -1.12000    -1.47900    -0.02650
C    -1.97670    -0.37640    -0.06340
H    0.41340    2.03040    -0.00130
H    2.69950    -0.53940    -0.78970
H    2.62640    0.99520    0.08060
H    2.64040    -0.52630    0.97630
H    -1.51040    -2.48680    -0.03210
H    -3.04830    -0.50350    -0.09830

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd3.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd3.xyz

0 1

