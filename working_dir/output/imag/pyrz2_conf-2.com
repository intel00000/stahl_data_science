%nprocs=16
%mem=16GB
%chk=pyrz2_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz2.xyz

0 1
Cl    -4.15180    1.29010    -0.29540
C    -2.48290    0.85750    -0.17140
C    -2.07660    -0.48040    -0.22270
C    -3.05660    -1.63290    -0.38860
N    -0.76430    -0.81710    -0.12500
C    0.16410    0.17000    0.02520
C    1.58060    -0.32790    0.12190
O    2.51140    0.64620    0.27640
C    3.89000    0.31240    0.38130
O    1.88190    -1.52360    0.06630
C    -0.26320    1.51550    0.07510
N    0.58900    2.60150    0.22460
N    -1.57230    1.85060    -0.02250
H    -2.55600    -2.60280    -0.40710
H    -3.61210    -1.52630    -1.32100
H    -3.77240    -1.64260    0.43410
H    4.48750    1.21610    0.50000
H    4.23510    -0.20890    -0.51330
H    4.07170    -0.33170    1.24330
H    0.20090    3.53080    0.25080
H    1.57960    2.41600    0.30370

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz2_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz2.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz2_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz2.xyz

0 1

