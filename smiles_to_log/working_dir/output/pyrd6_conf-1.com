%nprocs=16
%mem=16GB
%chk=pyrd6_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd6.xyz

0 1
Br    4.48030    -2.07190    -3.34020
C    3.64450    -1.15210    -1.89980
C    4.44430    -0.51970    -0.94760
C    3.83040    0.15770    0.10830
C    4.66230    0.84660    1.17170
C    2.43560    0.17610    0.18530
C    1.69200    -0.48150    -0.80700
N    0.29680    -0.53370    -0.84310
C    -0.63360    0.00480    0.02380
O    -1.86120    -0.31790    -0.43910
C    -3.09640    0.05630    0.18210
C    -3.27570    1.58800    0.21660
C    -3.23200    -0.55250    1.59300
C    -4.16710    -0.55490    -0.73440
O    -0.41450    0.66190    1.03900
N    2.28920    -1.13520    -1.83480
H    5.51950    -0.55830    -1.03590
H    4.20760    0.71430    2.15430
H    5.66670    0.42630    1.23020
H    4.73110    1.91570    0.96950
H    1.95700    0.69440    1.00180
H    -0.09460    -1.04040    -1.61680
H    -4.27150    1.86590    0.56120
H    -2.55810    2.06500    0.88420
H    -3.13170    2.02330    -0.77230
H    -4.22560    -0.38290    2.00730
H    -2.51350    -0.12220    2.29060
H    -3.05710    -1.62820    1.57580
H    -4.09060    -0.16320    -1.74930
H    -4.06050    -1.63820    -0.80090
H    -5.17380    -0.34340    -0.37360

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd6_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd6.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd6_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd6.xyz

0 1

