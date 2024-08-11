%nprocs=16
%mem=16GB
%chk=pyrd10_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    4.17620    -2.62250    -0.65250
C    3.78260    -2.09220    1.22720
C    3.20590    -0.66450    1.27550
C    1.86230    -0.48010    0.53780
C    0.73710    -1.33290    1.10510
C    -0.00190    -0.90550    2.21360
N    -1.00920    -1.63570    2.75690
C    -1.30440    -2.83170    2.18720
C    -0.61870    -3.33230    1.07850
C    0.41510    -2.56870    0.53270
H    4.71580    -2.15980    1.78770
H    3.09400    -2.82550    1.64980
H    3.08170    -0.36880    2.31840
H    3.94070    0.02540    0.85790
H    1.97750    -0.69090    -0.52700
H    1.55890    0.56730    0.58080
H    0.20710    0.04110    2.69020
H    -2.11160    -3.39030    2.63940
H    -0.88090    -4.28880    0.65060
H    0.96410    -2.93560    -0.32360

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

