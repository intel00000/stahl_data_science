%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd1.xyz

0 1
N    -3.92490    2.18450    0.21880
C    -2.77890    1.43550    0.09970
N    -2.92690    0.17200    -0.33880
C    -1.79110    -0.56580    -0.45530
C    -2.03790    -1.98100    -0.95580
C    -0.51490    -0.04810    -0.13690
C    0.76100    -0.83100    -0.25240
O    1.85420    -0.11850    0.11030
C    3.14580    -0.70670    0.06150
C    4.17080    0.33080    0.52500
O    0.82420    -2.00000    -0.63540
C    -0.50450    1.28600    0.30950
N    -1.62510    2.04140    0.43360
H    -3.83970    3.13130    0.54450
H    -4.79840    1.75170    -0.02540
H    -1.68880    -2.70490    -0.22010
H    -3.09410    -2.17800    -1.14200
H    -1.49850    -2.14850    -1.88750
H    3.37170    -1.03380    -0.95540
H    3.18190    -1.58810    0.70500
H    4.15540    1.21090    -0.11820
H    3.96530    0.65560    1.54520
H    5.17870    -0.08290    0.50210
H    0.41460    1.78760    0.58180

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd1.xyz

0 1

