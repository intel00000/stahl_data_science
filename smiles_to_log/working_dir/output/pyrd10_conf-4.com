%nprocs=16
%mem=16GB
%chk=pyrd10_conf-4.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    4.77530    -0.74460    -0.84300
C    4.60220    -1.77730    0.85140
C    3.15320    -2.25590    1.04190
C    2.14060    -1.11040    1.23500
C    0.72530    -1.62260    1.44640
C    0.20250    -1.78970    2.73360
N    -1.05340    -2.24880    2.96770
C    -1.82710    -2.55480    1.89550
C    -1.38440    -2.41740    0.57830
C    -0.08930    -1.94460    0.35530
H    5.29120    -2.62090    0.78980
H    4.93260    -1.13920    1.67230
H    2.86200    -2.86770    0.18610
H    3.11070    -2.91670    1.90940
H    2.14380    -0.44270    0.37140
H    2.42870    -0.48470    2.08170
H    0.79110    -1.55620    3.60880
H    -2.82230    -2.91650    2.11200
H    -2.02840    -2.67050    -0.25110
H    0.27720    -1.83020    -0.65500

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-4.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-4.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

