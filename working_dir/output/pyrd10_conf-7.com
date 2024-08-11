%nprocs=16
%mem=16GB
%chk=pyrd10_conf-7.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.88070    -2.12900    0.67950
C    4.29460    -1.01240    0.23810
C    3.03640    -1.66930    0.82300
C    1.75270    -0.87380    0.52500
C    0.51440    -1.53460    1.10850
C    0.01370    -1.15100    2.35770
N    -1.08660    -1.71480    2.91820
C    -1.72010    -2.69480    2.22560
C    -1.28720    -3.13870    0.97450
C    -0.15370    -2.54780    0.41190
H    4.46490    -0.01770    0.65250
H    4.24230    -0.92750    -0.84830
H    2.94160    -2.68100    0.42450
H    3.15660    -1.77890    1.90240
H    1.61450    -0.75890    -0.55180
H    1.83470    0.14210    0.91610
H    0.49520    -0.37730    2.93800
H    -2.59230    -3.12550    2.69680
H    -1.81630    -3.92280    0.45280
H    0.20120    -2.87690    -0.55410

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-7.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-7.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

