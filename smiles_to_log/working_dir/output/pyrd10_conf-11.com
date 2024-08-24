%nprocs=16
%mem=16GB
%chk=pyrd10_conf-11.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.92620    -0.98540    1.25630
C    4.34710    -1.29100    0.08560
C    3.10360    -1.49700    0.96210
C    1.82350    -1.74240    0.13520
C    0.53620    -1.95780    0.93010
C    0.48320    -1.94030    2.33240
N    -0.66200    -2.13200    3.03540
C    -1.80570    -2.35030    2.34020
C    -1.84720    -2.38460    0.94590
C    -0.65990    -2.18550    0.23880
H    4.24650    -0.42270    -0.56730
H    4.55860    -2.16390    -0.53370
H    3.28000    -2.34310    1.62870
H    2.97080    -0.61800    1.59550
H    1.96880    -2.61750    -0.50100
H    1.66020    -0.89570    -0.53420
H    1.36280    -1.77100    2.93480
H    -2.70080    -2.49940    2.92720
H    -2.77670    -2.56120    0.42480
H    -0.67230    -2.20850    -0.84150

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-11.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-11.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

