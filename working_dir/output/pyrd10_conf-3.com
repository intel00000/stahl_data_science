%nprocs=16
%mem=16GB
%chk=pyrd10_conf-3.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.07460    -0.09290    1.34640
C    3.62150    -1.45180    1.33670
C    2.96490    -1.47710    -0.05300
C    1.85550    -2.53820    -0.20030
C    0.68540    -2.32870    0.75020
C    0.48750    -3.17290    1.84890
N    -0.53690    -3.01780    2.72570
C    -1.40380    -1.99550    2.51450
C    -1.28210    -1.10910    1.44250
C    -0.22190    -1.28220    0.55030
H    4.06460    -2.41400    1.59640
H    2.90930    -1.18100    2.11820
H    2.56110    -0.48830    -0.27660
H    3.73440    -1.66270    -0.80400
H    2.27180    -3.53790    -0.06280
H    1.46980    -2.52970    -1.22140
H    1.15720    -3.99680    2.04860
H    -2.20950    -1.89680    3.22820
H    -1.99180    -0.30630    1.30580
H    -0.10630    -0.60740    -0.28570

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-3.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-3.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

