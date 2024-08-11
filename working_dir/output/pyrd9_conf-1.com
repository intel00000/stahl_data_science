%nprocs=16
%mem=16GB
%chk=pyrd9_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd9.xyz

0 1
C    2.23810    -1.14300    -0.01720
C    0.76530    -0.78490    0.01700
N    -0.08480    -1.84000    0.11010
C    -1.42400    -1.62350    0.14760
Cl    -2.45860    -2.99740    0.26610
C    -2.00940    -0.35670    0.09580
C    -1.12180    0.71810    0.00020
C    0.27920    0.53100    -0.04130
C    0.81480    1.85120    -0.13810
N    -0.13620    2.74970    -0.15480
N    -1.30380    2.08320    -0.07240
H    2.68850    -0.78170    -0.94140
H    2.40290    -2.22030    0.03940
H    2.75440    -0.67900    0.82290
H    -3.07990    -0.23240    0.12860
H    1.85370    2.14150    -0.19380
H    -2.17840    2.58420    -0.06890

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd9_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd9.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd9_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd9.xyz

0 1

