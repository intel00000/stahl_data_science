%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd5.xyz

0 1
O    -4.57930    -2.14320    0.20040
C    -4.09000    -1.01840    0.10310
C    -2.62830    -0.75250    0.06040
C    -1.71740    -1.81440    0.13220
N    -0.37310    -1.63300    0.09800
C    -2.07420    0.53520    -0.05100
C    -2.91180    1.80270    -0.14030
N    -0.73150    0.74110    -0.08710
C    0.07880    -0.35230    -0.01130
C    1.56810    -0.13060    -0.05080
C    2.10840    1.16250    -0.16250
C    3.48510    1.37090    -0.19940
C    4.34890    0.28330    -0.12460
C    3.84410    -1.00800    -0.01360
C    2.46650    -1.20960    0.02270
H    -4.72870    -0.13760    0.04180
H    -2.05480    -2.83960    0.21930
H    -3.53170    1.91150    0.74940
H    -2.29200    2.69720    -0.22110
H    -3.55750    1.76990    -1.01770
H    1.45830    2.02290    -0.22210
H    3.87880    2.37330    -0.28560
H    5.41750    0.44220    -0.15300
H    4.51620    -1.85190    0.04440
H    2.09930    -2.22180    0.10940

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd5.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd5.xyz

0 1

