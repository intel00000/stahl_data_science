%nprocs=16
%mem=16GB
%chk=pyrz1_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz1.xyz

0 1
Cl    3.38730    -2.17850    1.19710
C    2.04380    -1.27140    0.59600
C    2.23300    -0.04370    -0.05610
C    3.61010    0.56390    -0.27540
N    1.17860    0.67350    -0.53080
C    -0.07460    0.16820    -0.35710
C    -1.21130    0.98100    -0.89050
O    -2.41940    0.40480    -0.67460
C    -3.60880    1.04450    -1.11810
O    -1.07010    2.06250    -1.46690
C    -0.25430    -1.06230    0.29710
N    0.79630    -1.78040    0.77210
H    4.23010    -0.11060    -0.86660
H    3.56180    1.52030    -0.79900
H    4.10550    0.73070    0.68150
H    -4.47750    0.43920    -0.85920
H    -3.72480    2.02360    -0.65040
H    -3.59990    1.18090    -2.20090
H    -1.23440    -1.49190    0.45170

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz1_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz1_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz1.xyz

0 1

