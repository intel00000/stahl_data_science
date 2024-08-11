%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd3.xyz

0 1
Cl    4.28920    1.26210    -0.58420
C    2.69570    0.65920    -0.31670
N    2.14570    0.86560    0.90900
C    0.89050    0.38470    1.10450
C    0.18320    -0.29300    0.10570
C    -1.22200    -0.82080    0.34020
C    -2.29910    0.21160    -0.02130
C    0.83390    -0.45370    -1.12230
N    2.08830    0.01560    -1.34840
H    0.45990    0.55680    2.08040
H    -1.36600    -1.73460    -0.23950
H    -1.32480    -1.12310    1.38410
H    -3.29890    -0.18680    0.15410
H    -2.19560    1.11870    0.57560
H    -2.23750    0.49820    -1.07200
H    0.35750    -0.96050    -1.94910

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd3.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd3.xyz

0 1

