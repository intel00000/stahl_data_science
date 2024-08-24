%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd2.xyz

0 1
Cl    2.54790    2.90940    1.27280
C    2.44800    1.39490    0.45400
C    3.49680    1.00130    -0.38040
C    3.37780    -0.23290    -1.02370
Cl    4.64680    -0.76520    -2.06300
N    2.29110    -1.03090    -0.85630
C    1.32180    -0.55980    -0.02930
C    0.11840    -1.45830    0.20340
C    0.29370    -2.41730    1.37300
C    -0.11610    -2.03720    2.65380
C    0.02360    -2.92290    3.71940
C    0.56770    -4.18690    3.49950
C    0.97370    -4.57330    2.22360
C    0.83410    -3.68800    1.15770
N    1.34780    0.62050    0.64280
H    4.36550    1.62500    -0.52080
H    -0.09300    -2.02310    -0.70580
H    -0.76770    -0.84220    0.36330
H    -0.53310    -1.05510    2.82610
H    -0.28650    -2.63000    4.71220
H    0.67780    -4.87330    4.32680
H    1.39690    -5.55430    2.06190
H    1.15420    -3.98630    0.16940

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd2.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd2.xyz

0 1

