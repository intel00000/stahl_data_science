%nprocs=16
%mem=16GB
%chk=pyrz1_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz1.xyz

0 1
Cl    3.62090    -1.92270    1.07700
C    2.12840    -1.21210    0.57060
C    2.10430    0.01980    -0.10100
C    3.36890    0.80100    -0.42380
N    0.93230    0.58310    -0.50150
C    -0.22650    -0.08420    -0.23220
C    -1.53030    0.51110    -0.66070
O    -1.40230    1.70180    -1.29770
C    -2.55170    2.39350    -1.76630
O    -2.61710    -0.03380    -0.45200
C    -0.19240    -1.31520    0.43970
N    0.97600    -1.88000    0.84080
H    4.02940    0.20500    -1.05430
H    3.15700    1.73440    -0.94790
H    3.90500    1.04610    0.49370
H    -2.25530    3.32490    -2.24860
H    -3.10070    1.79330    -2.49370
H    -3.22560    2.63620    -0.94300
H    -1.09790    -1.86340    0.66470

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz1_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz1.xyz

0 1

