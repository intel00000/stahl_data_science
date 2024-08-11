%nprocs=16
%mem=16GB
%chk=pyrd11_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd11.xyz

0 1
Cl    -3.57010    0.63980    1.36390
C    -2.21290    0.84090    0.30890
C    -1.08900    -0.00640    0.36380
C    -1.06780    -1.18950    1.33680
C    -0.02120    0.25570    -0.53090
C    1.25840    -0.52370    -0.55630
O    1.77760    -0.70780    0.68150
C    3.01080    -1.39430    0.84570
O    1.80740    -0.91310    -1.59090
N    -0.07340    1.27690    -1.43180
C    -1.17010    2.06850    -1.46720
C    -2.25280    1.88620    -0.61270
H    -0.57790    -0.90880    2.26940
H    -0.54260    -2.04490    0.91160
H    -2.06870    -1.55650    1.56220
H    2.94760    -2.41240    0.45800
H    3.81800    -0.87820    0.32330
H    3.27020    -1.44800    1.90300
H    -1.16170    2.86360    -2.19980
H    -3.11240    2.53840    -0.66090

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd11_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd11.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd11_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd11.xyz

0 1

