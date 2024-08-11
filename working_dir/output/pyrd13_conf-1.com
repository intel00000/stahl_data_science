%nprocs=16
%mem=16GB
%chk=pyrd13_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd13.xyz

0 1
O    4.82140    -0.00380    -0.20850
C    3.80610    0.69550    -0.17050
C    2.41910    0.16310    -0.11120
C    2.16800    -1.21370    -0.08790
C    0.84520    -1.65590    -0.03760
C    -0.20760    -0.73350    -0.01180
C    -1.66590    -1.17230    -0.01110
C    -2.56510    -0.05350    0.54220
C    -2.24080    1.30970    -0.10030
N    -0.82860    1.63810    0.10840
C    0.13930    0.63170    -0.02240
N    1.41580    1.07910    -0.07590
H    3.87300    1.78380    -0.17870
H    2.98290    -1.92500    -0.11360
H    0.63940    -2.71650    -0.03160
H    -1.95730    -1.39330    -1.03910
H    -1.79870    -2.09290    0.55930
H    -2.43130    0.02270    1.62310
H    -3.61610    -0.29990    0.38850
H    -2.88800    2.10670    0.27000
H    -2.37340    1.25090    -1.18280
H    -0.53730    2.57900    -0.10870

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd13_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd13.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd13_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd13.xyz

0 1

