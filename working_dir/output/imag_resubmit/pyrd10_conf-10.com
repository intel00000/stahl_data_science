%nprocs=16
%mem=16GB
%chk=pyrd10_conf-10.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt=calcfc freq=noraman

 pyrd10.xyz

0 1
Br    5.58380    -1.83980    -0.78950
C    4.46410    -1.08710    0.67250
C    3.04830    -1.67240    0.57350
C    2.10340    -1.14300    1.67320
C    0.67110    -1.67580    1.64840
C    -0.23570    -1.23460    2.61970
N    -1.52650    -1.64940    2.68010
C    -1.95070    -2.53730    1.74720
C    -1.11910    -3.03230    0.74200
C    0.20620    -2.59460    0.69480
H    4.94400    -1.33840    1.61960
H    4.46500    -0.00140    0.56440
H    2.63810    -1.43570    -0.40990
H    3.11250    -2.76020    0.63540
H    2.05190    -0.05440    1.61010
H    2.52530    -1.37560    2.65280
H    0.06470    -0.52760    3.37950
H    -2.98270    -2.85030    1.81920
H    -1.49170    -3.73920    0.01510
H    0.85140    -2.97630    -0.08210

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-10.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-10.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

