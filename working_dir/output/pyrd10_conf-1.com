%nprocs=16
%mem=16GB
%chk=pyrd10_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    4.63420    -1.08840    1.87510
C    3.67620    -2.40540    0.72770
C    2.89410    -1.68140    -0.38480
C    1.78360    -0.73360    0.11650
C    0.69410    -1.44350    0.90630
C    0.70770    -1.44550    2.30570
N    -0.23700    -2.07180    3.05290
C    -1.23650    -2.71770    2.40070
C    -1.32650    -2.76320    1.00790
C    -0.34520    -2.11580    0.25400
H    4.42470    -3.07960    0.30970
H    3.02170    -2.99980    1.36740
H    2.45300    -2.42560    -1.04950
H    3.59740    -1.11620    -0.99830
H    2.21070    0.06820    0.72180
H    1.32050    -0.22410    -0.73030
H    1.49360    -0.94480    2.85400
H    -1.97500    -3.20440    3.02180
H    -2.13780    -3.28660    0.52360
H    -0.39310    -2.13720    -0.82520

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

