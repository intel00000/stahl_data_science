%nprocs=16
%mem=16GB
%chk=pyrmd6_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd6.xyz

0 1
C    -1.98050    0.00050    0.12680
C    -0.46570    -0.00050    0.04270
C    0.27130    1.18820    0.04950
N    1.62660    1.20430    -0.04030
C    2.24730    -0.00070    -0.14310
Cl    3.96810    -0.00090    -0.25750
N    1.61800    -1.20540    -0.16360
C    0.26280    -1.18890    -0.07210
H    -2.40250    -0.84540    -0.41750
H    -2.30900    -0.05220    1.16540
H    -2.39630    0.90100    -0.32700
H    -0.21240    2.15120    0.12840
H    -0.22770    -2.15150    -0.09180

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd6_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd6.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd6_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd6.xyz

0 1

