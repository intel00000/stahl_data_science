%nprocs=16
%mem=16GB
%chk=pyrd2_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd2.xyz

0 1
Cl    -1.40520    -2.49820    -0.41520
C    -0.19580    -1.27520    -0.23470
N    1.09970    -1.66250    -0.36510
C    2.07140    -0.72660    -0.22900
C    1.79500    0.61410    0.03880
C    0.46210    1.00590    0.17170
C    -0.55730    0.05630    0.03510
C    -2.01250    0.47690    0.17880
H    3.08860    -1.07530    -0.34000
H    2.59420    1.33390    0.14120
H    0.22710    2.04010    0.37920
H    -2.48570    -0.08000    0.98840
H    -2.12050    1.54180    0.38670
H    -2.56100    0.24890    -0.73590

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd2.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd2_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd2.xyz

0 1

