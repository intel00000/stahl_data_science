%nprocs=16
%mem=16GB
%chk=pyrd10_conf-6.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.14350    -0.23790    1.42060
C    3.31360    -0.64250    0.75310
C    3.21820    -2.13960    0.41840
C    1.86120    -2.56280    -0.17950
C    0.68670    -2.33500    0.76090
C    0.43100    -3.21850    1.81540
N    -0.59970    -3.05550    2.68350
C    -1.41300    -1.98340    2.50890
C    -1.23090    -1.05430    1.48250
C    -0.16560    -1.23670    0.59810
H    2.60580    -0.35490    1.53260
H    3.13320    -0.01660    -0.12160
H    4.01010    -2.39420    -0.28790
H    3.42140    -2.72190    1.31880
H    1.88950    -3.62200    -0.44190
H    1.68410    -2.04400    -1.12350
H    1.06030    -4.08040    1.98440
H    -2.22480    -1.87990    3.21490
H    -1.89980    -0.21300    1.37420
H    -0.00660    -0.52960    -0.20330

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-6.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-6.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

