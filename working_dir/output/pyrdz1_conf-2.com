%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-2.chk
# B3LYP gen 6D pseudo=read symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrdz1.xyz

0 1
C    -4.01330    -0.74770    1.06810
C    -3.40620    0.47440    0.35470
C    -3.77910    1.79120    1.06060
C    -1.89510    0.33430    0.17490
N    -1.43930    0.29080    -1.10630
N    -0.12360    0.16610    -1.34820
C    0.76980    0.07790    -0.32030
C    2.21110    -0.04350    -0.69790
C    2.71640    0.67800    -1.79220
C    4.05860    0.58590    -2.16040
C    4.92220    -0.23320    -1.43930
I    6.93460    -0.37470    -1.98650
C    4.44150    -0.95950    -0.35360
C    3.10000    -0.86920    0.01620
F    2.67700    -1.60660    1.06970
C    0.33820    0.12680    1.01070
C    -1.02600    0.25530    1.26580
H    -5.10150    -0.68620    1.10150
H    -3.75500    -1.67260    0.55020
H    -3.65730    -0.83530    2.09510
H    -3.85710    0.51270    -0.64090
H    -3.41220    1.81810    2.08710
H    -4.86040    1.92920    1.09430
H    -3.35720    2.65030    0.53710
H    2.06080    1.32230    -2.36120
H    4.42990    1.14920    -3.00430
H    5.10580    -1.60050    0.20520
H    1.04100    0.07080    1.82880
H    -1.39660    0.29290    2.27910

C H N F 0
6-31G(d,p)
****
I 0
LANL2DZ
****

I 0
LANL2DZ

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-2.chk
# M062X gen pseudo=read int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrdz1.xyz

0 1

C H N F 0
def2tzvp
****
I 0
SDD
****

I 0
SDD

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-2.chk
# M062X gen pseudo=read int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrdz1.xyz

0 1

C H N F 0
def2tzvp
****
I 0
SDD
****

I 0
SDD

