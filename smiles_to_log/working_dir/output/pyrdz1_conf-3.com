%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-3.chk
# B3LYP gen 6D pseudo=read symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrdz1.xyz

0 1
C    -4.12700    0.93330    -0.44580
C    -3.29610    0.58380    0.80360
C    -3.80270    -0.70920    1.47080
C    -1.80080    0.51580    0.49270
N    -1.39830    -0.43820    -0.39220
N    -0.09820    -0.55830    -0.70930
C    0.82980    0.27120    -0.15320
C    2.24640    0.07160    -0.57240
C    3.03270    1.18310    -0.91010
C    4.36120    1.03000    -1.30750
C    4.92610    -0.24090    -1.37070
I    6.91650    -0.47370    -1.96480
C    4.16230    -1.35570    -1.03630
C    2.83410    -1.20580    -0.63900
F    2.13560    -2.31500    -0.30210
C    0.45670    1.26420    0.75970
C    -0.89130    1.39010    1.09240
H    -5.18010    1.06950    -0.19810
H    -3.77610    1.85680    -0.90780
H    -4.06560    0.14660    -1.19930
H    -3.44600    1.39220    1.52210
H    -3.72800    -1.56300    0.79570
H    -4.84650    -0.61910    1.77270
H    -3.22190    -0.94660    2.36300
H    2.60530    2.17500    -0.87980
H    4.95340    1.89430    -1.57100
H    4.59520    -2.34320    -1.07950
H    1.19380    1.91270    1.20910
H    -1.21540    2.14100    1.79730

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
%chk=pyrdz1_conf-3.chk
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
%chk=pyrdz1_conf-3.chk
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

