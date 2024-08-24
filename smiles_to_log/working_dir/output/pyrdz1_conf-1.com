%nprocs=16
%mem=16GB
%chk=pyrdz1_conf-1.chk
# B3LYP gen 6D pseudo=read symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrdz1.xyz

0 1
C    -3.88550    1.81050    0.28310
C    -3.32210    0.47510    0.80550
C    -4.12870    -0.72290    0.26890
C    -1.82840    0.33320    0.51280
N    -1.45690    0.30200    -0.79700
N    -0.16060    0.17640    -1.12600
C    0.79760    0.07520    -0.15970
C    2.21080    -0.04580    -0.63190
C    2.64560    0.68520    -1.74990
C    3.96050    0.59380    -2.20600
C    4.86700    -0.23450    -1.55100
I    6.83880    -0.37530    -2.23010
C    4.45630    -0.97040    -0.44330
C    3.14210    -0.88070    0.01440
F    2.78670    -1.62750    1.08610
C    0.45390    0.11130    1.19700
C    -0.89020    0.24050    1.54300
H    -4.92710    1.94460    0.57610
H    -3.32210    2.65740    0.67680
H    -3.84040    1.86480    -0.80560
H    -3.44480    0.48040    1.89050
H    -4.09410    -0.77280    -0.82030
H    -5.17710    -0.66030    0.56210
H    -3.73650    -1.66580    0.65210
H    1.95600    1.33660    -2.26850
H    4.27740    1.16450    -3.06690
H    5.15380    -1.61840    0.06460
H    1.20900    0.04500    1.96640
H    -1.18810    0.26780    2.58050

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
%chk=pyrdz1_conf-1.chk
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
%chk=pyrdz1_conf-1.chk
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

