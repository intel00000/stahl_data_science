%nprocs=16
%mem=16GB
%chk=pyrd10_conf-5.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd10.xyz

0 1
Br    5.40840    -1.22240    1.30130
C    4.19140    -1.06720    -0.26790
C    2.72520    -1.01440    0.19320
C    2.24190    -2.30760    0.87760
C    0.78140    -2.22890    1.29060
C    0.41900    -1.77390    2.56320
N    -0.87120    -1.68130    2.97560
C    -1.84380    -2.05090    2.10440
C    -1.57000    -2.51560    0.81630
C    -0.23740    -2.60390    0.40770
H    4.46860    -0.15940    -0.80560
H    4.38110    -1.91450    -0.92840
H    2.09080    -0.81110    -0.67120
H    2.59110    -0.16670    0.86780
H    2.84710    -2.52450    1.75980
H    2.38000    -3.16620    0.21820
H    1.16920    -1.47050    3.27950
H    -2.86070    -1.96750    2.46080
H    -2.37090    -2.80120    0.15020
H    -0.00110    -2.96010    -0.58460

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-5.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd10.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd10_conf-5.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd10.xyz

0 1

