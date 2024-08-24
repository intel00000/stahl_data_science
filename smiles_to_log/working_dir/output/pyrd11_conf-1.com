%nprocs=16
%mem=16GB
%chk=pyrd11_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd11.xyz

0 1
Cl    -3.88890    0.56350    0.84090
C    -2.31630    0.83080    0.16870
C    -1.26350    -0.09390    0.31940
C    -1.47500    -1.37370    1.13550
C    -0.01510    0.21850    -0.27830
C    1.16440    -0.71930    -0.26380
O    2.35640    -0.06790    -0.21630
C    3.57440    -0.79680    -0.15580
O    1.07190    -1.94810    -0.18810
N    0.16410    1.37060    -0.98940
C    -0.86260    2.24240    -1.11120
C    -2.11290    2.00700    -0.55060
H    -0.62000    -1.55770    1.78650
H    -1.57910    -2.23380    0.47330
H    -2.33860    -1.33060    1.79710
H    3.61500    -1.41750    0.74070
H    3.68660    -1.44400    -1.02710
H    4.41930    -0.10870    -0.13200
H    -0.66370    3.14040    -1.67940
H    -2.91650    2.71880    -0.66980

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd11_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd11.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd11_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd11.xyz

0 1

