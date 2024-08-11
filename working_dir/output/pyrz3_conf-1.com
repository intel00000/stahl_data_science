%nprocs=16
%mem=16GB
%chk=pyrz3_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrz3.xyz

0 1
C    -3.76740    -0.46570    0.42790
C    -2.29670    -0.19890    0.18230
N    -1.46410    -0.38670    1.23920
C    -0.14070    -0.15200    1.03690
C    0.35720    0.26730    -0.20760
C    1.82070    0.51620    -0.39370
O    2.15020    0.91560    -1.64770
C    3.50330    1.19160    -1.98160
O    2.64710    0.37080    0.51090
N    -0.49100    0.45210    -1.26130
C    -1.81440    0.21790    -1.06100
H    -4.14310    0.18620    1.21720
H    -4.36460    -0.29190    -0.46740
H    -3.91710    -1.50000    0.73920
H    0.51310    -0.30460    1.88510
H    3.90270    1.99720    -1.36330
H    3.57460    1.49680    -3.02550
H    4.12910    0.30870    -1.84150
H    -2.46690    0.37010    -1.90700

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrz3.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrz3_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrz3.xyz

0 1

