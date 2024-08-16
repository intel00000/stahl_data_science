%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-4.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt=calcfc freq=noraman

 pyrmd1.xyz

0 1
N    -4.32330    1.51380    0.02600
C    -2.99640    1.15670    0.00690
N    -2.72270    -0.13660    -0.24620
C    -1.40920    -0.48710    -0.26540
C    -1.19150    -1.97110    -0.51740
C    -0.37460    0.45270    -0.04480
C    1.09940    0.18550    -0.04760
O    1.44890    -0.98750    -0.62230
C    2.81100    -1.39040    -0.72420
C    3.29080    -2.03750    0.58260
O    1.92490    1.00010    0.37070
C    -0.80020    1.76480    0.21940
N    -2.10550    2.13580    0.24620
H    -4.54990    2.47460    0.21500
H    -5.00830    0.79980    -0.14870
H    -0.54260    -2.39490    0.24820
H    -2.12580    -2.53360    -0.50080
H    -0.73040    -2.12040    -1.49300
H    3.45040    -0.54930    -1.00130
H    2.88730    -2.11510    -1.53530
H    3.26200    -1.32410    1.40690
H    2.66530    -2.88830    0.85250
H    4.31710    -2.39080    0.48670
H    -0.08750    2.55720    0.41250

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-4.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-4.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd1.xyz

0 1

