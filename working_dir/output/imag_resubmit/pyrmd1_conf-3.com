%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-3.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt=calcfc freq=noraman

 pyrmd1.xyz

0 1
N    -4.33680    1.46690    -0.24070
C    -3.01830    1.11500    -0.07880
N    -2.60370    0.01130    -0.72760
C    -1.29880    -0.33490    -0.56650
C    -0.92560    -1.60660    -1.31320
C    -0.41030    0.42330    0.23360
C    1.04250    0.14750    0.48340
O    1.56930    -0.81300    -0.30080
C    2.94070    -1.16380    -0.19100
C    3.23570    -2.28200    -1.19320
O    1.70920    0.78030    1.30550
C    -0.97620    1.54650    0.85940
N    -2.27600    1.90870    0.71370
H    -4.91240    0.88880    -0.82730
H    -4.66750    2.28720    0.23660
H    -0.48930    -2.33290    -0.62860
H    -1.78920    -2.07860    -1.78340
H    -0.20350    -1.38220    -2.09740
H    3.16440    -1.49670    0.82460
H    3.56880    -0.29480    -0.39710
H    2.62480    -3.16170    -0.98950
H    3.02960    -1.95740    -2.21330
H    4.28210    -2.58190    -1.14230
H    -0.38450    2.19040    1.49850

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-3.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd1.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd1_conf-3.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd1.xyz

0 1

