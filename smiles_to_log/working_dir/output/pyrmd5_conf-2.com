%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-2.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrmd5.xyz

0 1
O    -4.94930    -0.67770    0.08830
C    -3.99150    -1.44720    0.13610
C    -2.55830    -1.04530    0.08290
C    -1.57620    -2.04200    0.14860
N    -0.24740    -1.77270    0.10760
C    -2.08870    0.27880    -0.03030
C    -2.99870    1.49510    -0.11470
N    -0.76200    0.57110    -0.07320
C    0.11880    -0.46560    -0.00270
C    1.58990    -0.14550    -0.04990
C    2.04190    1.18100    -0.16270
C    3.40120    1.48100    -0.20660
C    4.33630    0.45360    -0.13810
C    3.91960    -0.86870    -0.02620
C    2.55870    -1.16220    0.01720
H    -4.15790    -2.52140    0.22480
H    -1.82640    -3.09060    0.23650
H    -3.62440    1.55880    0.77540
H    -2.44030    2.42810    -0.19790
H    -3.65020    1.41740    -0.98510
H    1.33510    1.99610    -0.21760
H    3.72610    2.50780    -0.29360
H    5.39160    0.68400    -0.17200
H    4.64730    -1.66550    0.02690
H    2.26080    -2.19660    0.10440

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrmd5.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrmd5_conf-2.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrmd5.xyz

0 1

