%nprocs=16
%mem=16GB
%chk=pyrd12_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd12.xyz

0 1
C    -3.12500    -1.55370    0.25440
O    -1.71080    -1.39420    0.25040
C    -1.17240    -0.12760    0.13510
C    0.22480    -0.02270    0.13650
C    1.18720    -1.19850    0.27380
C    2.52830    -0.93140    -0.47520
C    3.09000    0.49470    -0.26420
N    2.14850    1.42850    -0.07630
O    4.27950    0.78000    -0.31000
C    0.78010    1.25380    0.02170
N    0.05940    2.39110    -0.08460
C    -1.28950    2.27440    -0.07590
C    -1.93960    1.04430    0.02690
H    -3.58470    -1.03480    1.09700
H    -3.57160    -1.19680    -0.67500
H    -3.36700    -2.61230    0.34940
H    1.38800    -1.36400    1.33390
H    0.72730    -2.11950    -0.08890
H    3.29850    -1.63820    -0.16820
H    2.39450    -1.04430    -1.55050
H    2.52050    2.35570    0.02530
H    -1.84960    3.19620    -0.15980
H    -3.01650    1.01940    0.02430

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd12_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd12.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd12_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd12.xyz

0 1

