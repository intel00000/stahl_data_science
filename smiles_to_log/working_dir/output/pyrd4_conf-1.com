%nprocs=16
%mem=16GB
%chk=pyrd4_conf-1.chk
# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine) opt freq=noraman

 pyrd4.xyz

0 1
O    2.78370    0.94820    1.92220
C    1.93680    0.98400    2.81430
N    2.20360    1.15600    4.12250
C    3.41490    1.38740    4.78020
N    3.28300    1.44960    6.12860
C    4.39180    1.67300    6.87710
C    5.65840    1.84360    6.31820
C    5.79630    1.78340    4.92880
C    7.17100    1.96830    4.29530
C    7.61010    1.19940    3.06670
C    7.36990    2.69210    2.98000
C    4.65710    1.55270    4.14790
C    0.47420    0.78110    2.48130
C    -0.54690    1.35240    3.25990
C    -1.88380    1.14850    2.92700
C    -2.22590    0.37970    1.81350
B    -3.74700    0.15040    1.43900
O    -4.17680    -0.58820    0.35980
C    -5.57190    -0.77800    0.50800
C    -6.17450    -0.82060    -0.91080
C    -5.73350    -2.15140    1.20640
C    -5.97140    0.47480    1.37050
C    -6.16550    1.79530    0.58390
C    -7.14150    0.25410    2.35030
O    -4.80640    0.66490    2.15240
C    -1.21600    -0.17980    1.02950
C    0.12200    0.02160    1.35730
H    1.41480    1.08450    4.73720
H    4.24770    1.71350    7.94730
H    6.51470    2.01920    6.95190
H    7.95490    2.13740    5.03610
H    6.90350    0.51320    2.59620
H    8.63150    0.81820    3.03120
H    8.22160    3.36670    2.88300
H    6.49320    3.06490    2.44740
H    4.75540    1.50840    3.07600
H    -0.32450    1.97310    4.11640
H    -2.66280    1.59360    3.53020
H    -5.90520    0.05800    -1.49730
H    -7.26210    -0.87540    -0.87280
H    -5.81650    -1.68730    -1.46760
H    -6.78440    -2.41150    1.32990
H    -5.27140    -2.16070    2.19440
H    -5.26530    -2.94940    0.62880
H    -7.06940    1.76480    -0.02380
H    -5.32560    1.99770    -0.08210
H    -6.25210    2.65030    1.25540
H    -8.05920    0.01120    1.81530
H    -7.33100    1.14590    2.94870
H    -6.94130    -0.55490    3.05290
H    -1.47190    -0.77030    0.16090
H    0.88830    -0.41720    0.73110

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd4_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check

 pyrd4.xyz

0 1

--Link1--
%nprocs=16
%mem=16GB
%chk=pyrd4_conf-1.chk
# M062X/def2tzvp int=(grid=ultrafine) nmr=giao guess=read geom=check

 pyrd4.xyz

0 1

