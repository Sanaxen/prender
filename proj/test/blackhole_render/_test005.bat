set IMAGEDUMP=10

set IMAGE_X=4096*0.1
set IMAGE_Y=2048*0.1

set OMP_NUM_THREADS=10
:set SAMPLING=30
set SAMPLING=15
set SUPERSAMPLING=1
set NEXTEVENTESTIMATION=0

set XX=%1
set YY=%2
set ZZ=%3
set XXX=%4
set YYY=%5
set ZZZ=%6
set VX=%7
set VY=%8
set VZ=%9

shift
set ZDIR=%9

shift
set UPX=%9
shift
set UPY=%9
shift
set UPZ=%9

call render blackHole005.txt
