set IMAGEDUMP=100

set IMAGE_X=640*3
set IMAGE_Y=480*3

set SAMPLING=10
set SUPERSAMPLING=1
set NEXTEVENTESTIMATION=0

goto 2

call render blackHole000.txt
call render blackHole00.txt
call render blackHole0.txt
call render blackHole01.txt
call render blackHole02.txt

:2
set IMAGE_X=640*2
set IMAGE_Y=480*2
call render blackHole001.txt
goto end


:set IMAGE_X=640
:set IMAGE_Y=480
:set SAMPLING=40
set SUPERSAMPLING=2
call render blackHole1.txt

:end