set LDMPATH="C:\Users\vaio6\Desktop\レンダラ\レンダラ\パストレーシング\prender\x64\Release"
:set LDMPATH="C:\Users\vaio6\Desktop\レンダラ\レンダラ\パストレーシング\prender\x64\FullSpectralRelease"

set LDM=prender.exe

del image\*.* /Q
del *.ppm

set time2=%time: =0%
set time3=%time2:~0,2%%time2:~3,2%%time2:~6,2%

:%LDMPATH%\%LDM% %1 > log_%~n1_%time3%.txt

%LDMPATH%\%LDM% %1 > log.txt

:pause