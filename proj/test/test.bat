set NEXTEVENTESTIMATION=1
set SAMPLING=400
set SUPERSAMPLING=2

set TIMELIMIT=100


set IMAGE_X=1024
set IMAGE_Y=768

set IMAGE_X=320
set IMAGE_Y=240


set IMAGE_X=1024
set IMAGE_Y=768

set IMAGE_X=640
set IMAGE_Y=480


:call render sponza_���s��������2.txt
call render ������2.txt

goto end

del *.ppm

set LDMPATH="C:\Users\vaio6\Desktop\�����_��\�����_��\�p�X�g���[�V���O\prender\x64\FullSpectralRelease"
set LDM=prender.exe

del image\*.* /Q
del *.ppm
set SAMPLING=60
set SUPERSAMPLING=2
%LDMPATH%\%LDM% full�X�y�N�g��Test.txt > log.txt

:end
