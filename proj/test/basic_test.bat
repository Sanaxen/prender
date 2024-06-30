set OMP_NUM_THREADS=5

set IMAGE_X=640*0.25
set IMAGE_Y=480*0.25
set IMAGEDUMP=10
set NEXTEVENTESTIMATION=1

set SAMPLING=400

set SUPERSAMPLING=2
set TIMELIMIT=10000

set METROPOLIS=0
set MUTATION=30
set PRESAMPLE=0.1

:set ENVLIGHT_R=0.05
:set ENVLIGHT_G=0.025
:set ENVLIGHT_B=0.0005

set USE_BSSRDF=1
set SUFFIX_SYMBOL=_bssrdf
goto 1

call render CornellBox_����2.txt
call render CornellBox_����1.txt
goto end

goto 5
set IMAGE_X=804
set IMAGE_Y=1200
call render san-miguel3b.txt
goto end

:call render CornellBox���s����_�֗^�}������.txt
:call render CornellBox���s����.txt
:call render CornellBox_spot����.txt
:call render CornellBox_spot����_�֗^�}������.txt
:call render CornellBox�ʌ���.txt
:goto end

call render sponza_���s��������3_�{���̕��s����2.txt
goto end

:goto end


:1
call render CornellBox_SSSSWhite_Grapefruit_Juice_bunny_demo.txt
call render CornellBox_SSS�A�b�v��bunny_demo.txt
call render CornellBox_SSS�N���[��bunny_demo.txt
call render CornellBox_SSSRegular_Chocolate_Milk_bunny_demo.txt
call render CornellBox_SSSSkin2bunny_demo.txt
call render CornellBox_SSSSkin1bunny_demo.txt
call render CornellBox_SSS�|�e�gbunny_demo.txt
call render CornellBox_SSS�}�[�u��bunny_demo.txt
call render CornellBox_SSEspresso_bunny_demo.txt
call render CornellBox_SSSWholemilk_bunny_demo.txt
call render CornellBox_SSSSkimmilk_bunny_demo.txt
call render CornellBox_SSS�P�`���b�vbunny_demo.txt
:goto end

:goto end
call render CornellBox_SSSSpectralon_bunny_demo.txt
call render CornellBox_SSSorange_juice_bunny_demo.txt
call render CornellBox_SSS��2.txt
call render CornellBox_SSS����3.txt
call render CornellBox_SSS����3a.txt
call render CornellBox_SSS��.txt

goto end



:goto end
call render CornellBox_SSS�N���[��.txt
:goto end
call render CornellBox_SSS��.txt
:goto end
goto end

:goto end
goto end


set USE_BSSRDF=1
set SAMPLING=91
call render CornellBox_SSSSkin2bunny_demo.txt
call render CornellBox_SSS�}�[�u��bunny_demo.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS����3.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS�N���[��.txt
call render CornellBox_SSS����3a.txt
goto end

:call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt
:call render ����.txt
:call render CornellBox_����1.txt
:call render sponza_���s��������_�{���̕��s����2.txt
:call render California3.txt
:call render sponza_���s��������3_�{���̕��s����2.txt
:call render CornellBox_�K���X1.txt
:call render CornellBox_�K���X2.txt
:call render CornellBox_���ʔ���.txt
:goto end

set USE_BSSRDF=1
call render CornellBox_SSS����3a.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS�}�[�u��bunny_demo.txt
call render CornellBox_SSS�N���[��.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS����3.txt
call render CornellBox_SSS����4.txt

call render CornellBox_SSSbunny.txt
goto end

:call render bunny_test_sss.txt
goto end

:set BLACK_HOLE_Z=90
call render CornellBox_blackHole001.txt
goto end
call render CornellBox_blackHole00.txt
call render CornellBox_blackHole01.txt
call render CornellBox_blackHole02.txt
:call render CornellBox_spectrum2.txt
goto end

:call render ����.txt
:goto end
:goto 3

set USE_BSSRDF=1
call render CornellBox_SSS�N���[��.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS����3.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS����4.txt
call render CornellBox_SSS����3a.txt

:call render CornellBox_�g�U����.txt
goto end
set USE_BSSRDF=0
set SAMPLING=31
call render CornellBox_SSS����3a.txt
call render CornellBox_SSS����3.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS����4.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS�N���[��.txt
:call render CornellBox_�g�U����.txt

:goto end

:3
:goto end
:call render California3.txt
:call render CornellBox.txt
:call render CornellBox_SSS�N���[��.txt
call render CornellBox_spot����_�֗^�}������.txt
call render CornellBox_spot����.txt
call render CornellBox���s����_�֗^�}������.txt
call render CornellBox���s����.txt
:goto 4
goto end

set METROPOLIS=1
set MUTATION=40
set RESAMPLE=0.01
:call render CornellBox_spectrum.txt
:call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt
set SAMPLING=4000
set NEXTEVENTESTIMATION=1
set METROPOLIS=0
call render CornellBox_spectrum.txt
call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt
goto end

call render CornellBox_SSS��chk.txt
:call render texture_test5_bump.txt
:call render texture_test_bump.txt
goto end

:call render CornellBox3.txt
:call render CornellBox4.txt


goto aaa
set SAMPLING=200
set NEXTEVENTESTIMATION=0
call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt

set SAMPLING=100
set NEXTEVENTESTIMATION=1
call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt
:goto end

set SAMPLING=80
set ERPT=1
set METROPOLIS=0
set MUTATION=5
set MLTSAMPLE=20
set MLTRESAMPLE=0
call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt
:goto end

:aaa
set SAMPLING=1
set ERPT=0
set METROPOLIS=1
set MUTATION=50 -2
set MLTSAMPLE=5
set MLTRESAMPLE=0
call render CornellBox.txt
:call render CornellBoxMLT2_�Օ�����.txt

goto end


call render CornellBox_����.txt
:call render CornellBox_�v���Y��2.txt
:call render CornellBox_�v���Y��0.txt
:call render CornellBox_�v���Y��.txt
goto end




:goto end




:goto 2
:goto end

:4
:set SAMPLING=300
call render CornellBox.txt
call render CornellBox_����2.txt
call render CornellBox_����1.txt
call render CornellBox_�K���X1.txt
call render CornellBox_�K���X2.txt
call render CornellBox_���ʔ���.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS��.txt
call render CornellBox_SSS�N���[��.txt
call render CornellBox_�g�U����.txt
call render CornellBox_Phong2.txt
call render CornellBox_Phong1.txt
call render CornellBox_���ڌ����e�X�g.txt
call render CornellBox_���ڌ����e�X�g_yama.txt
call render CornellBox_participatingMedia�O���U��.txt
call render CornellBox_participatingMedia����U��.txt
call render CornellBox_participatingMedia�����U��.txt

call render CornellBox_spot����.txt
call render CornellBox_spot����_�֗^�}������.txt
call render CornellBox�ʌ���.txt
call render CornellBox���s����.txt
call render CornellBox���s����_�֗^�}������.txt

call render bunny_test.txt
call render bunny_test2.txt
call render bunny_test3.txt
call render texture_test_alp.txt
call render texture_test0.txt
call render texture_test1.txt
call render texture_test2.txt
call render texture_test3.txt
call render texture_test3b.txt
call render texture_test3c.txt
call render texture_test4.txt
call render texture_test5.txt
call render texture_test5_bump.txt
call render texture_test_bump.txt

goto end

:call render CornellBox.txt
:call render CornellBox_�v���Y��0.txt
:call render CornellBox_�v���Y��.txt
:call render CornellBox_participatingMedia2.txt
:call render CornellBox_participatingMedia.txt


:call render full�X�y�N�g��Test.txt
:goto end

:goto 2

:goto end


:goto 2
:goto end




:goto end


call render CornellBox_polygonTest2.txt
call render CornellBox.txt
call render CornellBox3.txt
call render CornellBox_SSS.txt
call render CornellBox_SSS2.txt
call render CornellBox_SSS3.txt
call render CornellBox_SSS6.txt
call render CornellBox_SSS5.txt
call render CornellBox_SSS4.txt

:goto end
:



set IMAGE_X=640
set IMAGE_Y=480
call render sponza_���s��������_�{���̕��s����.txt
:goto 3


call render sponza_���s��������3_�{���̕��s����.txt
:call render sponza_���s��������2.txt
:call render sponza_���s��������3.txt
goto end

:call render san-miguel2.txt
call render san-miguel5.txt
call render san-miguel6.txt
:goto end


set IMAGE_Y=596
set IMAGE_X=362
:call render san-miguel3.txt
call render san-miguel4.txt
goto end

:5
set IMAGE_X=640
set IMAGE_Y=480
call render ������3.txt
:call render ������2.txt
goto end






 


:3
set SAMPLING=50
set SUPERSAMPLING=2

set IMAGE_X=640
set IMAGE_Y=480

:set IMAGE_X=1024
:set IMAGE_Y=768
set TIMELIMIT=8



:call render ������2.txt
:call render ������.txt

call render sponza_���s��������.txt
call render sponza_���s��������2.txt
call render sponza_���s��������3.txt
goto end

:3
call render sponza_���s��������_���͌���Ȃ�.txt





set SAMPLING=400


set IMAGE_X=640
set IMAGE_Y=480
set SAMPLING=200

call render Venus.txt
call render sponza_���s����.txt
call render sponza_������.txt


:call render CornellBox_polygonTest1.txt > log.txt

:end
