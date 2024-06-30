set IMAGE_X=640
set IMAGE_Y=480


set IMAGEDUMP=5

set SUPERSAMPLING=2
set SAMPLING=10
set NEXTEVENTESTIMATION=1

set TIMELIMIT=10000

set METROPOLIS=0
set MUTATION=40
set RESAMPLE=0.01


call render CornellBox.txt
call render CornellBox_金属2.txt
call render CornellBox_金属1.txt
call render CornellBox_ガラス1.txt
call render CornellBox_ガラス2.txt
call render CornellBox_鏡面反射.txt
call render CornellBox_SSS黄.txt
call render CornellBox_SSS緑.txt
call render CornellBox_SSSクリーム.txt
call render CornellBox_SSS透明3.txt
call render CornellBox_SSS透明3a.txt
call render CornellBox_SSS透明4.txt
call render CornellBox_拡散反射.txt
call render CornellBox_Phong2.txt
call render CornellBox_Phong1.txt
call render CornellBox_直接光源テスト.txt
call render CornellBox_直接光源テスト_yama.txt
call render CornellBox_participatingMedia前方散乱.txt
call render CornellBox_participatingMedia後方散乱.txt
call render CornellBox_participatingMedia等方散乱.txt

call render CornellBox_spot光源.txt
call render CornellBox_spot光源_関与媒質あり.txt
call render CornellBox面光源.txt
call render CornellBox平行光源.txt
call render CornellBox平行光源_関与媒質あり.txt


call render CornellBoxMLT2_遮蔽あり.txt

call render CornellBox_SSSマーブルbunny_demo.txt
call render CornellBox_SSSbunny.txt

call render Venus.txt

call render sponza_平行光源実験.txt
call render sponza_平行光源実験2.txt
call render sponza_平行光源実験3.txt
call render sponza_平行光源実験_本当の平行光源.txt
call render sponza_平行光源実験3_本当の平行光源.txt
call render sponza_平行光源実験3_本当の平行光源2.txt


call render san-miguel2.txt
call render san-miguel5.txt
call render san-miguel6.txt


set IMAGE_Y=596
set IMAGE_X=362
call render san-miguel3.txt
call render san-miguel4.txt
goto end


set IMAGE_X=640
set IMAGE_Y=480
call render 博物館2.txt
call render 博物館.txt


call render CornellBox_spectrum.txt
:end
