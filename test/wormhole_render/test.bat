set IMAGEDUMP=5

set IMAGE_X=4000*0.1
set IMAGE_Y=2000*0.1

set SAMPLING=1
set SUPERSAMPLING=1
set NEXTEVENTESTIMATION=0

:call render wormHole001.txt

set IMAGE_X=640
set IMAGE_Y=480
set SAMPLING=50
:call render wormHole.txt
:call render wormHoleé©å»2x.txt
:call render wormHoleé©å».txt


cscript wormHole003.vbs
cscript wormHole002a.vbs
cscript wormHole002aa.vbs
cscript wormHole000.vbs
