ren %1 image%2.bmp
copy image%2.bmp toGif /v
del image%2.bmp 
:copy %1 toGif\image%2.bmp /v
:pause
