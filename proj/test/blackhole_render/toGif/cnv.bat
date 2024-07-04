:ffmpeg  -r 6 -y -i image%%d.bmp aaa.mp4
ffmpeg -r 6 -y -i image%%d.bmp -vcodec libx264 -pix_fmt yuv420p  aaa.mp4

:convert -geometry 45%% aaa.avi aaa.gif


