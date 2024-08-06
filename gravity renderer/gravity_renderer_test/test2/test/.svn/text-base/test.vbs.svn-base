Set objWshShell = WScript.CreateObject("WScript.Shell")

c = objWshShell.Run("cmd.exe /c del /Q *.bmp")
r = 10
step = 0.2
n = r/step

for i = 1 to n
	strCmd = "..\x64\Release\test.exe -back1 InterstellarWormhole_Fig10.jpg -t 10000 -thread 5 " & "-r " & Cstr(i*step)
	c = objWshShell.run (strCmd,1,1)
	c = objWshShell.Run("cmd.exe /c copy image.bmp toGif\image" & Cstr(n-i+1) & ".bmp")
Next
