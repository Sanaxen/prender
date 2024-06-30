Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

step = 15
n = 360/step

for i = 1 to n
	objUsrEnv.Item("XX") = Cstr(i*step)

	c = objWshShell.run ("_test003.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole003*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole003*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole003*.hdr",6,1)
	
	objUsrEnv.Item("XX") = ""

Next
