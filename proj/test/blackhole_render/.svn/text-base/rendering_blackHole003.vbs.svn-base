Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

step = 5
n = 360/step

for i = 0 to n
	objUsrEnv.Item("XX") = Cstr(i*step)

	c = objWshShell.run ("_test003.bat",6,1)
	c = objWshShell.Run("copyImage.bat blackHole003*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del blackHole003*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del blackHole003*.hdr",6,1)
	
	objUsrEnv.Item("XX") = ""

Next
