Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

n = 16

for i = 1 to n
	objUsrEnv.Item("XX") = Cstr(200+i)

	c = objWshShell.run ("_test000.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole000*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole000*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole000*.hdr",6,1)
	
	objUsrEnv.Item("XX") = ""

Next
