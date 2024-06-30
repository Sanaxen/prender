Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

i = 1
x = -15*4

While x < -1.3
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(0)

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 0.5
	i = i + 1
Wend
x = x -1

x = -1.3
for k = 1 to 20
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(10000000+x)

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 0.5
	i = i + 1
Next