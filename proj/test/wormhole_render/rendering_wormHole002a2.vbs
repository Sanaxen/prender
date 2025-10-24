Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

i = 1
x = -25

While x < -6
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(0)

	c = objWshShell.run ("_test002a2.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a2*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.hdr",6,1)
	
	x = x + 1
	i = i + 1
Wend
x = x -1


While x < -6
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(0)

	c = objWshShell.run ("_test002a2.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a2*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.hdr",6,1)
	
	x = x + 0.5
	i = i + 1
Wend
x = -6



j = i
step = 10
n = 360/step
a = 0
for k = 1 to n
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(k*step)
	objUsrEnv.Item("UU") = Cstr(k*step)
	objUsrEnv.Item("AA") = Cstr(0)

	c = objWshShell.run ("_test002a2.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a2*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.hdr",6,1)
	
	i = i + 1
Next

x = -6
step = 10
n = 180/step
for k = 1 to n
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(k*step)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(10000000+x)

	c = objWshShell.run ("_test002a2.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a2*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.hdr",6,1)
	
	x = x - 0.5
	i = i + 1
Next


While x < -6
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(180)
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(10000000+x)

	c = objWshShell.run ("_test002a2.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a2*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a2*.hdr",6,1)
	
	x = x + 1
	i = i + 1
Wend

x = -6
step = 10
n = 180/step
for k = 1 to n
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(180+k*step)
	objUsrEnv.Item("UU") = Cstr(0)
	'objUsrEnv.Item("AA") = Cstr(x)
	objUsrEnv.Item("AA") = Cstr(10000000+x)

	c = objWshShell.run ("_test002a.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x - 1
	i = i + 1
Next

for k = 1 to n
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	'objUsrEnv.Item("AA") = Cstr(x)
	objUsrEnv.Item("AA") = Cstr(10000000+x)

	c = objWshShell.run ("_test002a.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 1
	i = i + 1
Next

