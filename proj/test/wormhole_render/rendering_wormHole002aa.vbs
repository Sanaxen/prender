Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

i = 1
x = -15*4

side_change = 26
x = x + 100*0.5
'While i <= 119
'	x = x + 0.5
'	i = i + 1
'Wend

xx = 0	
While i < 18
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	If i >= side_change Then
		objUsrEnv.Item("AA") = Cstr(1000000)
	Else
		objUsrEnv.Item("AA") = Cstr(0)
	End if

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 0.5
	i = i + 1
Wend

While i < 38
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	If i >= side_change Then
		objUsrEnv.Item("AA") = Cstr(1000000)
	Else
		objUsrEnv.Item("AA") = Cstr(0)
	End if

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 0.2
	i = i + 1
Wend

nn = 1
n = 10
dt = 180/n
While i < 60

	objUsrEnv.Item("XX") = Cstr(x)
	If nn > n Then 
		objUsrEnv.Item("TT") = Cstr(180)
	Else
		objUsrEnv.Item("TT") = Cstr(dt*nn)
	End if
	objUsrEnv.Item("UU") = Cstr(0)
	objUsrEnv.Item("AA") = Cstr(1000000)

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	If nn > n Then
		x = x - 0.2
	Else
		x = x + 0.0
	End if
	i = i + 1
	nn = nn + 1
Wend
