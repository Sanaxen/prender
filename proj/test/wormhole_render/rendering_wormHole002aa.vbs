Set objWshShell = WScript.CreateObject("WScript.Shell")

Set objUsrEnv = objWshShell.Environment("Process")

i = 1
x = -15*4

'While i <= 120
'	x = x + 0.5
'	i = i + 1
'Wend
x = x + 100*0.5
	
While i < 132
	objUsrEnv.Item("XX") = Cstr(x)
	objUsrEnv.Item("TT") = Cstr(0)
	objUsrEnv.Item("UU") = Cstr(0)
	If i >= 118 Then
		objUsrEnv.Item("AA") = Cstr(-1)
	Else
		objUsrEnv.Item("AA") = Cstr(1)
	End if

	c = objWshShell.run ("_test002aa.bat",6,1)
	c = objWshShell.Run("copyImage.bat wormHole002a*.bmp " & Cstr(i),6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.bmp",6,1)
	c = objWshShell.Run("cmd.exe /c del wormHole002a*.hdr",6,1)
	
	x = x + 0.5
	i = i + 1
Wend
x = x -1
