Option Explicit



'===================================================

' Vector Utility (array [0..2])
Function V(x,y,z) : Dim a: a = Array(x,y,z): V = a : End Function
Function VAdd(a,b) : VAdd = V(a(0)+b(0), a(1)+b(1), a(2)+b(2)) : End Function
Function VSub(a,b) : VSub = V(a(0)-b(0), a(1)-b(1), a(2)-b(2)) : End Function
Function VScale(a,s): VScale = V(a(0)*s, a(1)*s, a(2)*s) : End Function
Function VDot(a,b) : VDot = a(0)*b(0) + a(1)*b(1) + a(2)*b(2) : End Function
Function VCross(a,b)
  VCross = V( a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0) )
End Function
Function VLen(a) : VLen = Sqr( VDot(a,a) ) : End Function
Function VNorm(a) : Dim n: n = VLen(a): If n=0 Then VNorm=V(0,0,0) Else VNorm=VScale(a,1/n) End If : End Function
Function VDist(a,b) : VDist = Sqr( VDot(VSub(a,b),VSub(a,b)) ) : End Function

'Trigonometric functions
Const PI = 3.1415926535897931
Function Atan2(y,x)
  If x > 0 Then
    Atan2 = Atn(y/x)
  ElseIf x < 0 And y >= 0 Then
    Atan2 = Atn(y/x) + PI
  ElseIf x < 0 And y < 0 Then
    Atan2 = Atn(y/x) - PI
  ElseIf x = 0 And y > 0 Then
    Atan2 = PI/2
  ElseIf x = 0 And y < 0 Then
    Atan2 = -PI/2
  Else
    Atan2 = 0
  End If
End Function

Function renorm(x)
  x = CInt(x*10000)
  renorm = x/10000.0
End Function

' Command Line Generation
Function BuildCmd(camPos, camDir)
  Dim cmd
  
  cmd = """" & "_test001a.bat" & """" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(2)) & """"

  BuildCmd = cmd
End Function

Function BuildCmd2(camPos, camDir)
  Dim cmd
  
  cmd = "call _test001a.bat" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(2)) 

  BuildCmd2 = cmd
End Function

Const ForAppending = 8
Dim logFile, batFile

Sub Log(msg)
  Dim line
  line = Now & " - " & msg
  logFile.WriteLine line
End Sub


' ==================== Main ====================
Dim fso, sh, i
Dim cmd


Set fso = CreateObject("Scripting.FileSystemObject")
Set sh  = CreateObject("WScript.Shell")

If fso.FileExists("run.log") Then
    fso.DeleteFile "run.log", True
End If
If fso.FileExists("output.bat") Then
    fso.DeleteFile "output.bat", True
End If

Set logFile = fso.OpenTextFile("run.log", ForAppending, True) 
Set batFile = fso.OpenTextFile("output.bat", ForAppending, True) 


Dim x, z
Dim xx, zz
Dim a,b,c,d
Dim camPos, camDir

a = (220 -180)/5
b = (5 - 70)/5
c = (-1 + 0.2)/5
d = (0 + 0.95)/5

For i = 0 To 5 Step 1

  x = 220 - a*i
  z = 5 - b*i
  

  xx = -1 - c*i
  zz = 0 - d*i
  

  ' コマンド生成 & 実行（待機=True）
  camPos = V(x,0, z)
  camdir = V(xx,0, zz)
  
  cmd = BuildCmd(camPos, camDir)
  'WScript.Echo cmd ' ←Debug
  
  batFile.WriteLine "rem " & Cstr(i)
  batFile.WriteLine BuildCmd2(camPos, camDir)
  batFile.WriteLine "call copyImage.bat blackHole001a*.bmp " & Cstr(i)
  batFile.WriteLine "del blackHole001a*.bmp"
  batFile.WriteLine "del blackHole001a*.hdr"
  batFile.WriteLine ""

  Log CStr(i) & " " & cmd

  'If i >= 725 Then
  '   sh.Run cmd, 6, True
  '   sh.Run "call copyImage.bat blackHole001a*.bmp " & Cstr(i),6,True
  '   sh.Run "cmd.exe /c del blackHole001a*.bmp",6,True
  '   sh.Run "cmd.exe /c del blackHole001a*.hdr",6,True
  'End If
  
  
  Log "---------------------------------------------------------------"
  
  'If camPos(0) < 140 Then WScript.Quit
Next

logFile.Close
batFile.Close
WScript.Echo "Done. Frames "
' ============================================================================

