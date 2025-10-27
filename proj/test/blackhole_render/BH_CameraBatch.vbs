' ====== Black Hole Camera Batch — flyby / plunge spiral path (VBScript) ======
Option Explicit

'====================== setting ======================
Const MODE         = "flyby"  ' "flyby" or "plunge"
Const Frames       = 1400      ' Total number of frames

' Black Hole Center (World Coordinate)
Dim X_BH, Y_BH, Z_BH
X_BH = 165 : Y_BH = 45 : Z_BH = 50

' Initial camera position (world coordinates)
Dim X0, Y0, Z0
X0 = 50 : Y0 = 45 : Z0 = 55 

' Initial Camera Direction (World)
Dim DX0, DY0, DZ0
DX0 = -1 : DY0 = 0 : DZ0 =0 

Dim r_p, r_h
r_p     = 3*Sqr(3)   ' Closest to the radius
r_h     = 2.0        ' Near the horizon


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
Function BuildCmd(camPos, camDir, V, zdir)
  Dim cmd
  
  cmd = """" & "_test005.bat" & """" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(1)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(1)) & " " & CStr(camDir(2)) _
      & " " & CStr(V(0)) & " " & CStr(V(1)) & " " & CStr(V(2)) & " " & CStr(zdir) & """"

  BuildCmd = cmd
End Function

Function BuildCmd2(camPos, camDir, V, zdir)
  Dim cmd
  
  cmd = "call _test005.bat" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(1)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(1)) & " " & CStr(camDir(2)) _
      & " " & CStr(V(0)) & " " & CStr(V(1)) & " " & CStr(V(2))  & " " & CStr(zdir)

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
Dim fso, sh, C, i, j, k, r
Dim camPos, camDir, cmd

Dim count
Dim xx,zz
Dim yy
Dim th,dist, dist1
Dim rot
Dim rot_start_idx
Dim rot_end_idx
Dim rr
Dim Velocity
Dim zdir
Dim stp

Const Rotation = 1       'Start Lap
Const Rotation_end = -1  'End Lap

Const ObitDivN = 250   'Number of divisions in circular motion
Const RotationN = 1    'Number of laps

Const StepDivN = 100.0 'Number of Approach Motion Segments

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

count = 0
xx = V(0,0,0)
yy = V(0,0,0)
zz = V(0,0,0)
th = 0.0
rot = 0
rot_end_idx = -1
zdir = 1

C = V(X_BH, Y_BH, Z_BH)
Velocity = V(0,0,0)

camPos = V(0,0,0)
camDir = V(0,0,0)
stp = 5
k = 1.0
For i = 0 To Frames Step stp

  camPos(0) = X0 +C(0) - X0*i*k/StepDivN
  camPos(1) = C(1)
  camPos(2) = C(2) + r_p
  

  camDir(0) = -1
  camDir(1) = 0
  camDir(2) = 0
  
  dist = VDist(camPos,C)
  dist1 = VDist(V(X0 +C(0) - X0*(i+stp)/StepDivN,C(1),C(2) + r_p),C)

  If dist < dist1 And VLen(xx) < 0.001 And rot = 0 Then rot = Rotation

  If rot = Rotation And VLen(xx) < 0.001 Then
    xx = camDir
    yy = VScale(VSub(camPos, C),-1)
    zz = VCross(xx, yy)
    yy = VCross(zz,xx)
    
    xx = VNorm(xx)
    yy = VNorm(yy)
    
    rot_start_idx = i
    count = 1
    k = 0.2
    Log "    xx: " & CStr(xx(0)) & "  " & CStr(xx(1)) & "  " & CStr(xx(2))
    Log "    yy: " & CStr(yy(0)) & "  " & CStr(yy(1)) & "  " & CStr(yy(2))
    Log "    zz: " & CStr(zz(0)) & "  " & CStr(zz(1)) & "  " & CStr(zz(2))
    
  End If
  

  'Log "    " & CStr(i) & " dist:" & CStr(VDist(camPos,C))

  th = RotationN*2*PI*count/ObitDivN
  If MODE="flyby" Then
     th = RotationN*PI*count/ObitDivN
  End If
  
  If MODE="plunge" Then
     th = RotationN*2*PI*count/ObitDivN
  End If
  

  If count > ObitDivN Then 
    th = 0
    rot = Rotation_end
    xx = V(0,0,0)
    
    If rot_end_idx = -1 Then rot_end_idx = i
    
    j =  rot_start_idx + (i-rot_end_idx)
    camPos(0) = X0 +C(0) - X0*j/StepDivN
    camPos(1) = C(1)
    camPos(2) = C(2) + r_p
  End If
  
  rr = r_h
  If VLen(xx) > 0.001 Then

    If MODE="flyby" Then
      rr = r_p
    End If

    If MODE="plunge" Then
      rr = r_p + (r_h - r_p)*count/ObitDivN
    End If

    'camPos(0) = VAdd(C, rr*xx(0)*Cos(th) + rr*yy(0)*Sin(th))
    'camPos(1) = VAdd(C, rr*xx(1)*Cos(th) + rr*yy(1)*Sin(th))
    'camPos(2) = VAdd(C, rr*xx(2)*Cos(th) + rr*yy(2)*Sin(th))
    

    camPos = VAdd(C, VAdd(VScale(xx, rr*Cos(th)), VScale(yy, rr*Sin(th))))
    
    'camDir(0) =  -rr*xx(0)*Sin(th) + rr*yy(0)*Cos(th)
    'camDir(1) =  -rr*xx(1)*Sin(th) + rr*yy(1)*Cos(th)
    'camDir(2) =  -rr*xx(2)*Sin(th) + rr*yy(2)*Cos(th)
    
     camDir = VSub(VScale(xx, -rr*Sin(th)), VScale(yy, rr*Cos(th)))
     camDir = VNorm(camDir)
     
     Log "    camPos: " & CStr(camPos(0)) & "  " & CStr(camPos(1)) & "  " & CStr(camPos(2))
     Log "    camDir: " & CStr(camDir(0)) & "  " & CStr(camDir(1)) & "  " & CStr(camDir(2))
  End If
  
  If ObitDivN/2 < count And VLen(xx) > 0.001 Then
    'zdir = -1
  End If
  
  dist = VDist(camPos,C)
  Velocity = VScale(camDir, Sqr(r_h/dist))
  If dist > r_p Then
     Velocity = VScale(Velocity, 1/((dist - r_p) + 1)^2)
  End If
  
  Log "    " & CStr(i) & " dist:" & CStr(VDist(camPos,C)) & " rot:" & CStr(rot) & " count:" & CStr(count) & " Velocity:" & CStr(VLen(Velocity))
  count = count + 1
  
  ' コマンド生成 & 実行（待機=True）
  
  cmd = BuildCmd(camPos, camDir, Velocity, zdir)
  'WScript.Echo cmd ' ←Debug
  
  batFile.WriteLine "rem " & Cstr(i)
  batFile.WriteLine BuildCmd2(camPos, camDir, Velocity, zdir)
  batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & Cstr(i)
  batFile.WriteLine "del blackHole005*.bmp"
  batFile.WriteLine "del blackHole005*.hdr"
  batFile.WriteLine ""

  Log CStr(i) & " " & cmd

  'If i >= 725 Then
  '   sh.Run cmd, 6, True
  '   sh.Run "call copyImage.bat blackHole005*.bmp " & Cstr(i),6,True
  '   sh.Run "cmd.exe /c del blackHole005*.bmp",6,True
  '   sh.Run "cmd.exe /c del blackHole005*.hdr",6,True
  'End If
  
  
  Log "---------------------------------------------------------------"
  
  'If camPos(0) < 140 Then WScript.Quit
Next

logFile.Close
batFile.Close
WScript.Echo "Done. Frames "
' ============================================================================

