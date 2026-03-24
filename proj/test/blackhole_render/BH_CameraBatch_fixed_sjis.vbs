' ====== Black Hole Camera Batch - plunge path (VBScript) ======
Option Explicit

'====================== setting ======================
Const MODE         = "plunge"
Const Frames       = 14000
Const ObitDivN     = 400
Const StepDivN     = 100.0
Const CamDirBlend  = 0.9   ' 0=BH方向, 1=接線方向, 0.7=接線70%+BH30%
Const dth_max_deg  = 0.5   ' 1フレームの最大角度変化（度）
Const PI           = 3.14159265358979

' Black Hole Center
Dim X_BH, Y_BH, Z_BH
X_BH = 165 : Y_BH = 45 : Z_BH = 50
Dim X0
X0 = 50
Dim r_p, r_h
r_h = 2.0
r_p = 3*Sqr(3)

Const Rotation     = 1
Const Rotation_end = -1

'===================================================
Function V(x,y,z) : Dim a: a = Array(x,y,z): V = a : End Function
Function VAdd(a,b) : VAdd = V(a(0)+b(0), a(1)+b(1), a(2)+b(2)) : End Function
Function VSub(a,b) : VSub = V(a(0)-b(0), a(1)-b(1), a(2)-b(2)) : End Function
Function VScale(a,s): VScale = V(a(0)*s, a(1)*s, a(2)*s) : End Function
Function VDot(a,b) : VDot = a(0)*b(0) + a(1)*b(1) + a(2)*b(2) : End Function
Function VCross(a,b)
  VCross = V( a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0) )
End Function
Function VLen(a) : VLen = Sqr( VDot(a,a) ) : End Function
Function VNorm(a)
  Dim n: n = VLen(a)
  If n = 0 Then VNorm = V(0,0,0) Else VNorm = VScale(a, 1/n)
End Function
Function VDist(a,b) : VDist = VLen(VSub(a,b)) : End Function

' ベクトルのSlerp（球面線形補間） t=0でa, t=1でb
Function VSlerp(a, b, t)
  Dim cosA, angle, sinA, ra, rb
  cosA = VDot(a, b)
  If cosA > 0.9999 Then
    ' ほぼ同じ方向なのでLinear補間
    VSlerp = VNorm(VAdd(VScale(a, 1-t), VScale(b, t)))
  ElseIf cosA < -0.9999 Then
    ' 逆方向の場合はそのままLinear
    VSlerp = VNorm(VAdd(VScale(a, 1-t), VScale(b, t)))
  Else
    angle = Atn(Sqr(1 - cosA*cosA) / cosA)  ' acos の代替
    sinA  = Sin(angle)
    ra    = Sin((1-t) * angle) / sinA
    rb    = Sin(t * angle) / sinA
    VSlerp = VNorm(VAdd(VScale(a, ra), VScale(b, rb)))
  End If
End Function

Function BuildCmd2(camPos, camDir, upVec, zdir, upvec3)
  Dim cmd
  cmd = "call _test005.bat" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(1)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(1)) & " " & CStr(camDir(2)) _
      & " " & CStr(upVec(0))  & " " & CStr(upVec(1))  & " " & CStr(upVec(2)) _
      & " " & CStr(zdir) _
      & " " & CStr(upvec3(0)) & " " & CStr(upvec3(1)) & " " & CStr(upvec3(2))
  BuildCmd2 = cmd
End Function

Const ForAppending = 8
Dim logFile, batFile
Sub WriteLog(msg)
  logFile.WriteLine Now & " - " & msg
End Sub

Function PlungeRadius(t_frac)
  Dim rr
  rr = r_p * Exp(Log(r_h / r_p) * t_frac)
  If rr < r_h * 1.001 Then rr = r_h * 1.001
  PlungeRadius = rr
End Function

'===================================================
Dim fso, C, i, k
Dim camPos, camDir, upVec, upvec3, Velocity
Dim xx, yy, zz
Dim th, th_accum, t_frac
Dim dist, dist1, rr
Dim rot, rot_start_idx
Dim count, stp, zdir
Dim n_sub, sub_i, dth_sub, dth_raw, dth_max_rad
Dim camDir_tangent, camDir_bh
Dim beta_vel2, approach_frac
Dim camDir_prev
Dim rotAxis, sinAngle, cosAngle, axv, adotv
Dim transFrames, transIdx, transT
Dim camDir_trans_start, upvec3_trans_start, camPos_trans_start
Dim camDir_trans_end, upvec3_trans_end, camPos_trans_end
Dim inTransition
Dim frameNo

Set fso = CreateObject("Scripting.FileSystemObject")
If fso.FileExists("run.log")    Then fso.DeleteFile "run.log",    True
If fso.FileExists("output.bat") Then fso.DeleteFile "output.bat", True
Set logFile = fso.OpenTextFile("run.log",    ForAppending, True)
Set batFile = fso.OpenTextFile("output.bat", ForAppending, True)

count       = 0
xx          = V(0,0,0)
yy          = V(0,0,0)
zz          = V(0,0,0)
th_accum    = 0.0
rot         = 0
zdir        = 1
dth_max_rad = dth_max_deg * PI / 180.0

C        = V(X_BH, Y_BH, Z_BH)
Velocity = V(0,0,0)
upVec    = V(0,1,0)
upvec3   = V(0,0,1)
camPos   = V(0,0,0)
camDir   = V(0,0,0)
camDir_prev = V(0,0,0)
stp      = 5
k        = 1.0
transFrames  = 10  ' 補完フレーム数
transIdx     = 0
inTransition = False
frameNo      = 0

For i = 0 To Frames Step stp

  '==================================================
  ' 接近フェーズ (rot=0): Z+方向オフセット、camDir=BH方向
  '==================================================
  If rot = 0 Then
    camPos(0) = X0 + C(0) - X0 * i * k / StepDivN
    camPos(1) = C(1)
    camPos(2) = C(2) + r_p
    camDir    = VNorm(VSub(C, camPos))
    ' upvec3: 前フレームのcamDirから現フレームへの回転をupvec3にも適用（ロドリゲス回転）
    If VLen(camDir_prev) > 0.001 Then
      rotAxis  = VCross(camDir_prev, camDir)
      sinAngle = VLen(rotAxis)
      cosAngle = VDot(camDir_prev, camDir)
      If sinAngle > 1.0e-8 Then
        rotAxis = VNorm(rotAxis)
        axv     = VCross(rotAxis, upvec3)
        adotv   = VDot(rotAxis, upvec3)
        upvec3  = VNorm(VAdd(VAdd( _
          VScale(upvec3, cosAngle), _
          VScale(axv, sinAngle)), _
          VScale(rotAxis, adotv * (1.0 - cosAngle))))
      End If
    End If
  End If

  dist  = VDist(camPos, C)
  dist1 = VDist(V(X0 + C(0) - X0*(i+stp)/StepDivN, C(1), C(2) + r_p), C)

  If dist < dist1 And VLen(xx) < 0.001 And rot = 0 Then rot = Rotation

  '==================================================
  ' 軌道フェーズ初期化（1回のみ）
  '==================================================
  If rot = Rotation And VLen(xx) < 0.001 Then
    xx     = V(1, 0, 0)
    yy     = V(0, 0, 1)
    zz     = V(0, 1, 0)
    upVec  = zz

    th_accum      = 0.0
    rot_start_idx = i
    count         = 1
    k             = 0.2

    ' 接近フェーズ最後のcamDir/upvec3/camPosを保存
    camDir_trans_start  = V(camDir(0), camDir(1), camDir(2))
    upvec3_trans_start  = V(upvec3(0), upvec3(1), upvec3(2))
    camPos_trans_start  = V(camPos(0), camPos(1), camPos(2))

    ' 軌道フェーズ最初のcamDir/upvec3/camPosを計算（補完終端）
    Dim th_first, rr_first
    rr_first = PlungeRadius(0.0)
    th_first = PI / 2
    camPos_trans_end = VAdd(C, VAdd(VScale(xx, rr_first * Cos(th_first)), VScale(yy, rr_first * Sin(th_first))))
    Dim cdt_end, cdb_end
    cdt_end = VNorm(VAdd(VScale(xx, -Sin(th_first)), VScale(yy, Cos(th_first))))
    cdb_end = VNorm(VSub(C, camPos_trans_end))
    camDir_trans_end  = VNorm(VAdd(VScale(cdt_end, CamDirBlend), VScale(cdb_end, 1 - CamDirBlend)))
    upvec3_trans_end  = V(upvec3(0), upvec3(1), upvec3(2))

    inTransition = True
    transIdx     = 0

    WriteLog "=== Orbit Start i=" & CStr(i) & " dist=" & CStr(dist)
  End If

  '==================================================
  ' 軌道フェーズ計算
  '==================================================
  If VLen(xx) > 0.001 Then

    If count <= ObitDivN Then
      t_frac = (count - 1) / CDbl(ObitDivN)
    Else
      t_frac = 1.0
    End If
    rr = PlungeRadius(t_frac)

    ' dthサブステップ
    dth_raw = (2 * PI / ObitDivN) * (r_p / rr) ^ 1.5
    If dth_raw > dth_max_rad Then
      n_sub = Int(dth_raw / dth_max_rad) + 1
    Else
      n_sub = 1
    End If
    dth_sub = dth_raw / n_sub
    For sub_i = 1 To n_sub
      th_accum = th_accum + dth_sub
    Next
    th = th_accum + PI / 2

    camPos = VAdd(C, VAdd(VScale(xx, rr * Cos(th)), VScale(yy, rr * Sin(th))))

    ' camDir = 接線方向*a + BH方向*(1-a) のブレンド
    camDir_tangent = VNorm(VAdd(VScale(xx, -Sin(th)), VScale(yy, Cos(th))))
    camDir_bh      = VNorm(VSub(C, camPos))
    camDir = VNorm(VAdd(VScale(camDir_tangent, CamDirBlend), VScale(camDir_bh, 1 - CamDirBlend)))
    upVec  = zz

    ' upvec3: 前フレームのcamDirから現フレームへの回転をupvec3にも適用（ロドリゲス回転）
    If VLen(camDir_prev) > 0.001 Then
      rotAxis  = VCross(camDir_prev, camDir)
      sinAngle = VLen(rotAxis)
      cosAngle = VDot(camDir_prev, camDir)
      If sinAngle > 1.0e-8 Then
        rotAxis = VNorm(rotAxis)
        axv     = VCross(rotAxis, upvec3)
        adotv   = VDot(rotAxis, upvec3)
        upvec3  = VNorm(VAdd(VAdd( _
          VScale(upvec3, cosAngle), _
          VScale(axv, sinAngle)), _
          VScale(rotAxis, adotv * (1.0 - cosAngle))))
      End If
    End If

    If count >= ObitDivN+1 Then
      WriteLog "=== Plunge complete i=" & CStr(i)
      dist = VDist(camPos, C)
      Velocity = VScale(camDir, Sqr(r_h / (2.0 * dist)))
      batFile.WriteLine "rem frame " & CStr(i) & " (plunge complete)"
      batFile.WriteLine BuildCmd2(camPos, camDir, Velocity, zdir, upvec3)
      batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & CStr(frameNo)
      batFile.WriteLine "del blackHole005*.bmp"
      batFile.WriteLine "del blackHole005*.hdr"
      batFile.WriteLine ""
      frameNo = frameNo + 1
      WScript.Quit
    End If

    count = count + 1
  End If

  '==================================================
  ' Velocity (beta)
  ' 接近フェーズ: 遠方では速度小さく、光子球付近で軌道速度に合わせる
  ' 軌道フェーズ: 円軌道速度 sqrt(r_h/(2*dist))
  '==================================================
  dist = VDist(camPos, C)
  If VLen(xx) > 0.001 Then
    ' 軌道フェーズ: 円軌道速度
    beta_vel2 = Sqr(r_h / (2.0 * dist))
  Else
    ' 接近フェーズ: 光子球までの距離に応じてゼロから軌道速度へ
    ' dist=遠方→betaapprox0, dist=r_p→beta=sqrt(r_h/(2*r_p))
    If dist > r_p Then
      approach_frac = r_p / dist   ' 0(遠方)→1(光子球)
    Else
      approach_frac = 1.0
    End If
    beta_vel2 = Sqr(r_h / (2.0 * dist)) * approach_frac
  End If
  If beta_vel2 > 0.9999 Then beta_vel2 = 0.9999

  ' 今フレームのcamDirを保存（次フレームのロドリゲス回転に使用）
  camDir_prev = V(camDir(0), camDir(1), camDir(2))

  Velocity = VScale(camDir, beta_vel2)

  WriteLog CStr(i) & " dist=" & CStr(dist) & " rot=" & CStr(rot)

  ' 補完フェーズ中はSlerp補完フレームを出力
  If inTransition Then
    Dim jj
    For jj = 0 To transFrames - 1
      transT = CDbl(jj) / CDbl(transFrames)
      Dim tDir, tUp, tPos, tVel
      tDir = VSlerp(camDir_trans_start, camDir_trans_end, transT)
      tUp  = VSlerp(upvec3_trans_start, upvec3_trans_end, transT)
      tPos = VAdd(VScale(camPos_trans_start, 1-transT), VScale(camPos_trans_end, transT))
      tVel = VScale(tDir, beta_vel2)
      batFile.WriteLine "rem frame_trans " & CStr(frameNo)
      batFile.WriteLine BuildCmd2(tPos, tDir, tVel, zdir, tUp)
      batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & CStr(frameNo)
      batFile.WriteLine "del blackHole005*.bmp"
      batFile.WriteLine "del blackHole005*.hdr"
      batFile.WriteLine ""
      frameNo = frameNo + 1
    Next
    inTransition = False
  End If

  batFile.WriteLine "rem frame " & CStr(i)
  batFile.WriteLine BuildCmd2(camPos, camDir, Velocity, zdir, upvec3)
  batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & CStr(frameNo)
  batFile.WriteLine "del blackHole005*.bmp"
  batFile.WriteLine "del blackHole005*.hdr"
  batFile.WriteLine ""
  frameNo = frameNo + 1

Next

logFile.Close
batFile.Close
WScript.Echo "Done."
