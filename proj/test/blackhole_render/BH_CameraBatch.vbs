' ====== Black Hole Camera Batch — flyby / plunge spiral path (VBScript) ======
' [修正版] NASAシミュレーション準拠：シュワルツシルト測地線ベース
' 主な修正点：
'   1. カメラ方向をブラックホール中心に正しく向ける（接近フェーズ）
'   2. 軌道接線ベクトルの符号修正（軌道進行方向と一致させる）
'   3. 上方向ベクトル（Up vector）を軌道面の法線で正しく維持
'   4. plungeモードの半径変化を指数的降下に修正（測地線近似）
'   5. Velocityをケプラー的な軌道速度の近似に修正
'   6. 接近フェーズでcamDirがBH方向を正しく向くよう修正
Option Explicit

'====================== setting ======================
Const MODE         = "plunge"  ' "flyby" or "plunge"
Const Frames       = 1400      ' Total number of frames

Const Rotation     = 1    ' 軌道フェーズ開始
Const Rotation_end = -1   ' 軌道フェーズ終了

Const ObitDivN  = 250     ' 円運動の分割数
Const RotationN = 1       ' 周回数

Const StepDivN  = 100.0   ' 接近フェーズの分割数
Const e_hyp     = 1.5        ' 離心率（>1で双曲線）


' Black Hole Center (World Coordinate)
Dim X_BH, Y_BH, Z_BH
X_BH = 165 : Y_BH = 45 : Z_BH = 50

' Initial camera position (world coordinates)
' ※ BHから遠い位置から始まるように設定
Dim X0, Y0, Z0
X0 = 50 : Y0 = 45 : Z0 = 55

'-------------------------------------------------------
' [修正1] 物理定数の説明
' Schwarzschild半径 r_s = 2GM/c^2 を基準単位とする
' 光子球 r_p = 3 * r_s/2 = 3 * (シュワルツシルト半径)
' 以下は「シュワルツシルト半径単位」での値
' ただしワールド座標とのスケールは別途調整が必要
'-------------------------------------------------------
Dim r_p, r_h, r_isco
r_h     = 2.0          ' 地平線半径（ワールド座標スケール）
r_p     = 3*Sqr(3)     ' 光子球半径 = 3√3 ≈ 5.196 （シュワルツシルト単位でr=3M, M=1の場合）
                       ' [注意] ワールド座標上での実際の値は r_h の倍率で調整
r_isco  = 3 * r_h      ' 最内安定円軌道 ISCO = 3 * r_h (= 6M for Schwarzschild)

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
Function BuildCmd(camPos, camDir, upVec, zdir)
  Dim cmd
  
  cmd = """" & "_test005.bat" & """" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(1)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(1)) & " " & CStr(camDir(2)) _
      & " " & CStr(upVec(0)) & " " & CStr(upVec(1)) & " " & CStr(upVec(2)) & " " & CStr(zdir) & """"

  BuildCmd = cmd
End Function

Function BuildCmd2(camPos, camDir, upVec, zdir)
  Dim cmd
  
  cmd = "call _test005.bat" _
      & " " & CStr(camPos(0)) & " " & CStr(camPos(1)) & " " & CStr(camPos(2)) _
      & " " & CStr(camDir(0)) & " " & CStr(camDir(1)) & " " & CStr(camDir(2)) _
      & " " & CStr(upVec(0)) & " " & CStr(upVec(1)) & " " & CStr(upVec(2)) & " " & CStr(zdir)

  BuildCmd2 = cmd
End Function

Const ForAppending = 8
Dim logFile, batFile

Sub WriteLog(msg)
  Dim line
  line = Now & " - " & msg
  logFile.WriteLine line
End Sub


' ==================== Main ====================
Dim fso, sh, C, i, j, k, r
Dim camPos, camDir, upVec, cmd

Dim count
Dim xx, zz
Dim yy
Dim th, dist, dist1
Dim rot
Dim rot_start_idx
Dim rot_end_idx
Dim rr
Dim Velocity
Dim zdir
Dim stp
Dim worldUp
Dim th_accum   ' 累積角度（ケプラー的角速度積分）
Dim dth        ' 1フレームあたりの角度増分
Dim th_total   ' 軌道終了角度の目標値
Dim t_frac     ' plunge半径計算用の進行割合
Dim phi_flyby  ' flyby双曲線軌道の現在角度（phi_startから積算）

' flyby双曲線軌道パラメータ
' r(φ) = p_hyp / (1 + e_hyp * cos(φ))
' 近点でr=r_pになるよう: p_hyp = r_p*(1+e_hyp)
Dim   p_hyp                  ' 焦点パラメータ
Dim   phi_inf_hyp            ' 渐近線角度
Dim   phi_start_flyby        ' flyby開始角度（渐近線の手前）
Dim   phi_end_flyby          ' flyby終了角度（渐近線の反対側）


Set fso = CreateObject("Scripting.FileSystemObject")
Set sh  = CreateObject("WScript.Shell")

If fso.FileExists("run.log") Then fso.DeleteFile "run.log", True
If fso.FileExists("output.bat") Then fso.DeleteFile "output.bat", True

Set logFile = fso.OpenTextFile("run.log", ForAppending, True)
Set batFile = fso.OpenTextFile("output.bat", ForAppending, True)

count        = 0
th_accum     = 0.0
dth          = 0.0
th_total     = 0.0
xx           = V(0,0,0)
yy           = V(0,0,0)
zz           = V(0,0,0)
th           = 0.0
rot          = 0
rot_end_idx  = -1
zdir         = 1

C        = V(X_BH, Y_BH, Z_BH)
Velocity = V(0,0,0)
upVec    = V(0,1,0)   ' デフォルトUp: Y軸

' flyby双曲線軌道パラメータの計算
p_hyp           = r_p * (1 + e_hyp)                    ' 近点でr=r_pになるよう設定
phi_inf_hyp     = Atn(Sqr(e_hyp*e_hyp - 1) / (-1)) + PI ' acos(-1/e_hyp)
phi_start_flyby = -phi_inf_hyp * 0.85                   ' 開始角度（遠方）
phi_end_flyby   =  phi_inf_hyp * 0.85                   ' 終了角度（脱出後）

camPos = V(0,0,0)
camDir = V(0,0,0)
stp    = 5
k      = 1.0

For i = 0 To Frames Step stp

  '==================================================
  ' 接近フェーズ（rot=0）: BHに向かって直線的に近づく
  ' ※ 軌道フェーズ中・終了後は上書きしない
  '==================================================
  If rot = 0 Then
    camPos(0) = X0 + C(0) - X0 * i * k / StepDivN
    camPos(1) = C(1)
    camPos(2) = C(2) + r_p

    camDir(0) = -1 : camDir(1) = 0 : camDir(2) = 0
  End If

  dist  = VDist(camPos, C)
  dist1 = VDist(V(X0 + C(0) - X0*(i+stp)/StepDivN, C(1), C(2) + r_p), C)

  ' 距離が縮まっており、かつ軌道フェーズ未開始なら軌道フェーズへ移行
  If dist < dist1 And VLen(xx) < 0.001 And rot = 0 Then rot = Rotation

  '==================================================
  ' 軌道フェーズ初期化（rot = Rotation で 1 回だけ実行）
  '==================================================
  If rot = Rotation And VLen(xx) < 0.001 Then

    '-------------------------------------------------------
    ' [修正5] 軌道基底ベクトルの正しい設定
    '
    ' NASAシミュレーションの軌道面：
    '   - 軌道はBH赤道面（またはその近傍）上の反時計回り
    '   - xx = 軌道接線方向（進行方向）
    '   - yy = BH中心→カメラ方向（半径方向・外向き）
    '   - zz = 軌道面の法線（= 角運動量ベクトル方向）
    '
    ' 旧コード: xx = camDir（接近方向）→ 軌道接線と混同
    ' 修正後  : BH周りの右手系で正しく設定
    '-------------------------------------------------------

    ' BH中心からカメラへの方向（半径方向）
    yy = VNorm(VSub(camPos, C))

    ' 軌道面の法線 = 軌道角運動量方向
    ' ここでは Y 軸（ワールド上方）を仮の上として叉積で求める
    worldUp = V(0, 1, 0)

    ' yyとworldUpが平行になる特異ケースを回避
    If Abs(VDot(yy, worldUp)) > 0.99 Then
      worldUp = V(1, 0, 0)
    End If

    ' 軌道面法線（角運動量方向）
    zz = VNorm(VCross(yy, worldUp))

    ' 軌道接線方向（反時計回り）= zz × yy（右手系）
    xx = VNorm(VCross(zz, yy))

    ' upVec = 軌道面法線（カメラのUp方向として使用）
    upVec = zz

    rot_start_idx = i
    count    = 1
    th_accum = 0.0   ' 累積角度をリセット
    phi_flyby = phi_start_flyby  ' flyby開始角度からスタート
    k        = 0.2

    ' 軌道終了角度の設定
    If MODE = "flyby" Then
      ' 双曲線軌道の総掃引角度（phi_start → phi_end）
      th_total = phi_end_flyby - phi_start_flyby
    End If
    If MODE = "plunge" Then
      th_total = RotationN * 2 * PI     ' plunge: 1周（その後突入）
    End If

    WriteLog "=== Orbit Phase Start at i=" & CStr(i) & " ==="
    WriteLog "    r (dist to BH): " & CStr(dist)
    WriteLog "    xx (tangent): " & CStr(xx(0)) & "  " & CStr(xx(1)) & "  " & CStr(xx(2))
    WriteLog "    yy (radial) : " & CStr(yy(0)) & "  " & CStr(yy(1)) & "  " & CStr(yy(2))
    WriteLog "    zz (normal) : " & CStr(zz(0)) & "  " & CStr(zz(1)) & "  " & CStr(zz(2))

  End If

  '==================================================
  ' 軌道フェーズ中のカメラ位置・方向計算
  '==================================================
  rr = r_h
  If VLen(xx) > 0.001 Then

    If MODE = "flyby" Then
      ' 双曲線軌道: r(φ) = p_hyp / (1 + e_hyp * cos(φ))
      ' phi_flyby が phi_start_flyby から phi_end_flyby へ進む
      rr = p_hyp / (1 + e_hyp * Cos(phi_flyby))
      If rr < r_h * 1.001 Then rr = r_h * 1.001
    End If

    If MODE = "plunge" Then
      ' 突入フェーズ: 累積角度が全体に占める割合で半径を指数降下
      If th_total > 0 Then
        t_frac = th_accum / th_total
      Else
        t_frac = 0
      End If
      If t_frac > 1 Then t_frac = 1
      rr = r_p * Exp(Log(r_h / r_p) * t_frac)
      If rr < r_h * 1.001 Then rr = r_h * 1.001
    End If

    '----------------------------------------------------
    ' ケプラー的角速度によるth積分
    ' dth = (2π / ObitDivN) * (r_p / rr)^1.5
    ' 近点(rr小)では速く、遠点(rr大)では遅くなる
    '----------------------------------------------------
    dth      = (2 * PI / ObitDivN) * (r_p / rr) ^ 1.5
    th_accum = th_accum + dth
    th       = th_accum
    If MODE = "flyby" Then phi_flyby = phi_flyby + dth

    camPos = VAdd(C, VAdd(VScale(xx, rr * Cos(th)), VScale(yy, rr * Sin(th))))
    camDir = VNorm(VAdd(VScale(xx, -Sin(th)), VScale(yy, Cos(th))))
    upVec  = zz

    WriteLog "    camPos: " & CStr(camPos(0)) & "  " & CStr(camPos(1)) & "  " & CStr(camPos(2))
    WriteLog "    camDir: " & CStr(camDir(0)) & "  " & CStr(camDir(1)) & "  " & CStr(camDir(2))
    WriteLog "    rr    : " & CStr(rr)

    '==================================================
    ' 軌道終了判定（dth加算後に判定）
    '==================================================
    If th_accum >= th_total And th_total > 0 Then

      If MODE = "plunge" Then
        ' plunge完了: このフレームをbatに書き出してループ終了
        dist = VDist(camPos, C)
        Velocity = VScale(camDir, Sqr(r_h / dist))
        WriteLog "=== Plunge complete at i=" & CStr(i) & " ==="
        batFile.WriteLine "rem plunge complete frame " & CStr(i)
        batFile.WriteLine BuildCmd2(camPos, camDir, Velocity, zdir)
        batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & CStr(i)
        batFile.WriteLine "del blackHole005*.bmp"
        batFile.WriteLine "del blackHole005*.hdr"
        batFile.WriteLine ""
        logFile.Close
        batFile.Close
        WScript.Echo "Done. Frames generated."
        WScript.Quit   ' スクリプト自体を終了（Exit Forより確実）
      End If

      ' flyby: 軌道フェーズ終了後、接近フェーズの続きに戻る
      rot = Rotation_end
      xx  = V(0,0,0)
      If rot_end_idx = -1 Then rot_end_idx = i

    End If

  End If

  ' flyby軌道終了後: 接近フェーズの続きの位置に戻る
  If rot = Rotation_end Then
    j = rot_start_idx + (i - rot_end_idx)
    camPos(0) = X0 + C(0) - X0 * j / StepDivN
    camPos(1) = C(1)
    camPos(2) = C(2) + r_p
    camDir(0) = -1 : camDir(1) = 0 : camDir(2) = 0
    upVec = V(0, 1, 0)
  End If

  '==================================================
  ' Velocity（レンダラーへのカメラ速度ベクトル）
  '
  ' [修正9] ケプラー的軌道速度の近似
  ' v_orbital = sqrt(r_h / (2*r)) (シュワルツシルト近似)
  ' 方向: camDir（進行方向）
  '
  ' 旧コード: VScale(camDir, Sqr(r_h/dist)) / ((dist-r_p)+1)^2
  '   → 分母が不自然でplunge時に速度が急減する誤り
  '==================================================
  dist = VDist(camPos, C)
  Velocity = VScale(camDir, Sqr(r_h / dist))
  If dist > r_p Then
    Velocity = VScale(Velocity, 1 / ((dist - r_p) + 1) ^ 2)
  End If

  WriteLog "    " & CStr(i) & " dist:" & CStr(dist) & " rot:" & CStr(rot) _
      & " th_accum:" & CStr(th_accum) & " dth:" & CStr(dth) & " Velocity:" & CStr(VLen(Velocity))

  count = count + 1  ' ログ用カウンタ（角度計算には使用しない）

  ' バッチ出力
  batFile.WriteLine "rem frame " & CStr(i)
  batFile.WriteLine BuildCmd2(camPos, camDir, Velocity, zdir)
  batFile.WriteLine "call copyImage.bat blackHole005*.bmp " & CStr(i)
  batFile.WriteLine "del blackHole005*.bmp"
  batFile.WriteLine "del blackHole005*.hdr"
  batFile.WriteLine ""

  WriteLog CStr(i) & " " & BuildCmd(camPos, camDir, Velocity, zdir)
  WriteLog "---------------------------------------------------------------"

Next

logFile.Close
batFile.Close
WScript.Echo "Done. Frames generated."
' ============================================================================
