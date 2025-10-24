
for /l %%i in (0,20,410) do (
	set XX=%%i
	del /Q image\*.bmp
	call _test004.bat
	call copyImage.bat blackHole004*.bmp %%i
	del blackHole004*.bmp
	del blackHole004*.hdr
)
for /l %%i in (411,1,500) do (
	set XX=%%i
	del /Q image\*.bmp
	call _test004.bat
	call copyImage.bat blackHole004*.bmp %%i
	del blackHole004*.bmp
	del blackHole004*.hdr
)
