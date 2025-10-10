
for /l %%i in (0,1,72) do (
	set XX=%%i
	del /Q image\*.bmp
	call _test003.bat
	call copyImage.bat blackHole003*.bmp %%i
	del blackHole003*.bmp
)
