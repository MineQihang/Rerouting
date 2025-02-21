@echo off
setlocal enabledelayedexpansion

set totalScore=0
set count=0

for %%i in (01 02 03 04 05 06 07 08 09 10) do (
    echo Running test case %%i
    .\build\main.exe < .\testcases\%%i > .\output\%%i 2> .\output\%%i.score
    powershell -Command "Get-Content .\output\%%i.score | Select-Object -Last 2"
    for /f "tokens=3" %%j in ('findstr "Total score:" .\output\%%i.score') do (
        set /a totalScore+=%%j
        set /a count+=1
    )
)

if !count! gtr 0 (
    set /a averageScore=totalScore/count
    echo Average Total Score: !averageScore!
) else (
    echo No scores found.
)

endlocal