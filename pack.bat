@echo off
setlocal

:: 获取当前日期和时间
for /f "tokens=2 delims==." %%i in ('wmic os get localdatetime /value') do set datetime=%%i

:: 格式化日期和时间
set year=%datetime:~0,4%
set month=%datetime:~4,2%
set day=%datetime:~6,2%
set hour=%datetime:~8,2%
set minute=%datetime:~10,2%
set second=%datetime:~12,2%

set formatted_datetime=%year%%month%%day%_%hour%%minute%%second%

cp src/rerouting.cpp main.cpp
zip submit/rerouting_%formatted_datetime%.zip main.cpp
rm main.cpp