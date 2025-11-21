@echo off
REM Wrapper script to run P2RANK with PyCharm's bundled Java
REM Usage: run_p2rank.bat <p2rank_command>

set JAVA_HOME=C:\Program Files\JetBrains\PyCharm Community Edition 2024.2.1\jbr
set PATH=%JAVA_HOME%\bin;%PATH%

REM Change to P2RANK directory
cd /d "%~dp0p2rank_2.4.2"

REM Run P2RANK using the prank script
call prank.bat %*
