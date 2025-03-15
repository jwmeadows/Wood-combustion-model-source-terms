@echo off
echo.
echo Copyright 1987-2023 ANSYS, Inc. All Rights Reserved.

set PYTHONHOME=%FLUENT_INC%\..\commonfiles\CPython\3_10\winx64\Release\python
set PYTHONPATH=%FLUENT_INC%\..\commonfiles\CPython\3_10\winx64\Release\python

:: check if VC++ exists
where cl > nul 2>&1
if not errorlevel 1 (goto nmake_c) else (goto clang_c)

:clang_c
echo Compiler and linker: Clang (builtin)
set FLUENT_UDF_COMPILER=clang
set FLUENT_UDF_CLANG=builtin
"%PYTHONPATH%\Scripts\scons.exe" -s
goto Exit

:nmake_c
echo Compiler and linker: Microsoft Visual C++
set FLUENT_UDF_COMPILER=msvc
nmake /S /C
goto Exit

:Exit: 
