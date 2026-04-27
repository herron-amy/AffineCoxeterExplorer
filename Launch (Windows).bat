@echo off
cd /d "%~dp0"

WHERE java >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    echo Java not found. Please install Java 11+ from https://adoptium.net
    pause & exit /b 1
)

echo Starting Affine Coxeter Explorer...
java --source 11 AffineCoxeterExplorer.java
IF %ERRORLEVEL% NEQ 0 ( echo Error occurred. & pause )
