@echo off
cd /d "%~dp0"

IF NOT EXIST ".git" (
    git init
    git branch -M main
)

git add .
git commit -m "v1.1: add 2nd order ionospheric correction (Hawarey, Hobiger & Schuh 2005, GRL 32, L11304) + 5 validation tests"

git remote get-url origin >nul 2>&1
IF ERRORLEVEL 1 (
    gh repo create mhawarey/vlbi-delay-model --public --source=. --remote=origin --push
) ELSE (
    git push -u origin main
)

echo [DONE] https://github.com/mhawarey/vlbi-delay-model
pause
