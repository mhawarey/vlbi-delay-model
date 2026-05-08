@echo off
cd /d "%~dp0"

IF NOT EXIST ".git" (
    git init
    git branch -M main
)

git add .
git commit -m "Initial release: VLBI Delay Model v1.0 with screenshots"

git remote get-url origin >nul 2>&1
IF ERRORLEVEL 1 (
    gh repo create mhawarey/vlbi-delay-model --public --source=. --remote=origin --push
) ELSE (
    git push -u origin main
)

echo [DONE] https://github.com/mhawarey/vlbi-delay-model
pause
