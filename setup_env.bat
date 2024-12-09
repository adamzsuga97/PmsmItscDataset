@echo off
REM Check if requirements.txt exists
IF NOT EXIST requirements.txt (
    echo requirements.txt not found.
    exit /b 1
)

REM Create virtual environment named venv
python -m venv venv

REM Activate the virtual environment
call venv\Scripts\activate

REM Upgrade pip
python -m pip install --upgrade pip

REM Install the dependencies from requirements.txt
pip install -r requirements.txt

REM Deactivate the virtual environment
activate

REM Download FEMM installer
echo Downloading FEMM installer...
curl -L -o femm42.exe http://www.femm.info/wiki/uploads/Archives/femm42bin_x64.exe

REM Install FEMM
echo Installing FEMM...
start /wait femm42.exe /SILENT

REM Clean up installer
del femm42.exe

echo Virtual environment setup complete and FEMM installed.