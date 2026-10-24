"""Replay frozen native supply, controller, compatibility and fibre certificates."""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
flags=['-B']+(['-O'] if sys.flags.optimize else [])
for relative in ('supply/check.py','controller/check.py','controller/padding_probe.py','compatibility/check.py','fibres/reproduce.py'):
    subprocess.run([sys.executable,*flags,str(HERE/relative)],check=True)
