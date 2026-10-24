"""Replay all finite censuses and native checks, with or without -O."""
import subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
flags=['-O']if sys.flags.optimize else []
for k in range(4,9):subprocess.run([sys.executable,*flags,str(HERE/'census.py'),str(k)],check=True)
for name in ('stabilization.py','check.py'):subprocess.run([sys.executable,*flags,str(HERE/name)],check=True)
