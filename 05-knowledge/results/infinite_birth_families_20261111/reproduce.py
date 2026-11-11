"""Replay producers and the incoming-engine comparison, then check payload hashes."""
from pathlib import Path
import subprocess,sys

HERE=Path(__file__).resolve().parent
for name in ['enumerate.py','audit.py','families.py','arithmetic.py','compatibility.py','compare_incoming.py']:
    subprocess.run([sys.executable,'-B',str(HERE/name)],check=True,cwd=HERE)
if (HERE/'SHA256SUMS').exists():
    subprocess.run([sys.executable,'-B',str(HERE/'seal.py'),'--check'],check=True,cwd=HERE)
print('PASS: all birth-family, arithmetic, source and child checks')
