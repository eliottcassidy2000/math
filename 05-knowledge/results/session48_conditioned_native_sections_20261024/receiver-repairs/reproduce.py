#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
P=Path(__file__).resolve().parent
for script in ('verify_manifest.py','verify_native.py','check.py','replay_filtration.py'):
 for optimized in (False,True):
  cmd=[sys.executable]+(['-O'] if optimized else [])+['-B',str(P/script)]
  subprocess.run(cmd,cwd=P,check=True)
