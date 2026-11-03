#!/usr/bin/env python3
"""Replay all order-nine blocks and hostile controls in scratch space."""
import shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='s146-carrier9-') as d:
 d=Path(d)
 for name in ['verify9.cpp','hostiles.py']+[f'block{i}.{ext}' for i in range(1,5) for ext in ('contexts','reps')]:shutil.copyfile(H/name,d/name)
 subprocess.run(['g++','-O3','-std=c++17','-o',str(d/'verify9'),str(d/'verify9.cpp')],check=True)
 for i in range(1,5):
  r=subprocess.run([str(d/'verify9'),str(i),str(d/f'block{i}.contexts'),str(d/f'block{i}.reps')],check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
 r=subprocess.run([sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(d/'hostiles.py')],check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
