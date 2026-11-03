#!/usr/bin/env python3
"""Replay all order-ten blocks and hostile controls in scratch space."""
import hashlib,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='s146-carrier10-') as d:
 d=Path(d)
 for name in ['verify10.cpp','hostiles.py']+[f'block{i}.{ext}' for i in range(1,5) for ext in ('contexts','reps')]+['block5.reps']:shutil.copyfile(H/name,d/name)
 data=(H/'block5.contexts.part00').read_bytes()+(H/'block5.contexts.part01').read_bytes()
 if hashlib.sha256(data).hexdigest()!='3c20444ab7aa8325c90adf38d83430b6cdb5cd7ac11979337aba7253860977d8':raise RuntimeError('joined bank hash mismatch')
 (d/'block5.contexts').write_bytes(data)
 subprocess.run(['g++','-O3','-std=c++17','-o',str(d/'verify10'),str(d/'verify10.cpp')],check=True)
 for i in range(1,6):
  r=subprocess.run([str(d/'verify10'),str(i),str(d/f'block{i}.contexts'),str(d/f'block{i}.reps')],check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
 r=subprocess.run([sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(d/'hostiles.py')],check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
