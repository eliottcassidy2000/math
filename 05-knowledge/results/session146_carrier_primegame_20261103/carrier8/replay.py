#!/usr/bin/env python3
"""Replay order-eight quotient rank and hostile controls in scratch space."""
import hashlib,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='s146-carrier8-') as d:
 d=Path(d)
 for name in ('verify8.cpp','controls.py','reps6.txt','reps8.txt'):shutil.copyfile(H/name,d/name)
 data=(H/'contexts8.part00').read_bytes()+(H/'contexts8.part01').read_bytes()
 if hashlib.sha256(data).hexdigest()!='94df48ba7e434651698928b156d81ef253c39cf76f9df8c20a4c87b6adcad57b':raise RuntimeError('joined bank hash mismatch')
 (d/'contexts8.txt').write_bytes(data)
 subprocess.run(['g++','-O3','-std=c++17','-o',str(d/'verify'),str(d/'verify8.cpp')],check=True)
 for args in ([str(d/'verify'),'8',str(d/'contexts8.txt'),str(d/'reps8.txt')],[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(d/'controls.py')]):
  r=subprocess.run(args,check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
