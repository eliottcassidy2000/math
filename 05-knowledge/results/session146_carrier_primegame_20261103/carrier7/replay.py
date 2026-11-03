#!/usr/bin/env python3
"""Replay certificate in scratch space; do not modify immutable source."""
import hashlib,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='s146-carrier-') as d:
 d=Path(d)
 for name in ('verify.cpp','controls.py','reps6.txt','reps7.txt'):
  shutil.copyfile(H/name,d/name)
 data=(H/'contexts7.part00').read_bytes()+(H/'contexts7.part01').read_bytes()
 if hashlib.sha256(data).hexdigest()!='ee0c332d4dd580eb3073cf95eba3ca098d005cd1f3f26e14c0aaaa9c6671baef':raise RuntimeError('joined bank hash mismatch')
 (d/'contexts7.txt').write_bytes(data)
 subprocess.run(['g++','-O3','-std=c++17','-o',str(d/'verify'),str(d/'verify.cpp')],check=True)
 for args in ([str(d/'verify'),'7',str(d/'contexts7.txt'),str(d/'reps7.txt')],[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(d/'controls.py')]):
  r=subprocess.run(args,check=True,capture_output=True);sys.stdout.buffer.write(r.stdout)
