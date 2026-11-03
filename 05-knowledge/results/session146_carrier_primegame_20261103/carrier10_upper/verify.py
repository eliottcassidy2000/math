#!/usr/bin/env python3
"""Reconstruct chunks and run the independent exact integer audit."""
from pathlib import Path
from tempfile import TemporaryDirectory
from hashlib import sha256
import subprocess
H=Path(__file__).resolve().parent
raw=(H/'kernel.part00').read_bytes()+(H/'kernel.part01').read_bytes()
expected='91c68094d7c5d74467a9de29e949b6bbabe844edd0e2c0193bca9d5b241522fb'
if sha256(raw).hexdigest()!=expected:raise RuntimeError('kernel chunk hash mismatch')
with TemporaryDirectory() as d:
 d=Path(d);(d/'kernel.txt').write_bytes(raw)
 subprocess.run(['g++','-O3','-std=c++17',str(H/'independent.cpp'),'-o',str(d/'audit')],check=True)
 p=subprocess.run([str(d/'audit'),str(H/'reps4.txt'),str(d/'kernel.txt')],check=True,text=True,capture_output=True)
 reference=(H/'independent.log').read_text()
 if p.stdout!=reference:raise RuntimeError('independent audit output mismatch')
 print('joined kernel SHA-256',expected)
 print(p.stdout,end='')
