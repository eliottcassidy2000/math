#!/usr/bin/env python3
"""Replay the complete bounded search without writing into frozen packet bytes."""
import os,shutil,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='session73-full-depth-replay-') as temp:
 script=Path(temp)/'search_depth.py';shutil.copy2(root/'session73-ports-depth/search_depth.py',script)
 env={k:v for k,v in os.environ.items() if k not in ('PYTHONPATH','PYTHONHOME')};env['PYTHONDONTWRITEBYTECODE']='1'
 p=subprocess.run([sys.executable,'-I','-B',*(['-O'] if sys.flags.optimize else []),str(script)],cwd=temp,capture_output=True,env=env,timeout=1100)
 if p.returncode:raise RuntimeError(p.stderr.decode())
 if (Path(temp)/'search_result.json').read_bytes()!=p.stdout:raise RuntimeError('search output file differs from stdout')
 if p.stdout!=(root/'session73-ports-depth/receipt.json').read_bytes():raise RuntimeError('full depth search receipt differs')
 sys.stdout.buffer.write(p.stdout)
