#!/usr/bin/env python3
"""Verify this portable packet and replay its exact finite certificates."""
from hashlib import sha256
from pathlib import Path
import json,subprocess,sys
HERE=Path(__file__).resolve().parent
def need(ok,message):
 if not ok:raise RuntimeError(message)
def main():
 manifest=json.loads((HERE/'MANIFEST.json').read_text())
 actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
 need(actual==set(manifest),'manifest is the complete portable packet, excluding itself')
 for name,row in manifest.items():
  data=(HERE/name).read_bytes()
  need(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'],'manifest file '+name)
 for script in ('theory_check.py','native_check.py'):
  command=[sys.executable]
  if sys.flags.optimize:command.append('-O')
  command+=['-B',str(HERE/script)]
  subprocess.run(command,check=True,cwd=HERE)
 print('PASS complete portable manifest and both exact checkers')
if __name__=='__main__':main()
