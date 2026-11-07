#!/usr/bin/env python3
"""Portable replay. Checks are active under -O; generated receipts stay temporary."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent

def need(p, message):
    if not p:
        raise RuntimeError(message)

def main():
    for row in (ROOT / 'SHA256SUMS').read_text().splitlines():
        digest, name = row.split('  ', 1)
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name)
    with tempfile.TemporaryDirectory(prefix='session222-') as d:
        dst = Path(d)/'packet'
        shutil.copytree(ROOT, dst)
        for script in ('shells/verify.py', 'shells/probe_monotonicity.py', 'audit/audit.py'):
            vals=[]
            for flag in ([], ['-O']):
                p=subprocess.run([sys.executable,*flag,str(dst/script)],cwd=dst,text=True,capture_output=True,check=True)
                value=json.loads(p.stdout)
                value.pop('elapsed_seconds',None)
                vals.append(value)
            need(vals[0]==vals[1], script+' differs under optimization')
            expected='shells/verification.json' if script=='shells/verify.py' else 'shells/monotonicity.json' if script=='shells/probe_monotonicity.py' else 'audit/RESULT.json'
            saved=json.loads((ROOT/expected).read_text());saved.pop('elapsed_seconds',None)
            need(vals[0]==saved, script+' differs from saved receipt')
            print('PASS',script)
    print('PASS hashes, normal/optimized replay, and saved exact receipts')

if __name__=='__main__':
    main()
