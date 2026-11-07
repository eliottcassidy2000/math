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
        jobs=(
            ('shells/verify.py','shells/verification.json',False),
            ('shells/probe_monotonicity.py','shells/monotonicity.json',False),
            ('doubling/verify.py','doubling/RESULT.json',True),
            ('plane/oddplane.py','plane/checks.json',True),
            ('audit/audit.py','audit/RESULT.json',True),
            ('audit/doubling_audit.py','audit/DOUBLING_RESULT.json',True),
            ('audit/selector_audit.py','audit/SELECTOR_RESULT.json',True),
            ('audit/plane_audit.py','audit/PLANE_RESULT.json',True),
        )
        for script,expected,writes_receipt in jobs:
            vals=[]
            for flag in ([], ['-O']):
                p=subprocess.run([sys.executable,*flag,str(dst/script)],cwd=dst,text=True,capture_output=True,check=True)
                value=json.loads((dst/expected).read_text() if writes_receipt else p.stdout)
                value.pop('elapsed_seconds',None)
                vals.append(value)
            need(vals[0]==vals[1], script+' differs under optimization')
            saved=json.loads((ROOT/expected).read_text());saved.pop('elapsed_seconds',None)
            need(vals[0]==saved, script+' differs from saved receipt')
            print('PASS',script)
        need((dst/'plane/oddplane.svg').read_bytes()==(ROOT/'plane/oddplane.svg').read_bytes(),'diagram regeneration')
    print('PASS hashes, normal/optimized replay, and saved exact receipts')

if __name__=='__main__':
    main()
