#!/usr/bin/env python3
"""Portable exact replay of the audited session packet."""
from pathlib import Path
import hashlib, json, subprocess, sys
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def main():
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        h,name=line.split('  ',1)
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,'hash '+name)
    expected=json.loads((ROOT/'bridges/VERIFY.json').read_text())
    for flag in ([], ['-O']):
        p=subprocess.run([sys.executable,'-B',*flag,str(ROOT/'bridges/verify.py')],capture_output=True,text=True,check=True)
        need(json.loads(p.stdout)==expected,'bridge replay '+str(flag))
    need(expected==json.loads((ROOT/'bridges/VERIFY_OPT.json').read_text()),'optimized receipt')
    need(expected==json.loads((ROOT/'audit/BRIDGE_VERIFY.json').read_text()),'audit receipt')
    print('PASS hashes and normal/optimized exact bridge replay')

if __name__=='__main__': main()
