#!/usr/bin/env python3
"""Portable exact replay. No network, repository mutation or disabled assertions."""
from pathlib import Path
import hashlib, json, subprocess, sys
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def main():
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        h,name=line.split('  ',1)
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,'hash '+name)
    jobs=[
        ('bridges/verify.py','bridges/VERIFY.json',True),
        ('cutoff/verify_drift.py','cutoff/VERIFY.txt',False),
        ('cutoff/verify_blocks.py','cutoff/BLOCK_VERIFY.txt',False),
        ('threeadic/verify.py','threeadic/VERIFY.txt',False),
        ('threeadic/exact_higher.py','threeadic/HIGHER.txt',False),
        ('threeadic/exact_81.py','threeadic/HIGHER81.txt',False),
        ('audit/audit_mod3.py','audit/MOD3_VERIFY.txt',False),
        ('audit/audit_higher.py','audit/HIGHER_VERIFY.txt',False),
        ('audit/audit_threeadic.py','audit/THREEADIC_CHECKS.json',True),
    ]
    for script,receipt,as_json in jobs:
        expected=(ROOT/receipt).read_text()
        for flags in ([], ['-O']):
            p=subprocess.run([sys.executable,'-B',*flags,str(ROOT/script)],cwd=ROOT/script.split('/')[0],capture_output=True,text=True,check=True)
            need((json.loads(p.stdout)==json.loads(expected)) if as_json else p.stdout==expected,'replay '+script+' '+str(flags))
        print('PASS',script,flush=True)
    pairs=[
        ('bridges/VERIFY.json','bridges/VERIFY_OPT.json'),
        ('bridges/VERIFY.json','audit/BRIDGE_VERIFY.json'),
        ('audit/BRIDGE_VERIFY.json','audit/BRIDGE_OPT_VERIFY.json'),
        ('threeadic/VERIFY.txt','threeadic/VERIFY_OPT.txt'),
        ('threeadic/HIGHER.txt','threeadic/HIGHER_OPT.txt'),
        ('threeadic/HIGHER81.txt','threeadic/HIGHER81_OPT.txt'),
        ('audit/MOD3_VERIFY.txt','audit/MOD3_OPT_VERIFY.txt'),
        ('audit/HIGHER_VERIFY.txt','audit/HIGHER_OPT_VERIFY.txt'),
        ('audit/CUTOFF_VERIFY.txt','audit/CUTOFF_OPT_VERIFY.txt'),
        ('cutoff/BLOCK_VERIFY.txt','audit/CUTOFF_VERIFY.txt'),
        ('audit/THREEADIC_CHECKS.json','audit/THREEADIC_CHECKS_OPT.json'),
    ]
    for a,b in pairs:need((ROOT/a).read_bytes()==(ROOT/b).read_bytes(),'receipt pair '+a+' '+b)
    print('PASS all hashes, nine normal/optimized replays, and exact saved receipts')

if __name__=='__main__': main()
