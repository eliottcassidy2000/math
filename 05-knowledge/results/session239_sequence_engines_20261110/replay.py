#!/usr/bin/env python3
"""Replay all load-bearing correctness checks; archived benchmark files untouched."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
SCRIPTS = ['audit/verify.py', 'audit/strong_gap_audit.py', 'audit/palette_audit.py',
           'proof_bridges/verify.py', 'proof_bridges/prime_gap.py',
           'proof_bridges/c7_low_fibres.py', 'precision/precision.py', 'verify_cli.py']
results=[]
for rel in SCRIPTS:
    cmd=[sys.executable] + (['-O'] if sys.flags.optimize else []) + [str(BASE/rel)]
    start=time.perf_counter()
    p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError((rel,p.returncode,p.stdout,p.stderr))
    data=json.loads(p.stdout)
    results.append({'script':rel,'pass':True,'seconds':time.perf_counter()-start,
                    'source_sha256':hashlib.sha256((BASE/rel).read_bytes()).hexdigest(),
                    'result':data})
print(json.dumps({'status':'PASS','optimized':bool(sys.flags.optimize),'checks':results},indent=2))
