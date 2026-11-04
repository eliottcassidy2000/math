#!/usr/bin/env python3
"""Replay the exact research certificates; benchmark timings are not replayed.

Each job runs normally and optimized, compares both byte streams with the
published receipt, and reports hashes. Independent audit universes remain
those declared in their own receipts; replay does not enlarge their scope.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOBS=(
    ('counting/check.py','counting/check.json'),
    ('counting/full_checkpoint.py','counting/full256.json'),
    ('bridge/check.py','bridge/check.json'),
    ('counting/bridge_audit.py','counting/bridge_audit.json'),
    ('words/check.py','words/check.json'),
    ('words/independent_audit.py','words/independent_audit.json'),
    ('coverage/check.py','coverage/check.json'),
    ('coverage/independent_audit.py','coverage/independent_audit.json'),
)

def digest(data): return hashlib.sha256(data).hexdigest()

def run_job(job):
    script,receipt=job;expected=(HERE/receipt).read_bytes()
    parsed=json.loads(expected)
    if parsed.get('status')!='PASS': raise RuntimeError('Non-PASS receipt: '+receipt)
    for optimized in (False,True):
        command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(HERE/script)]
        process=subprocess.run(command,cwd=ROOT,capture_output=True,timeout=240)
        if process.returncode:
            raise RuntimeError(script+' failed: '+process.stderr.decode(errors='replace')[-4000:])
        if process.stdout!=expected:
            raise RuntimeError(script+' did not reproduce '+receipt+(' with -O' if optimized else ''))
    return {'script':script,'script_sha256':digest((HERE/script).read_bytes()),
            'receipt':receipt,'receipt_sha256':digest(expected),
            'normal_and_optimized_byte_identical':True}

def main():
    pins=json.loads((HERE/'sources.json').read_text())
    for entry in pins['sources']:
        actual=digest((ROOT/entry['path']).read_bytes())
        if actual!=entry['sha256_current']:
            raise RuntimeError('Changed direct source: '+entry['path'])
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(run_job,JOBS))
    print(json.dumps({'status':'PASS','jobs':results,
        'direct_source_pins_checked':len(pins['sources']),
        'scope':'Eight exact jobs replayed in normal and optimized Python; timing measurements excluded.',
        'runner_sha256':digest(Path(__file__).read_bytes()),
        'sources_sha256':digest((HERE/'sources.json').read_bytes())},indent=2,sort_keys=True))

if __name__=='__main__': main()
