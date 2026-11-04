#!/usr/bin/env python3
"""Replay the declared exact certificates in normal and optimized Python."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOBS=(
    ('interior/check.py','interior/check.json'),
    ('interior/graded.py','interior/graded.json'),
    ('interfaces/check.py','interfaces/check.json'),
    ('interfaces/exceptional.py','interfaces/exceptional.json'),
    ('interfaces/audit_transport.py','interfaces/audit_transport.json'),
    ('conditional/check.py','conditional/check.json'),
    ('conditional/audit_interfaces.py','conditional/audit_interfaces.json'),
    ('transport/check.py','transport/check.json'),
    ('transport/audit_interior.py','transport/audit_interior.json'),
)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def run(job):
    script,receipt=job
    expected=(HERE/receipt).read_bytes()
    if json.loads(expected).get('status')!='PASS':
        raise RuntimeError('Receipt is not PASS: '+receipt)
    for optimized in (False,True):
        command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(HERE/script)]
        result=subprocess.run(command,cwd=ROOT,capture_output=True,timeout=240)
        if result.returncode:
            raise RuntimeError(script+': '+result.stderr.decode(errors='replace')[-4000:])
        if result.stdout!=expected:
            raise RuntimeError(script+' output mismatch'+(' under -O' if optimized else ''))
    return {'script':script,'script_sha256':sha((HERE/script).read_bytes()),
            'receipt':receipt,'receipt_sha256':sha(expected),
            'normal_and_optimized_byte_identical':True}

def main():
    pins=json.loads((HERE/'sources.json').read_text())
    for item in pins['sources']:
        if sha((ROOT/item['path']).read_bytes())!=item['sha256']:
            raise RuntimeError('Changed direct source: '+item['path'])
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs=list(pool.map(run,JOBS))
    print(json.dumps({'status':'PASS','jobs':jobs,
        'direct_source_pins_checked':len(pins['sources']),
        'scope':'Nine scoped exact producers and auditors; all normal/optimized outputs match saved receipts. Universal source acceptance remains OPEN.',
        'runner_sha256':sha(Path(__file__).read_bytes()),
        'sources_sha256':sha((HERE/'sources.json').read_bytes())},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
