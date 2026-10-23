#!/usr/bin/env python3
"""Replay frozen sources and receipts; never refresh expectations on failure."""
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
manifest=json.loads((HERE/'manifest.json').read_text())
for name,digest in manifest['files'].items():
    if sha256((HERE/name).read_bytes()).hexdigest()!=digest:
        raise RuntimeError('packet hash mismatch: '+name)
for row in manifest['sources']:
    data=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=ROOT)
    if sha256(data).hexdigest()!=row['sha256']:
        raise RuntimeError('immutable source hash mismatch: '+row['path'])
for row in manifest['working_inputs']:
    data=(ROOT/row['path']).read_bytes()
    if row.get('normalize_text'):data=data.replace(b'\r\n',b'\n')
    if sha256(data).hexdigest()!=row['sha256']:
        raise RuntimeError('working input hash mismatch: '+row['path'])
print('Frozen packet, immutable Git sources and working inputs: PASS',flush=True)
for name in ('routing','audit','electrical','ports','carry'):
    expected=(HERE/(name+'_receipt.json')).read_bytes()
    for mode in ([],['-O']):
        actual=subprocess.check_output([sys.executable,'-B',*mode,str(HERE/(name+'.py'))],cwd=ROOT)
        if actual!=expected:raise RuntimeError('receipt mismatch: '+name+' '+str(mode))
    print(name+': normal/-O exact replay PASS',flush=True)
print('All five finite receipts PASS; all-size claims remain scoped by their proofs.')
