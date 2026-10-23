#!/usr/bin/env python3
"""Replay immutable inputs and exact normal/optimized receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
manifest=json.loads((HERE/'manifest.json').read_text())
for item in manifest['external_sources']:
    data=subprocess.check_output(['git','show',item['commit']+':'+item['path']],cwd=ROOT)
    if hashlib.sha256(data).hexdigest()!=item['sha256']:raise RuntimeError('source changed: '+item['path'])
for name,digest in manifest['packet_files'].items():
    if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('packet changed: '+name)
jobs=(('native_cube.py','cube_receipt.json'),('pointed_atlas.py','atlas_receipt.json'),
      ('hall_components.py','hall_receipt.json'),('stabilization.py','stabilization_receipt.json'),
      ('one_bit_probe.py','one_bit_receipt.json'),('audit.py','audit_receipt.json'))
for script,receipt in jobs:
    expected=json.loads((HERE/receipt).read_text())
    for flags in (['-B'],['-O','-B']):
        actual=json.loads(subprocess.check_output([sys.executable,*flags,str(HERE/script)],cwd=ROOT))
        if actual!=expected:raise RuntimeError('receipt mismatch: '+script+' '+str(flags))
    print(script+': normal and optimized exact receipts match')
print('PASS: pinned sources and all six exact replays; general statements require the accompanying proofs.')
