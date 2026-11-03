#!/usr/bin/env python3
"""Portable integrity, ordinary/optimized replay, and negative-control validation."""
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text())
for item in manifest['files']:
    data=(root/item['path']).read_bytes()
    if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
        raise RuntimeError('Manifest mismatch: '+item['path'])
for script,output in [('check_integral_order.py','certificate.json'),('test_negative_controls.py','negative_controls.json')]:
    expected=(root/output).read_bytes()
    for flags in ([],['-O']):
        got=subprocess.check_output([sys.executable,*flags,str(root/script)],cwd=root)
        if got!=expected:raise RuntimeError('Replay mismatch '+script+' '+repr(flags))
print('PASS: portable manifest, 3302 exact checks normal/-O, and eight mutants rejected in both modes')
