#!/usr/bin/env python3
"""Verify portable hashes and replay the exact census in ordinary and -O modes."""
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text())
for item in manifest['files']:
    p=root/item['path']
    data=p.read_bytes()
    if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
        raise RuntimeError('Manifest mismatch: '+item['path'])
expected=(root/'classification.json').read_bytes()
for flags in ([],['-O']):
    actual=subprocess.check_output([sys.executable,*flags,str(root/'check_tournaments.py')],cwd=root)
    if actual!=expected:
        raise RuntimeError('Replay mismatch: '+repr(flags))
subprocess.run([sys.executable,str(root/'catalogue_crosscheck.py')],check=True,cwd=root)
print('PASS: portable manifest and exact normal/-O replay; '+str(json.loads(expected)['checks_passed'])+' checks per run')
