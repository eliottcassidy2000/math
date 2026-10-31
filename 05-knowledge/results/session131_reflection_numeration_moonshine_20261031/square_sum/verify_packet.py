#!/usr/bin/env python3
"""Verify the explicit release manifest and content/source bindings."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def require(x,message):
    if not x:raise ValueError(message)

def main():
    names=json.loads((HERE/'RELEASE_FILES.json').read_text())
    hashes=json.loads((HERE/'MANIFEST.json').read_text())
    expected=set(names)-{'MANIFEST.json'}
    require(set(hashes)==expected,'manifest whitelist mismatch')
    for name,digest in hashes.items():
        require(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,
                'hash mismatch: '+name)
    a=(HERE/'verification.normal.json').read_bytes()
    b=(HERE/'verification.optimized.json').read_bytes()
    require(a==b,'normal/optimized mismatch')
    receipt=json.loads(a)
    for name,digest in receipt['input_sha256'].items():
        require(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,
                'receipt source mismatch: '+name)
    require(receipt['status']=='PASS','receipt did not pass')
    print(json.dumps({'status':'PASS','release_files':len(names),
                      'checked_conditions':receipt['checked_conditions']},sort_keys=True))
if __name__=='__main__':main()
