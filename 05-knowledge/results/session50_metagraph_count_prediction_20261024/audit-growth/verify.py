#!/usr/bin/env python3
"""Manifest-closed verification with optional independent recomputation."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def load(name):
    return json.loads((ROOT/name).read_text())

def normalized(receipt):
    return {k:v for k,v in receipt.items() if k!='elapsed_seconds'}

def main(deep):
    manifest=load('manifest.json')
    expected=manifest['files']
    actual=set()
    for path in ROOT.rglob('*'):
        need(not path.is_symlink(),'symlink in sealed packet')
        if path.is_file():
            actual.add(path.relative_to(ROOT).as_posix())
    need(actual==set(expected)|{'manifest.json'},'manifest file closure')
    for name,sha in expected.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'digest mismatch: '+name)
    normal=load('replay-normal-order8.json'); optimized=load('replay-optimized-order8.json')
    need(normal==optimized,'full normal/optimized receipts')
    need(normal['total_checks']==2641571,'full replay explicit check count')
    need((ROOT/'replay-normal-order8.out').read_bytes()==(ROOT/'replay-optimized-order8.out').read_bytes(),
         'full normal/optimized displayed output')
    independent=load('independent-receipt.json')
    need(normalized(independent)==normalized(load('independent-optimized.json')),
         'independent normal/optimized receipts')
    need(independent['total_checks']==10503,'independent explicit check count')
    pinned=hashlib.sha256((ROOT/'source/transfer.py.txt').read_bytes()).hexdigest()
    need(independent['producer_sha256']==pinned,'independent producer pin')
    for mode in ('normal','optimized'):
        timing=load('replay-'+mode+'-order8.time.json')
        need(timing['returncode']==0 and timing['source_sha256']==pinned,'replay source and successful exit')
    for pin in load('source-pins.json')['sources']:
        need(hashlib.sha256((ROOT/pin['archive']).read_bytes()).hexdigest()==pin['sha256'],
             'inherited archive pin')
    if deep:
        with tempfile.TemporaryDirectory(prefix='growth-audit-check-') as tmp:
            receipt=Path(tmp)/'receipt.json'
            command=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(ROOT/'independent_audit.py'),'--receipt',str(receipt)]
            result=subprocess.run(command,cwd=tmp,text=True,capture_output=True)
            need(result.returncode==0,'independent rerun failed: '+result.stderr)
            need(normalized(json.loads(receipt.read_text()))==normalized(independent),'independent rerun receipt')
            need(result.stdout==(ROOT/'independent.out').read_text(),'independent rerun displayed output')
    print(json.dumps({'status':'PASS','manifest_files':len(expected),'deep':deep,
                      'independent_checks':10503,'producer_checks':2641571},sort_keys=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--deep',action='store_true')
    main(parser.parse_args().deep)
