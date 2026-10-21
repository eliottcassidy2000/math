#!/usr/bin/env python3
"""Authenticate this audit and replay its independent normal/-O checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=HERE.parent/'session18-zero-heavy')
    p.add_argument('--repo',type=Path,default=HERE.parent/'math')
    p.add_argument('--authenticate-only',action='store_true',help='Authenticate artifacts and inputs without recomputation.')
    args=p.parse_args()
    manifest=json.loads((HERE/'manifest.json').read_text())
    for name,digest in manifest['files'].items():
        need(sha((HERE/name).read_bytes())==digest,'audit artifact '+name)
    source_manifest=args.source/'manifest.json'
    need(sha(source_manifest.read_bytes())==manifest['source_manifest_sha256'],'frozen source manifest')
    src=json.loads(source_manifest.read_text())
    for name,digest in src['files'].items():
        need(sha((args.source/name).read_bytes())==digest,'source artifact '+name)
    for name,digest in src['dependencies'].items():
        need(sha((args.repo/name).read_bytes())==digest,'dependency '+name)
    expected=(HERE/'receipt.normal.json').read_bytes()
    need(expected==(HERE/'receipt.optimized.json').read_bytes(),'included normal/optimized identity')
    if args.authenticate_only:
        print(json.dumps({'status':'PASS authentication only; no recomputation','audit_files':len(manifest['files']),
                          'source_files':len(src['files']),'dependencies':len(src['dependencies'])},indent=2,sort_keys=True))
        return
    def run(optimized):
        command=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'independent_audit.py'),
                  '--source',str(args.source.resolve()),'--repo',str(args.repo.resolve())]
        done=subprocess.run(command,capture_output=True)
        need(done.returncode==0,'independent execution failed: '+done.stderr.decode())
        need(done.stdout==expected,'independent receipt mismatch, optimized='+str(optimized))
        return dict(optimized=optimized,bytes=len(done.stdout),sha256=sha(done.stdout),returncode=done.returncode)
    with ThreadPoolExecutor(max_workers=2) as pool:
        runs=list(pool.map(run,(False,True)))
    result=dict(status='PASS',frozen_manifest_sha256=manifest['source_manifest_sha256'],runs=runs)
    need(result==json.loads((HERE/'verification.json').read_text()),'verification receipt')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
