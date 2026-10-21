#!/usr/bin/env python3
"""Portable artifact authentication and independent normal/-O replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise RuntimeError(message)

def sha(data):return hashlib.sha256(data).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=HERE.parent/'session18-zero-heavy-extension')
    p.add_argument('--base',type=Path,default=HERE.parent/'session18-zero-heavy')
    p.add_argument('--repo',type=Path,default=HERE.parent/'math')
    p.add_argument('--authenticate-only',action='store_true')
    args=p.parse_args();m=json.loads((HERE/'manifest.json').read_text())
    for name,digest in m['files'].items():need(sha((HERE/name).read_bytes())==digest,'audit artifact '+name)
    for folder,expected in [(args.source,m['extension_manifest_sha256']),(args.base,m['base_manifest_sha256'])]:
        data=(folder/'manifest.json').read_bytes();need(sha(data)==expected,'frozen source manifest')
        source=json.loads(data)
        for name,digest in source['files'].items():need(sha((folder/name).read_bytes())==digest,'frozen source '+name)
        for name,digest in source['dependencies'].items():
            legacy='05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md'
            target=HERE/'frozen_base_zero_heavy.md' if folder==args.base and name==legacy else args.repo/name
            need(sha(target.read_bytes())==digest,'dependency '+name)
    expected=(HERE/'receipt.normal.json').read_bytes()
    need(expected==(HERE/'receipt.optimized.json').read_bytes(),'included receipt identity')
    if args.authenticate_only:
        print(json.dumps(dict(status='PASS authentication only; no recomputation',audit_files=len(m['files']),
                              historical_base_source_capture=True),indent=2,sort_keys=True));return
    def run(optimized):
        cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'independent_extension_audit.py'),
             '--source',str(args.source.resolve()),'--base',str(args.base.resolve()),'--repo',str(args.repo.resolve())]
        result=subprocess.run(cmd,capture_output=True)
        need(result.returncode==0,'independent checker: '+result.stderr.decode())
        need(result.stdout==expected,'independent receipt mismatch')
        return dict(optimized=optimized,returncode=result.returncode,bytes=len(result.stdout),sha256=sha(result.stdout))
    with ThreadPoolExecutor(max_workers=2) as pool:runs=list(pool.map(run,(False,True)))
    record=dict(status='PASS',extension_manifest_sha256=m['extension_manifest_sha256'],
                base_manifest_sha256=m['base_manifest_sha256'],runs=runs)
    need(record==json.loads((HERE/'verification.json').read_text()),'verification receipt identity')
    print(json.dumps(record,indent=2,sort_keys=True))

if __name__=='__main__':main()
