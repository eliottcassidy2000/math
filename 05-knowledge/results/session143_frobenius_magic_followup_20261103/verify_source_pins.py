#!/usr/bin/env python3
"""Verify frozen snapshots offline; optionally check all commit/path memberships."""
import argparse
from hashlib import sha1, sha256
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def check_bytes(data, row):
    require(len(data) == row['bytes'], 'source byte count: '+row['path'])
    require(sha256(data).hexdigest() == row['sha256'], 'source SHA256: '+row['path'])
    require(sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() == row['git_blob'],
            'source Git blob: '+row['path'])

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, help='optional local Git checkout for exact published commit/path membership')
    args=parser.parse_args()
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())['repository_sources']
    snapshots=0
    for row in pins:
        for name in row['snapshots']:
            path=(ROOT/name).resolve()
            require(ROOT in path.parents and not (ROOT/name).is_symlink(), 'unsafe snapshot path')
            check_bytes(path.read_bytes(), row)
            snapshots+=1
        if args.repo:
            data=subprocess.check_output(['git','-C',str(args.repo),'show',row['commit']+':'+row['path']])
            check_bytes(data,row)
    print(json.dumps({'status':'PASS','source_identities':len(pins),'snapshots_verified':snapshots,
                      'commit_path_membership_verified':bool(args.repo),
                      'scope':'Offline byte identity only unless --repo is supplied; metadata-only contextual pins are not replayed as mathematics.'},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
