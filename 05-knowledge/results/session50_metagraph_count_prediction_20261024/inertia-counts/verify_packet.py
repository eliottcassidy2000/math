#!/usr/bin/env python3
"""Verify byte pins and replay ONLY a manifest-listed fresh-directory packet."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import json
import os
import shutil
import subprocess
import sys


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    root=Path(__file__).resolve().parent
    manifest=json.loads((root/'MANIFEST.json').read_text())
    with TemporaryDirectory(prefix='session50-inertia-manifest-') as tmp:
        fresh=Path(tmp)
        for entry in manifest['files']:
            rel=Path(entry['path'])
            require(not rel.is_absolute() and '..' not in rel.parts,'unsafe manifest path')
            source=root/rel
            data=source.read_bytes()
            require(len(data)==entry['bytes'],'manifest size mismatch: '+str(rel))
            require(sha(data)==entry['sha256'],'manifest hash mismatch: '+str(rel))
            target=fresh/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(data)
        runs=[]
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal'
            args=[sys.executable]+(['-O'] if optimized else [])+['check.py','--output','fresh-'+mode+'.json']
            completed=subprocess.run(args,cwd=fresh,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False,
                                     env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
            require(completed.returncode==0,'fresh replay failed: '+completed.stderr.decode())
            actual=(fresh/('fresh-'+mode+'.json')).read_bytes()
            require(actual==(fresh/(mode+'.json')).read_bytes(),'fresh JSON mismatch: '+mode)
            require(completed.stdout==(fresh/(mode+'.out')).read_bytes(),'fresh stdout mismatch: '+mode)
            require(not completed.stderr,'unexpected replay stderr: '+mode)
            runs.append({'mode':mode,'json_sha256':sha(actual),'stdout_sha256':sha(completed.stdout),'byte_identical':True})
        print(json.dumps({'manifest_only_copy':True,'files_verified':len(manifest['files']),
                          'normal_and_optimized_replayed':True,'runs':runs},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
