#!/usr/bin/env python3
"""Verify exact curated packet bytes; no dependencies outside the standard library."""
import argparse, hashlib, json, subprocess
from pathlib import Path

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def verify(root, repo=None):
    manifest=json.loads((root/'FILE_MANIFEST.json').read_text())
    expected={row['path']:row for row in manifest['files']}
    require(len(expected)==len(manifest['files']),'duplicate manifest paths')
    actual=set()
    for p in root.rglob('*'):
        require(not p.is_symlink(),'symbolic link: '+str(p))
        if p.is_file():actual.add(p.relative_to(root).as_posix())
    require(actual==set(expected)|{'FILE_MANIFEST.json'},'exact file set mismatch')
    for name,row in expected.items():
        require(not Path(name).is_absolute() and '..' not in Path(name).parts,'unsafe path')
        raw=(root/name).read_bytes()
        require(len(raw)==row['bytes'],'size mismatch: '+name)
        require(hashlib.sha256(raw).hexdigest()==row['sha256'],'SHA-256 mismatch: '+name)
    pins=json.loads((root/'PROVENANCE.json').read_text())['inherited_inputs']
    for pin in pins:
        raw=(root/pin['local_path']).read_bytes()
        require(hashlib.sha256(raw).hexdigest()==pin['sha256'],'local source pin mismatch')
        if repo:
            frozen=subprocess.check_output(['git','-C',str(repo),'show',pin['commit']+':'+pin['path']])
            require(frozen==raw,'repository source pin mismatch')
    return {'status':'PASS','files':len(expected),'exact_file_set':True,'inherited_source_pins':len(pins),'repository_source_bytes_checked':bool(repo)}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args()
    print(json.dumps(verify(a.root,a.repo),sort_keys=True))
