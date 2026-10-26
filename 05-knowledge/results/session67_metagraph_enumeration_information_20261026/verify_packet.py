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
        if 'git_blob' in pin:
            require(hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==pin['git_blob'],'local source Git blob mismatch')
        if repo:
            frozen=subprocess.check_output(['git','-C',str(repo),'show',pin['commit']+':'+pin['path']])
            require(frozen==raw,'repository source pin mismatch')
    native_checks=0
    def native(path, wanted):
        nonlocal native_checks
        require(hashlib.sha256(path.read_bytes()).hexdigest()==wanted,'native manifest mismatch: '+str(path))
        native_checks+=1
    for lane, filename in [('session67-local-audit','DELIVERABLES.json'),('session67-metagraph','artifact_manifest.json'),('session67-recursion','manifest.json')]:
        file=root/lane/filename
        if not file.exists():continue
        m=json.loads(file.read_text())
        if filename=='DELIVERABLES.json':rows=[(r['path'],r['sha256']) for r in m['files']]
        elif filename=='artifact_manifest.json':rows=[(k,v['sha256']) for k,v in m.items()]
        else:rows=list(m['files'].items())
        for name,wanted in rows:native(root/lane/name,wanted)
    for lane in ['session67-information','session67-cycle-audit','session67-incoming','session67-incoming/klein-packet']:
        for sums in (root/lane).glob('*SHA256SUMS'):
            for line in sums.read_text().splitlines():
                if line.strip():
                    wanted,name=line.split(maxsplit=1);native(root/lane/name.lstrip('*'),wanted)
    return {'status':'PASS','files':len(expected),'exact_file_set':True,'inherited_source_pins':len(pins),'repository_source_bytes_checked':bool(repo),'native_manifest_hashes_checked':native_checks}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args()
    print(json.dumps(verify(a.root,a.repo),sort_keys=True))
