#!/usr/bin/env python3
"""Freeze exactly the producer manifest, replay, and test hostile mutations.
Run: python3 verify_replay.py [candidate directory].
No producer files are edited. All mutation cases live in temporary copies.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
FROZEN=HERE/'frozen_candidate'
SOURCE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else (FROZEN if FROZEN.exists() else HERE.parent/'session52-fifteen-kernels')

def require(x,message):
    if not x:
        raise RuntimeError(message)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

manifest=json.loads((SOURCE/'delivery-manifest.json').read_text())
reuse=FROZEN.exists()
FROZEN.mkdir(exist_ok=True)
listed=[]
for item in manifest['files']:
    p=Path(item['path'])
    require(not p.is_absolute() and '..' not in p.parts,'unsafe manifest path')
    require(digest(SOURCE/p)==item['sha256'],'source digest mismatch: '+str(p))
    require((SOURCE/p).stat().st_size==item['bytes'],'source size mismatch: '+str(p))
    require(not (SOURCE/p).is_symlink(),'source is symlink: '+str(p))
    if reuse:
        require(digest(FROZEN/p)==item['sha256'],'existing frozen file mismatch')
    else:
        (FROZEN/p).parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(SOURCE/p,FROZEN/p)
    listed.append(p.as_posix())
if reuse:
    require(digest(SOURCE/'delivery-manifest.json')==digest(FROZEN/'delivery-manifest.json'),'existing frozen manifest mismatch')
else:
    shutil.copyfile(SOURCE/'delivery-manifest.json',FROZEN/'delivery-manifest.json')
expected_files=sorted(listed+['delivery-manifest.json'])
require(sorted(p.relative_to(FROZEN).as_posix() for p in FROZEN.rglob('*') if p.is_file())==expected_files,'exact manifest file set')
require(not any(p.is_symlink() for p in FROZEN.rglob('*')),'frozen tree symlink')

env=dict(os.environ)
env.pop('PYTHONPATH',None)
def run(root,optimized,out):
    command=[sys.executable,'-I']+(['-O'] if optimized else [])+[str(root/'check.py'),'--receipt',str(out)]
    return subprocess.run(command,cwd='/tmp',env=env,capture_output=True)

with tempfile.TemporaryDirectory(prefix='independent-fifteen-realcopy-') as t:
    relocated=Path(t)/'packet'
    shutil.copytree(FROZEN,relocated)
    require(not any(p.is_symlink() for p in relocated.rglob('*')),'relocated symlink')
    require(sorted(p.relative_to(relocated).as_posix() for p in relocated.rglob('*') if p.is_file())==expected_files,'relocation extra files')
    normal=run(relocated,False,HERE/'replay_receipt.json')
    optimized=run(relocated,True,HERE/'replay_receipt_optimized.json')
for result in [normal,optimized]:
    require(result.returncode==0,result.stderr.decode())
require(normal.stdout==optimized.stdout==(FROZEN/'receipt.json').read_bytes(),'normal/-O/original mismatch')
require(normal.stdout==(HERE/'replay_receipt.json').read_bytes()==(HERE/'replay_receipt_optimized.json').read_bytes(),'stdout/receipt mismatch')
(HERE/'replay_normal.out').write_bytes(normal.stdout)
(HERE/'replay_optimized.out').write_bytes(optimized.stdout)

def reseal(root,name):
    p=root/'source-manifest.json'
    m=json.loads(p.read_text())
    for entry in m['files']:
        if entry['path']==name:
            entry['sha256']=digest(root/name)
    p.write_text(json.dumps(m,indent=2)+'\n')

cases=[]
with tempfile.TemporaryDirectory(prefix='independent-fifteen-hostiles-') as t:
    for name,expected in [('tampered_input','input digest mismatch'),('missing_input','input missing'),('resealed_wrong_orbits','published Bolza kernel orbit partition'),('resealed_overlapping_cycles','disjoint cycle transcription')]:
        root=Path(t)/name
        shutil.copytree(FROZEN,root)
        target=root/'inputs/bolza-kernel-action.json'
        if name=='tampered_input':
            target.write_bytes(target.read_bytes()+b'\n')
        elif name=='missing_input':
            target.unlink()
        else:
            data=json.loads(target.read_text())
            if name=='resealed_wrong_orbits':
                data['expected_orbits'][0][0]=5
            else:
                data['generators_cycles']['sigma'].append([8,10])
            target.write_text(json.dumps(data,indent=2)+'\n')
            reseal(root,'inputs/bolza-kernel-action.json')
        outcomes=[]
        for opt in [False,True]:
            result=run(root,opt,root/'hostile_receipt.json')
            require(result.returncode!=0,'hostile accepted: '+name)
            require(expected.encode() in result.stderr,'wrong rejection reason: '+name+' '+result.stderr.decode())
            require(not (root/'hostile_receipt.json').exists(),'hostile emitted a success receipt')
            outcomes.append({'optimized':opt,'exit_code':result.returncode,'matched_reason':expected})
        cases.append({'case':name,'outcomes':outcomes})

# Original source and frozen bytes must remain unchanged after all replay tests.
for item in manifest['files']:
    require(digest(SOURCE/item['path'])==item['sha256'],'source mutated during audit')
    require(digest(FROZEN/item['path'])==item['sha256'],'frozen candidate mutated during audit')
result={'status':'PASS','declared_files_copied':len(listed),'actual_copy_not_symlinks':True,
        'manifest_only_copy':True,'isolated_normal_and_optimized_identical':True,
        'candidate_manifest_sha256':digest(FROZEN/'delivery-manifest.json'),
        'candidate_checker_sha256':digest(FROZEN/'check.py'),
        'candidate_readme_sha256':digest(FROZEN/'README.md'),
        'receipt_sha256':hashlib.sha256(normal.stdout).hexdigest(),
        'candidate_gates':json.loads(normal.stdout)['checks'],
        'hostile_cases':cases,'source_and_frozen_bytes_unchanged':True}
print(json.dumps(result,sort_keys=True,indent=2))
