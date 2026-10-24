#!/usr/bin/env python3
"""Verify the outer manifest, copy only its payloads, and replay in isolation.

The outer manifest is the package root. It does not hash itself. The temporary
package receives exactly that manifest plus its listed payload files, no cache,
source checkout, probe script, inherited note, or unlisted producer artifact.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
MANIFEST='AUDIT_MANIFEST.json'

def require(ok,message):
    if not ok:raise ValueError(message)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(root):
    doc=json.loads((root/MANIFEST).read_text())
    listed={r['path']:r for r in doc['files']}
    require(len(listed)==len(doc['files']),'duplicate outer manifest entry')
    for name,row in listed.items():
        rel=Path(name)
        require(not rel.is_absolute() and '..' not in rel.parts,'unsafe manifest path')
        require('__pycache__' not in rel.parts and rel.suffix!='.pyc','cache in deliverable')
        p=root/rel
        require(p.is_file(),'missing outer payload: '+name)
        require(p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],'outer hash mismatch: '+name)
    needed={'independent_check.py','independent_results.json','verify_replay.py',
            'frozen_candidate/MANIFEST.json','frozen_candidate/results.json',
            'frozen_candidate/inputs/prime_offsets_certificates.json'}
    require(needed<=set(listed),'required audit dependency absent from outer manifest')
    inner=root/'frozen_candidate/MANIFEST.json'
    require(digest(inner)==doc['candidate_manifest_sha256'],'wrong frozen candidate manifest')
    candidate=json.loads(inner.read_text())
    for row in candidate['files']:
        name='frozen_candidate/'+row['path']
        require(name in listed,'unlisted nested candidate payload: '+name)
        require(listed[name]['sha256']==row['sha256'] and listed[name]['bytes']==row['bytes'],'inner/outer disagreement: '+name)
    return doc

def main():
    doc=verify(ROOT)
    commands=[
        ('candidate normal',['frozen_candidate/check.py'],'replay_normal.out'),
        ('candidate optimized',['-O','frozen_candidate/check.py'],'replay_optimized.out'),
        ('candidate hostiles normal',['frozen_candidate/test_rejections.py'],'rejections_normal.out'),
        ('candidate hostiles optimized',['-O','frozen_candidate/test_rejections.py'],'rejections_optimized.out'),
        ('independent normal',['independent_check.py'],'independent_normal.out'),
        ('independent optimized',['-O','independent_check.py'],'independent_optimized.out')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='fortunate_square_manifest_replay_') as tmp:
        work=Path(tmp);package=work/'package';package.mkdir();unrelated=work/'unrelated-cwd';unrelated.mkdir()
        shutil.copyfile(ROOT/MANIFEST,package/MANIFEST)
        for row in doc['files']:
            src=ROOT/row['path'];dst=package/row['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        expected={MANIFEST}|{r['path'] for r in doc['files']}
        actual={p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()}
        require(actual==expected,'copy contains unlisted files or lacks payloads')
        verify(package)
        for label,args,saved in commands:
            argv=[sys.executable,'-B']+[arg if arg.startswith('-') else str(package/arg) for arg in args]
            result=subprocess.run(argv,cwd=unrelated,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
            require(result.returncode==0,label+' failed: '+result.stderr.decode(errors='replace'))
            require(result.stdout==(package/saved).read_bytes(),label+' stdout differs from saved artifact')
            runs.append({'check':label,'status':'PASS','stdout_sha256':hashlib.sha256(result.stdout).hexdigest()})
        verify(package)
        actual_after={p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()}
        require(actual_after==expected,'replay introduced unlisted artifacts')
    out={'status':'PASS','package_payload_files':len(doc['files']),
         'copied_files_including_root_manifest':len(doc['files'])+1,
         'copy_policy':'Only the root manifest and its final listed payloads; no caches or unlisted dependencies.',
         'working_directory':'Unrelated empty directory, outside the copied package.',
         'all_listed_hashes_match_before_and_after':True,
         'candidate_manifest_sha256':doc['candidate_manifest_sha256'],
         'independent_results_sha256':digest(ROOT/'independent_results.json'),
         'runs':runs}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':main()
