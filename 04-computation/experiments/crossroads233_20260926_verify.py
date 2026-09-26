"""Reproduce all four 233-session programs in normal and optimized modes."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREFIX = 'crossroads233_20260926_'
CASES = ['carry', 'flow', 'graph', 'scale_audit']


def one(stem):
    source = Path('04-computation/experiments') / (PREFIX+stem+'.py')
    output = Path('05-knowledge/results') / (PREFIX+stem+'.out')
    results = []
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    for flags in ([], ['-O']):
        p = subprocess.run([sys.executable, *flags, str(source)], cwd=ROOT,
                           env=env, capture_output=True, text=True,
                           encoding='utf-8', timeout=240)
        if p.returncode:
            raise RuntimeError(f'{stem}: {p.stderr}')
        results.append(p.stdout.replace('\r\n', '\n'))
    if results[0] != results[1]:
        raise RuntimeError(f'{stem}: normal/-O mismatch')
    (ROOT/output).write_text(results[0], encoding='utf-8', newline='\n')
    print('PASS normal/-O:', stem, flush=True)
    return {'source': source.as_posix(), 'output': output.as_posix(),
            'normal_equals_optimized': True}


def main():
    with ThreadPoolExecutor(max_workers=3) as pool:
        cases = [f.result() for f in as_completed([pool.submit(one,s) for s in CASES])]
    paths = []
    for folder in ('04-computation/experiments','05-knowledge/results'):
        paths.extend(p for p in (ROOT/folder).glob(PREFIX+'*')
                     if p.suffix in ('.py','.md','.out'))
    for number in (4492,4493):
        paths.extend((ROOT/'01-canon/theorems').glob(f'THM-{number}-*.md'))
    artifacts = []
    for p in sorted(paths):
        raw = p.read_bytes().replace(b'\r\n',b'\n')
        p.write_bytes(raw)
        artifacts.append({'path':p.relative_to(ROOT).as_posix(),
                          'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
    manifest = {'status':'REPRODUCTION VERIFIED; scope is stated in the audited proofs',
                'generated_utc':datetime.now(timezone.utc).isoformat(),
                'python':sys.version.split()[0], 'hash_basis':'raw UTF-8 LF bytes',
                'cases':sorted(cases,key=lambda x:x['source']), 'artifacts':artifacts}
    target=ROOT/'05-knowledge/results'/(PREFIX+'manifest.json')
    target.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(f'PASS: {len(cases)} programs in both modes; {len(artifacts)} artifact hashes')


if __name__ == '__main__':
    main()
