#!/usr/bin/env python3
"""Portable release verification, with normal/-O and actual relocation.

All paths derive from this source file; run from any working directory.
Only Python stdlib and bundled inputs are used. No network or repository.
"""
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
PROGRAMS=(('verify_balance_dynamics.py','normal.json','optimized.json',5784),
          ('seed_balance_window.py','seed_normal.json','seed_optimized.json',271140),
          ('independent_audit.py','audit_normal.json','audit_optimized.json',8776),
          ('cross_observers.py','cross_normal.json','cross_optimized.json',305))
GENERATED=('classification.json','transitions.tsv','seed_balance_window.json',
           'seed_balance_trace.json','cross_observers.json')


def require(condition,detail):
    if not condition:
        raise RuntimeError(detail)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_inputs(root):
    manifest=json.loads((root/'inputs/SOURCE_PINS.json').read_text())
    for entry in manifest['inputs']:
        path=root/'inputs'/entry['local']
        require(path.is_file(),f'missing bundled input {entry["local"]}')
        require(digest(path)==entry['sha256'],f'input hash mismatch {entry["local"]}')
        require(path.stat().st_size==entry['bytes'],f'input size mismatch {entry["local"]}')
    return len(manifest['inputs'])


def run_suite(root,cwd):
    records=[]
    inputs=check_inputs(root)
    for program,normal,optimized,checks in PROGRAMS:
        outputs=[]
        for mode,receipt in (([],normal),(['-O'],optimized)):
            result=subprocess.run([sys.executable,'-B',*mode,str(root/program)],
                                  cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
            require(not result.stderr,f'{program} emitted stderr')
            payload=json.loads(result.stdout)
            require(payload['checks_total']==checks,f'{program} check count changed')
            (root/receipt).write_bytes(result.stdout)
            outputs.append(result.stdout)
        require(outputs[0]==outputs[1],f'{program} normal/-O mismatch')
        records.append({'program':program,'checks':checks,'normal_optimized_identical':True,
                        'stdout_sha256':sha256(outputs[0]).hexdigest()})
    return {'bundled_inputs_checked':inputs,'programs':records,
            'generated_sha256':{name:digest(root/name) for name in GENERATED}}


def main():
    # Even the local suite uses an unrelated working directory.
    with tempfile.TemporaryDirectory(prefix='balance-release-local-cwd-') as directory:
        local=run_suite(ROOT,Path(directory))
    with tempfile.TemporaryDirectory(prefix='balance-release-relocation-') as directory:
        temporary=Path(directory)
        relocated=temporary/'copied-package'
        shutil.copytree(ROOT,relocated,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        unrelated=temporary/'unrelated-working-directory'
        unrelated.mkdir()
        moved=run_suite(relocated,unrelated)
    require(local==moved,'relocated outputs differ from local outputs')
    report={'status':'PASS','scripts':len(PROGRAMS),
            'checks_per_mode':sum(p[3] for p in PROGRAMS),
            'normal_optimized_identical':True,
            'actual_relocation_performed':True,
            'working_directory_separate_from_package':True,
            'relocated_outputs_identical':True,
            'network_or_repository_required':False,
            'verification':local}
    encoded=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (ROOT/'RELEASE_CHECK.json').write_text(encoded)
    print(encoded,end='')


if __name__=='__main__':
    main()
