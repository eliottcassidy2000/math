#!/usr/bin/env python3
"""Actual relocation: copy only manifest bytes into a clean temp directory.
Run independent audit, frozen producer, and comparison under normal and -O.
No source checkout, sibling directory, network, or installed package is used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from verify_manifest import verify

H=Path(__file__).resolve().parent

def require(ok,why):
    if not ok:
        raise ValueError(why)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    files_verified=verify(H)
    m=json.loads((H/'MANIFEST.json').read_text())
    runs=[]
    env=dict(os.environ)
    env.pop('PYTHONPATH',None)
    env['PYTHONDONTWRITEBYTECODE']='1'
    with tempfile.TemporaryDirectory(prefix='marked-regularization-audit-') as tmp:
        root=Path(tmp)/'portable'
        root.mkdir()
        for name in [r['path'] for r in m['files']]+['MANIFEST.json']:
            dst=root/name
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(H/name,dst)
        require(verify(root)==files_verified,'relocated digest verification')
        for options,label in [([], 'normal'),(['-O'], 'optimized')]:
            for script,output,reference in [
                ('audit.py','fresh-audit-'+label+'.json','normal.json'),
                ('producer/check.py','fresh-producer-'+label+'.json','producer/normal.json'),
                ('compare.py','fresh-compare-'+label+'.json','comparison-normal.json')]:
                cmd=[sys.executable,*options,str(root/script),'--output',str(root/output)]
                done=subprocess.run(cmd,cwd=tmp,env=env,capture_output=True,check=True)
                data=(root/output).read_bytes()
                require(data == (root/reference).read_bytes(),'relocated receipt differs: '+script+' '+label)
                require(done.stdout == data,'stdout differs from receipt')
                require(not done.stderr,'unexpected replay stderr')
                runs.append({'script':script,'mode':label,'status':'PASS','sha256':hashlib.sha256(data).hexdigest()})
        # Integrity negative control against the actual relocated packet.
        corrupt=root/'producer/normal.json'
        corrupt.write_bytes(corrupt.read_bytes()+b' ')
        try:
            verify(root)
        except ValueError:
            corruption_rejected=True
        else:
            raise ValueError('manifest accepted a corrupted producer receipt')
    result={'status':'PASS','copied_only_manifest_entries_plus_manifest':True,
            'files_verified_before_relocation':files_verified,'runs':runs,
            'corrupted_relocated_producer_receipt_rejected':corruption_rejected,
            'producer_manifest_sha256':hashlib.sha256((H/'producer/MANIFEST.json').read_bytes()).hexdigest()}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    Path(args.output).write_text(payload)
    print(payload,end='')

if __name__ == '__main__':
    main()
