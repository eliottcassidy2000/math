#!/usr/bin/env python3
"""Verify literal hashes and rerun using only manifest-listed relocated files."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest-only',action='store_true',required=True)
    parser.add_argument('--max-order',type=int,choices=[7,8],default=7)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    manifest=json.loads((here/'manifest.json').read_text())
    for name,expected in manifest['sha256'].items():
        require(not Path(name).is_absolute() and '..' not in Path(name).parts,'safe relative manifest path')
        require(digest(here/name)==expected,'source hash: '+name)
    for pin in json.loads((here/'source-pins.json').read_text())['sources']:
        require(digest(here/pin['archive_path'])==pin['sha256'],'archive pin')
    suffix='-order8' if args.max_order==8 else ''
    expected_receipt=(here/('receipt'+suffix+'.json')).read_bytes()
    expected_stdout=(here/('order8.out' if args.max_order==8 else 'normal.out')).read_bytes()
    outcomes=[]
    with tempfile.TemporaryDirectory(prefix='marked-growth-manifest-only-') as directory:
        root=Path(directory)
        for name in list(manifest['sha256'])+['manifest.json']:
            target=root/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(here/name,target)
        for optimized in [False,True]:
            output=root/('fresh-optimized.json' if optimized else 'fresh-normal.json')
            command=[sys.executable]+(['-O'] if optimized else [])+[str(root/'transfer.py'),
                     '--max-order',str(args.max_order),'--receipt',str(output)]
            completed=subprocess.run(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
            require(completed.returncode==0,completed.stderr.decode())
            require(completed.stdout==expected_stdout,'relocated stdout differs')
            require(output.read_bytes()==expected_receipt,'relocated receipt differs')
            outcomes.append({'optimized':optimized,'receipt_sha256':digest(output)})
        hostile=subprocess.run([sys.executable,str(root/'refinement_hostile.py')],cwd=root,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
        require(hostile.returncode==0,hostile.stderr.decode())
        require(hostile.stdout==(root/'refinement-hostile.out').read_bytes(),'independent hostile stdout')
    print(json.dumps({'status':'PASS','mode':'manifest-only relocation','max_order':args.max_order,
                      'manifest_sha256':digest(here/'manifest.json'),
                      'files_verified':len(manifest['sha256']),'replays':outcomes},sort_keys=True))

if __name__=='__main__':
    main()
