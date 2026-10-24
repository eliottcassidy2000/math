#!/usr/bin/env python3
"""Manifest-only copy and replay, with all execution in an unrelated cwd."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import independent as a
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    args = parser.parse_args()
    manifest_hash = a.verify_manifest(HERE)
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    runs = []
    with tempfile.TemporaryDirectory(prefix='s51-score-audit-relocation-') as temp:
        packet = Path(temp)/'packet'
        work = Path(temp)/'unrelated-working-directory'
        packet.mkdir()
        work.mkdir()
        for name in ['MANIFEST.json']+[r['path'] for r in manifest['files']]:
            dest = packet/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(HERE/name,dest)
        def run(script, flags, expected=None):
            result = subprocess.run([sys.executable,*flags,str(packet/script)],cwd=work,capture_output=True)
            a.require(result.returncode == 0,'relocated '+script+': '+result.stderr.decode())
            if expected:
                a.require(result.stdout == (packet/expected).read_bytes(),'relocated exact receipt '+script)
            runs.append(dict(script=script,optimized=bool(flags),result='PASS'))
        for flags in ([],['-O']):
            run('independent.py',flags,'normal.json')
            run('hostiles.py',flags,'hostile-normal.json')
            run('crosscheck.py',flags,'crosscheck-normal.json')
            run('verify_audit.py',flags)
            run('producer/check.py',flags,'producer/normal.json')
            run('producer/verify_manifest.py',flags)
            run('nonrigid/check.py',flags,'nonrigid/normal.json')
            run('nonrigid/verify_manifest.py',flags)
        run('producer/replay.py',[])
        run('nonrigid/replay.py',[])
        integrity = []
        def reject(script, name):
            for flags in ([],['-O']):
                result = subprocess.run([sys.executable,*flags,str(packet/script)],cwd=work,capture_output=True)
                a.require(result.returncode != 0,'integrity mutation accepted: '+name)
                integrity.append(dict(name=name,optimized=bool(flags),result='REJECTED'))
        for name,script in [('normal.json','verify_audit.py'),('producer/check.py','producer/verify_manifest.py'),('nonrigid/normal.json','independent.py')]:
            path = packet/name
            original = path.read_bytes()
            path.write_bytes(original+b' ')
            reject(script,'changed byte count '+name)
            path.write_bytes(original)
        original_manifest = (packet/'MANIFEST.json').read_bytes()
        duplicate = json.loads(original_manifest)
        duplicate['files'].append(duplicate['files'][0])
        (packet/'MANIFEST.json').write_text(json.dumps(duplicate))
        reject('verify_audit.py','duplicate manifest path')
        traversal = json.loads(original_manifest)
        traversal['files'][0]['path'] = '../unlisted-outside-file'
        (packet/'MANIFEST.json').write_text(json.dumps(traversal))
        reject('verify_audit.py','traversing manifest path')
        (packet/'MANIFEST.json').write_bytes(original_manifest)
        run('verify_audit.py',[])
    receipt = dict(status='PASS',manifest_sha256=manifest_hash,
                   manifest_only_physical_copy=True,unrelated_working_directory=True,
                   replayed_steps=runs,independent_semantic_rejections_per_mode=14,
                   integrity_rejections=integrity,
                   producer_and_nonrigid_own_replays_passed=True)
    text = json.dumps(receipt,sort_keys=True,indent=2)+'\n'
    if args.output:
        Path(args.output).write_text(text)
    print(text,end='')
if __name__ == '__main__':
    main()
