#!/usr/bin/env python3
"""Verify this packet's pinned sources, receipts, and exact replays."""
import hashlib
import json
from pathlib import Path
import subprocess


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    base=Path(__file__).resolve().parent
    root=base.parents[2]
    manifest=json.loads((base/'MANIFEST.json').read_text())
    for name,sha in manifest['files'].items():
        if digest((base/name).read_bytes()) != sha:
            raise RuntimeError('manifest mismatch: '+name)
    sources=json.loads((base/'SOURCES.json').read_text())
    for source in sources['source_pins']:
        raw=subprocess.check_output(['git','show',source['commit']+':'+source['path']],cwd=root)
        if digest(raw) != source['sha256']:
            raise RuntimeError('source pin mismatch: '+source['path'])
    for script,label in [('check.py','check'),('audit_monoid.py','monoid'),
                         ('check_prime.py','prime'),('check_native.py','native'),
                         ('check_histogram.py','histogram')]:
        receipt=(base/(label+'.json')).read_bytes()
        if json.loads(receipt)['script_sha256'] != digest((base/script).read_bytes()):
            raise RuntimeError('script/receipt mismatch: '+script)
        for options,suffix in [([],'.json'),(['-O'],'.optimized.json')]:
            replay=subprocess.check_output(['python3','-B',*options,str(base/script)],cwd=root)
            if replay != receipt or replay != (base/(label+suffix)).read_bytes():
                raise RuntimeError('replay mismatch: '+script+' '+suffix)
    print(json.dumps({'status':'PASS','manifest_files':len(manifest['files']),
                      'pinned_repo_sources':len(sources['source_pins']),
                      'normal_optimized_pairs':5},sort_keys=True))


if __name__ == '__main__':
    main()
