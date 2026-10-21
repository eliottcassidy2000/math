#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def need(ok,message):
    if not ok:raise RuntimeError(message)


def main():
    m=json.loads((HERE/'manifest.json').read_text())
    for name,sha in m['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'frozen file '+name)
    runs=[]
    for script,receipt in [('extension_check.py','extension-check.json'),
                           ('local_ports_check.py','local-ports-check.json'),
                           ('first_owner_check.py','first-owner-check.json'),
                           ('owner_four_check.py','owner-four-check.json')]:
        expected=(HERE/receipt).read_bytes()
        for opt in ([],['-O']):
            got=subprocess.check_output([sys.executable]+opt+[str(HERE/script)],cwd=HERE)
            need(got==expected,'replay '+script+' '+str(opt))
            runs.append(dict(checker=script,optimized=bool(opt),bytes=len(got),
                             sha256=hashlib.sha256(got).hexdigest()))
    print(json.dumps(dict(status='PASS',runs=runs),indent=2,sort_keys=True))


if __name__=='__main__':main()
