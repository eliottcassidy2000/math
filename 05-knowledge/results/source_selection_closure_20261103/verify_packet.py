#!/usr/bin/env python3
"""Reproduce five fixed universes, both with and without Python assertions."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = ('grammar_inverse', 'branch_probe', 'native_port_selector',
        'graded_carrier', 'independent_audit')


def main():
    results = []
    for name in JOBS:
        script, receipt = HERE/(name+'.py'), HERE/(name+'.json')
        expected = receipt.read_bytes()
        for optimized in (False, True):
            command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [str(script)]
            actual = subprocess.check_output(command)
            if actual != expected:
                raise RuntimeError(f'{name}: {"optimized" if optimized else "normal"} receipt differs')
        results.append(dict(script=script.name, receipt=receipt.name,
                            script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
                            receipt_sha256=hashlib.sha256(expected).hexdigest(),
                            normal_exact=True, optimized_exact=True))
        print('PASS', name, flush=True)
    (HERE/'verification.json').write_text(json.dumps(dict(status='PASS',
        scope='Five frozen jobs reproduced in normal and optimized Python; no enlarged universe',
        jobs=results), indent=2, sort_keys=True)+'\n')


if __name__ == '__main__': main()
