#!/usr/bin/env python3
"""Replay the five frozen control jobs normally and with optimization."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = [('junction.py', 'junction.json'), ('gluing/check.py', 'gluing/check.json'),
        ('axiom/counter_check.py', 'axiom/counter_check.json'),
        ('machines/check.py', 'machines/normal.json'),
        ('machines/hankel.py', 'machines/hankel.normal.json')]


def digest(data):
    return sha256(data).hexdigest()


def run(job):
    script, receipt = job
    expected = (HERE / receipt).read_bytes()
    for optimization in (False, True):
        command = [sys.executable, '-B'] + (['-O'] if optimization else []) + [str(HERE / script)]
        env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')
        proc = subprocess.run(command, cwd=HERE, env=env, capture_output=True, timeout=180)
        if proc.returncode or proc.stdout != expected:
            raise RuntimeError(f'{script} optimized={optimization}: return={proc.returncode}; '
                               f'stdout_sha={digest(proc.stdout)}; expected_sha={digest(expected)}; '
                               f'stderr={proc.stderr.decode(errors="replace")[:1000]}')
    return dict(script=script, script_sha256=digest((HERE / script).read_bytes()),
                receipt=receipt, receipt_sha256=digest(expected), normal=True, optimized=True)


def main():
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs = list(pool.map(run, JOBS))
    result = dict(status='PASS', jobs=jobs,
                  runner_sha256=digest(Path(__file__).read_bytes()),
                  scope='Five frozen local jobs; byte-identical normal and optimized outputs, not spectrum coverage')
    (HERE / 'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
