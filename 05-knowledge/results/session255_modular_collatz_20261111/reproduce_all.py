#!/usr/bin/env python3
"""Replay this packet's new-claim controls, not inherited research censuses."""
import hashlib, json, os, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = ['central_table/check.py', 'collatz_dictionary/check_dictionary.py',
           'finite_infinite/reproduce.py', 'cross_frontiers/probe.py']
results = []
for rel in SCRIPTS:
    started = time.monotonic()
    path = HERE / rel
    p = subprocess.run([sys.executable, '-B', str(path)], cwd=ROOT,
                       text=True, capture_output=True)
    row = {'script': rel, 'exit_code': p.returncode,
           'status': 'PASS' if p.returncode == 0 else 'FAIL',
           'script_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
           'stdout_sha256': hashlib.sha256(p.stdout.encode()).hexdigest(),
           'stdout': p.stdout, 'stderr': p.stderr}
    results.append(row)
    print(f"{row['status']}: {rel}", flush=True)
    if p.returncode:
        print(p.stdout, p.stderr, flush=True)
        break
out = {'scope': 'Only Session255 new constructions', 'engines': results,
       'all_passed': len(results) == len(SCRIPTS) and all(x['exit_code'] == 0 for x in results)}
(HERE/'REPLAY.json').write_text(json.dumps(out, indent=2) + '\n')
raise SystemExit(0 if out['all_passed'] else 1)
