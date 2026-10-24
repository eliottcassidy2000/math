"""Replay frozen certificates; optimization never disables a required check."""
from pathlib import Path
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = (
    'coordinates.py',
    'observer.py',
    'coverage/check.py',
    'erasure/check.py',
    'erasure/audit.py',
    'native/route.py',
    'theory/check.py',
    'theory/audit.py',
    'theory/repair_labels.py',
    'manifest.py',
)

for script in SCRIPTS:
    command = [sys.executable] + (['-O'] if sys.flags.optimize else []) + ['-B', str(HERE/script)]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(f'FAILED {script}: {result.returncode}')
    print(json.dumps({'script': script, 'optimized': bool(sys.flags.optimize),
                      'output': result.stdout.strip()}), flush=True)
print('PASS complete external-pin packet')
