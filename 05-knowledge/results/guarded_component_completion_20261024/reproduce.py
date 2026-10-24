"""Replay every frozen verifier, propagating optimization mode explicitly."""
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = (
    'exchange/check.py',
    'native/check.py',
    'native/audit.py',
    'sections/pair_sections.py',
    'sections/native_2125.py',
    'sections/check.py',
    'native/xor_tree.py',
    'sections/phase_sections.py',
    'sections/admission.py',
    'sections/phase_route.py',
    'theory/check.py',
    'theory/check_reset.py',
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
print('PASS complete guarded-component packet')
