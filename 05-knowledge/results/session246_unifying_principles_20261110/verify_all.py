"""Replay Session246 producers and independent controls; stdlib only."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
JOBS = (
    ('descent/verify.py', 'descent/RECEIPT.json'),
    ('bridges/verify_reset.py', 'bridges/RESET_CHECK.json'),
    ('switches/verify.py', 'switches/RESULT.json'),
    ('periodic/verify.py', 'periodic/RECEIPT.json'),
    ('audit/independent.py', 'audit/INDEPENDENT.json'),
    ('audit/periodic_independent.py', 'audit/PERIODIC_INDEPENDENT.json'),
    ('audit/switches_independent.py', 'audit/SWITCHES_INDEPENDENT.json'),
)

def require(test, message):
    if not test:
        raise RuntimeError(message)

def run(path, optimized=False, fail=False):
    cmd = [sys.executable, '-I', '-B']
    if optimized:
        cmd.append('-O')
    cmd.append(str(ROOT / path))
    if fail:
        cmd.append('--fail')
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=180)

if '--fail' in sys.argv:
    require(False, 'intentional aggregate sentinel')

receipts = {}
for program, receipt in JOBS:
    ordinary, optimized = run(program), run(program, optimized=True)
    require(ordinary.returncode == 0, (program, ordinary.stderr))
    require(optimized.returncode == 0, (program, optimized.stderr))
    require(ordinary.stdout == optimized.stdout, (program, 'optimized mismatch'))
    expected = (ROOT / receipt).read_text()
    require(ordinary.stdout == expected, (program, 'stored receipt mismatch'))
    hostile = run(program, optimized=True, fail=True)
    require(hostile.returncode != 0, (program, 'failing sentinel unexpectedly passed'))
    receipts[program] = {
        'sha256': hashlib.sha256(ordinary.stdout.encode()).hexdigest(),
        'normal_optimized': 'MATCH',
        'failure_sentinel': 'FAILS_AS_REQUIRED',
    }

# Audited proof sources are sealed separately from the executable receipts.
first_seals = json.loads((ROOT / 'audit/INDEPENDENT.json').read_text())['source_sha256']
second_seals = json.loads((ROOT / 'audit/SECOND_CHECKPOINT_SHA256.json').read_text())
for relative, expected in {**first_seals, **second_seals}.items():
    require(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected,
            (relative, 'audited source changed'))

for lane in ('descent', 'periodic'):
    for line in (ROOT / lane / 'SHA256SUMS').read_text().splitlines():
        digest, filename = line.split(None, 1)
        require(hashlib.sha256((ROOT / lane / filename).read_bytes()).hexdigest() == digest,
                (lane, filename, 'lane hash mismatch'))

sources = {
    str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(ROOT.rglob('*'))
    if p.is_file() and p.suffix in ('.py', '.md')
}
print(json.dumps({
    'status': 'PASS',
    'proof_scope': 'Structural theorems and scoped physical selectors; Chan-Pak OPEN.',
    'programs': receipts,
    'audited_sources_checked': len(set(first_seals) | set(second_seals)),
    'source_sha256': sources,
}, indent=2, sort_keys=True))
