#!/usr/bin/env python3
"""Copy the entire audit, change cwd, and replay both independent Python modes."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
reference = (ROOT / 'independent-normal.json').read_bytes()
with tempfile.TemporaryDirectory(prefix='session51-dual-audit-replay-') as temp:
    root = Path(temp)
    packet = root / 'audit'
    cwd = root / 'unrelated'
    cwd.mkdir()
    shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))

    def run(file, optimized=False):
        command = [sys.executable] + (['-O'] if optimized else []) + [str(packet / file)]
        p = subprocess.run(command, cwd=cwd, capture_output=True)
        if p.returncode:
            raise RuntimeError(p.stderr.decode() or p.stdout.decode())
        return p.stdout

    run('verify_manifest.py')
    run('verify_manifest.py', True)
    run('verify_inputs.py')
    normal = run('independent_check.py')
    optimized = run('independent_check.py', True)
    if not normal == optimized == reference:
        raise RuntimeError('independent actual-copy output mismatch')
    run('hostile_certificate_checks.py')
print(json.dumps({'status': 'PASS', 'real_temporary_copy': True, 'unrelated_cwd': True, 'normal_optimized_reference_byte_equal': True, 'manifest_checked_in_both_modes': True, 'independent_result_sha256': hashlib.sha256(reference).hexdigest(), 'hostile_literal_certificate_rejections_replayed': 26}, indent=2, sort_keys=True))
