"""Recompile and reproduce the incoming independent complete k17 audit."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'target_directed_cores_20261112/ports/complete_audit.cpp'

with tempfile.TemporaryDirectory(prefix='ready_hole_replay_') as temporary:
    executable = Path(temporary) / 'audit.exe'
    output = Path(temporary) / 'result.json'
    compiler = shutil.which('g++') or shutil.which('c++')
    if compiler is None:
        raise RuntimeError('A C++17 compiler is required')
    subprocess.run([compiler, '-O3', '-std=c++17', str(SOURCE), '-o', str(executable)], check=True)
    subprocess.run([str(executable), '53084231', str(output)], check=True)
    actual = json.loads(output.read_text(encoding='utf-8'))
    expected = json.loads((SOURCE.parent / 'complete_audit_53084231.json').read_text(encoding='utf-8'))
    if actual != expected:
        raise RuntimeError('Incoming all-order receipt mismatch')
    report = dict(status='PASS', target=53084231, incoming_audit_reproduced_exactly=True,
                  total_cores=actual['total_cores'], direct_outer_tests=actual['direct_outer_tests'],
                  found=actual['found'], source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  scope='Independent incoming depth-first/direct-word implementation recompiled and replayed; its complete reduction and two implementations were reviewed separately.')
    suffix = '' if __debug__ else '.optimized'
    (HERE / ('incoming_hole_replay' + suffix + '.json')).write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report))
