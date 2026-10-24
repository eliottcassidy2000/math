#!/usr/bin/env python3
"""Actual manifest-only replay in a fresh directory."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parent
manifest = json.loads((root/'MANIFEST.json').read_text())
for record in manifest['files']:
    data = (root/record['path']).read_bytes()
    if len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
        raise ValueError('manifest mismatch: '+record['path'])
with tempfile.TemporaryDirectory(prefix='session50-odd-mass-replay-') as raw:
    target = Path(raw)
    for record in manifest['files']:
        dst = target/record['path']
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root/record['path'], dst)
    commands = [
        [sys.executable,'check.py','--output','replay-normal.json'],
        [sys.executable,'-O','check.py','--output','replay-optimized.json'],
        [sys.executable,'predictor.py','--moments','eight_moments.json','--output','replay-blind.json'],
        [sys.executable,'predictor.py','--moments','twelve_residues.json','--modular-universe-size','2097152','--output','replay-modular.json'],
    ]
    logs = []
    for command in commands:
        result = subprocess.run(command, cwd=target, text=True, capture_output=True, check=True)
        logs.append(result.stdout.strip())
    for actual, expected in [('replay-normal.json','normal.json'),('replay-optimized.json','optimized.json'),('replay-blind.json','blind_prediction.json'),('replay-modular.json','modular_prediction.json')]:
        if (target/actual).read_bytes() != (target/expected).read_bytes():
            raise ValueError('replay result mismatch: '+actual)
print(json.dumps({'status':'PASS','manifest_files':len(manifest['files']),
                  'fresh_directory_manifest_only':True,'normal_optimized_identical':True,
                  'blind_predictor_replay':True,'modular_predictor_replay':True,'logs':logs},indent=2))
