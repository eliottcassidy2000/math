#!/usr/bin/env python3
"""Replay normal/-O from any cwd and an isolated relocated copy; stdlib only."""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
JOBS = [('check_five_families.py','checks.json'),
        ('check_inherited_metagraph.py','metagraph_checks.json'),
        ('check_a5_balance_meeting.py','a5_balance_normal.json')]

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def check_hashes(root):
    path = root/'SHA256SUMS'
    if not path.exists():
        return
    named = []
    for line in path.read_text().splitlines():
        digest,relative = line.split('  ',1)
        need(not Path(relative).is_absolute() and '..' not in Path(relative).parts,
             'relative manifest path')
        need(hashlib.sha256((root/relative).read_bytes()).hexdigest() == digest,
             'SHA256SUMS mismatch '+relative)
        named.append(relative)
    files = sorted(str(p.relative_to(root)) for p in root.rglob('*')
                   if p.is_file() and p.name != 'SHA256SUMS' and '__pycache__' not in p.parts)
    need(sorted(named) == files, 'SHA256SUMS covers every packet file')

def replay(root,cwd):
    check_hashes(root)
    for script,expected in JOBS:
        need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((root/script).read_text()))),
             'optimized-safe explicit gates in '+script)
        expected_bytes = (root/expected).read_bytes()
        for optimize in (False,True):
            command = [sys.executable,'-B']+(['-O'] if optimize else [])+[str(root/script)]
            result = subprocess.run(command,cwd=cwd,check=True,capture_output=True)
            need(result.stdout == expected_bytes, 'byte-exact replay '+script+' optimized='+str(optimize))
    need((root/'a5_balance_normal.json').read_bytes() ==
         (root/'a5_balance_optimized.json').read_bytes(), 'saved A5 receipts agree')
    for entry in json.loads((root/'SOURCE_PINS.json').read_text())['files']:
        need(hashlib.sha256((root/entry['published_path']).read_bytes()).hexdigest() ==
             entry['published_sha256'], 'publication source pin '+entry['published_path'])
        if not entry['modified_for_publication']:
            need(entry['original_sha256'] == entry['published_sha256'], 'unchanged dependency pin')

with tempfile.TemporaryDirectory(prefix='session126_five_families_') as temporary:
    base = Path(temporary)
    other_cwd = base/'unrelated_working_directory'
    other_cwd.mkdir()
    replay(ROOT,other_cwd)
    relocated = base/'relocated_packet'
    shutil.copytree(ROOT,relocated,ignore=shutil.ignore_patterns('__pycache__'))
    replay(relocated,other_cwd)

print(json.dumps({'status':'PASS','scripts':[s for s,_ in JOBS],
                  'normal_optimized_byte_identical':True,
                  'relocated_replay_byte_identical':True,
                  'working_directory_independent':True,
                  'all_test_gates_survive_optimization':True,
                  'source_pin_digests_verified':True,
                  'SHA256SUMS_verified_when_present':True},indent=2,sort_keys=True))
