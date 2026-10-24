#!/usr/bin/env python3
"""Audit frozen copies and every inherited Git source pin without importing producers."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
parser = argparse.ArgumentParser()
parser.add_argument('--source-workspace', type=Path, help='Additionally compare original packet directories and exact inherited Git blobs in this session workspace')
args = parser.parse_args()
REPOSITORY = args.source_workspace / 'session50-publication/clean-repo' if args.source_workspace else None


def sha(data):
    return hashlib.sha256(data).hexdigest()


def need(ok, name):
    if not ok:
        raise RuntimeError(name)


receipt = {'status': 'PASS', 'packet_payloads': {}, 'inherited_sources': [], 'source_copy_pins': {}}
for name, original in [('prime-offsets', 'session51-prime-offsets'), ('implication-intervals', 'session51-implication-intervals')]:
    copy = ROOT / 'copies' / name
    manifest = json.loads((copy / 'MANIFEST.json').read_text())
    entries = manifest['files']
    if isinstance(entries, list):
        entries = {r['path']: r['sha256'] for r in entries}
    for filename, expected in entries.items():
        got = sha((copy / filename).read_bytes())
        need(got == expected, 'manifest hash ' + name + '/' + filename)
        if args.source_workspace:
            need((copy / filename).read_bytes() == (args.source_workspace / original / filename).read_bytes(), 'actual-copy equality ' + name + '/' + filename)
    if name == 'implication-intervals':
        actual = {p.relative_to(copy).as_posix() for p in copy.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p != copy / 'MANIFEST.json'}
        need(actual == set(entries), 'strict interval payload inventory')
    receipt['packet_payloads'][name] = {'manifest_sha256': sha((copy / 'MANIFEST.json').read_bytes()), 'declared_payload_count': len(entries), 'declared_payloads_verified': True}
    receipt['source_copy_pins'][name] = {p.relative_to(copy).as_posix(): sha(p.read_bytes()) for p in sorted(copy.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}

prime_sources = json.loads((ROOT / 'copies/prime-offsets/SOURCES.json').read_text())
interval_sources = json.loads((ROOT / 'copies/implication-intervals/SOURCEPINS.json').read_text())
for row in prime_sources['inherited_files']:
    target = prime_sources['repository_commit'] + ':' + row['path']
    if REPOSITORY:
        blob = subprocess.check_output(['git', 'show', target], cwd=REPOSITORY)
        need(sha(blob) == row['sha256'], 'prime inherited Git source ' + target)
    receipt['inherited_sources'].append({'commit_path': target, 'sha256': row['sha256'], 'verified_against_git_blob': bool(REPOSITORY)})
for row in interval_sources['files']:
    target = interval_sources['commit'] + ':' + row['repository_path']
    local = (ROOT / 'copies/implication-intervals' / row['local_path']).read_bytes()
    need(sha(local) == row['sha256'], 'interval inherited local source ' + target)
    if REPOSITORY:
        blob = subprocess.check_output(['git', 'show', target], cwd=REPOSITORY)
        need(blob == local, 'interval inherited Git source ' + target)
    blob_sha = hashlib.sha1(b'blob ' + str(len(local)).encode() + b'\0' + local).hexdigest()
    need(blob_sha == row['git_blob_sha1'], 'Git object hash agrees')
    receipt['inherited_sources'].append({'commit_path': target, 'sha256': row['sha256'], 'git_blob_sha1': blob_sha, 'verified_against_git_blob': bool(REPOSITORY)})
print(json.dumps(receipt, indent=2, sort_keys=True))
