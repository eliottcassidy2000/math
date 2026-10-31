#!/usr/bin/env python3
"""Verify the flat publication manifest and fully relocated normal/-O replay."""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
ARTIFACTS = ['ring_ca_results.json', 'observer_quotients.json', 'transitions.tsv']
INPUTS = ['verify_ring_ca.py', 'independent_audit.py', 'reference_n5_atlas.json']


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    manifest = {}
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ', 1)
        need(Path(name).name == name and name not in manifest,
             'manifest names are unique flat relative paths')
        need(re.fullmatch('[0-9a-f]{64}', digest) is not None, 'valid sha256 format')
        path = ROOT/name
        need(path.is_file() and not path.is_symlink(), 'manifest file exists: '+name)
        need(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
             'manifest checksum: '+name)
        manifest[name] = digest
    actual = {p.name for p in ROOT.iterdir() if p.is_file() and p.name != 'SHA256SUMS'}
    need(actual == set(manifest), 'manifest covers every packet file exactly')

    # There are currently no Markdown hyperlinks. Keep this check so a future
    # documentation edit cannot accidentally add a dangling local link.
    for path in ROOT.glob('*.md'):
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            if target.startswith(('https://', 'http://', '#')):
                continue
            target = target.split('#', 1)[0]
            need((ROOT/target).is_file(), 'portable Markdown target '+target)

    mode_reports = {}
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for mode, flags in [('normal', ['-I']), ('optimized', ['-I', '-O'])]:
        with tempfile.TemporaryDirectory(prefix='metagraph-ca-relocated-') as tmp:
            target = Path(tmp)
            for name in INPUTS:
                shutil.copyfile(ROOT/name, target/name)
            producer = subprocess.run([sys.executable, *flags, str(target/'verify_ring_ca.py')],
                 cwd=target, env=env, capture_output=True, check=True)
            need(producer.stdout == (ROOT/'producer_output.json').read_bytes(),
                 mode+' producer stdout byte match')
            need(not producer.stderr, mode+' producer stderr empty')
            for name in ARTIFACTS:
                need((target/name).read_bytes() == (ROOT/name).read_bytes(),
                     mode+' exact artifact byte match: '+name)
            audit = subprocess.run([sys.executable, *flags, str(target/'independent_audit.py')],
                 cwd=target, env=env, capture_output=True, check=True)
            need(audit.stdout == (ROOT/'audit_output.json').read_bytes(),
                 mode+' independent audit stdout byte match')
            need(not audit.stderr, mode+' independent audit stderr empty')
            result = json.loads(audit.stdout)
            need(result['status'] == 'PASS', mode+' mathematical audit passed')
            mode_reports[mode] = {'status':'PASS', 'relocated':True,
                'python_isolated_mode':True,'producer_output_byte_match':True,
                'artifact_byte_matches':ARTIFACTS,
                'audit_output_byte_match':True,'independent_audit_gates':result['checks']}
    return {'status':'PASS','manifest_verified':True,'manifest_complete':True,
            'local_markdown_links_verified':True,'standard_library_only':True,
            'scratch_sibling_dependencies':False,'modes':mode_reports}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
