#!/usr/bin/env python3
"""Closed manifest, source lineage and exact normal/-O replay; standard library."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = {str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file()}
    require(actual == set(manifest['files']) | {'MANIFEST.json'}, 'closed file universe')
    require(not any(p.is_symlink() for p in HERE.rglob('*')), 'no symlink substitution')
    for relative, digest in manifest['files'].items():
        require(sha(HERE / relative) == digest, 'packet hash: ' + relative)
    for relative, digest in manifest['inherited_sources'].items():
        require(sha(ROOT / relative) == digest, 'inherited source hash: ' + relative)
    for child, name in [('aliquot', 'MANIFEST.json'), ('chains', 'SOURCE_MANIFEST.json'),
                        ('tournaments', 'MANIFEST.json'), ('tournament_audit', 'MANIFEST.json')]:
        child_manifest = json.loads((HERE / child / name).read_text())
        for relative, digest in child_manifest['files'].items():
            require(sha(HERE / child / relative) == digest, child + ' manifest: ' + relative)
        for relative, digest in child_manifest.get('inherited', {}).items():
            require(sha(ROOT / relative) == digest, child + ' source: ' + relative)
        for relative, digest in child_manifest.get('producer_files', {}).items():
            require(sha(HERE / 'tournaments' / relative) == digest, 'audited producer: ' + relative)
        if 'independent_audit' in child_manifest:
            audit = child_manifest['independent_audit']
            require(sha(HERE / child / audit['path']) == audit['sha256'], child + ' audit hash')
    corrections = json.loads((HERE / 'CORRECTION_PREIMAGES.json').read_text())
    for entry in corrections['changes']:
        path = ROOT / entry['path']
        require(sha(path) == entry['postimage_sha256'], 'corrected source: ' + entry['path'])
        current = path.read_text()
        require(current.count(entry['replacement']) == 1, 'unique correction replacement')
        restored = current.replace(entry['replacement'], entry['removed'])
        require(hashlib.sha256(restored.encode()).hexdigest() == entry['preimage_sha256'],
                'exact correction preimage: ' + entry['path'])
    suites = [
        ('aliquot/check_aliquot.py', 'aliquot/results.normal.json'),
        ('chains/check.py', 'chains/normal.json'),
        ('chains/audit_aliquot.py', 'chains/aliquot_audit.normal.json'),
        ('tournaments/check_tournaments_124.py', 'tournaments/normal.json'),
        ('build_microcosm.py', 'MICROCOSM.json'),
    ]
    if (HERE / 'tournament_audit/check_independent.py').is_file():
        suites.append(('tournament_audit/check_independent.py',
                       'tournament_audit/independent.normal.json'))
    results = []
    for script, expected in suites:
        frozen = (HERE / expected).read_bytes()
        for optimized in (False, True):
            command = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(HERE / script)]
            run = subprocess.run(command, cwd=HERE, capture_output=True, check=False)
            require(run.returncode == 0, script + ': ' + run.stderr.decode(errors='replace'))
            require(run.stdout == frozen, script + ' frozen receipt mismatch')
        receipt = json.loads(frozen)
        expected_status = ('FINITE-EXACT; all-parameter conclusions use RESULT.md proofs'
                           if script == 'chains/check.py' else 'PASS')
        require(receipt['status'] == expected_status, script + ' status')
        results.append({'script': script, 'normal_and_optimized': 'identical to frozen',
                        'gates': receipt.get('checks', receipt.get('gates'))})
    print(json.dumps({'status': 'PASS', 'manifest_files': len(manifest['files']),
                      'inherited_sources': len(manifest['inherited_sources']),
                      'source_corrections': len(corrections['changes']),
                      'replays': results}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
