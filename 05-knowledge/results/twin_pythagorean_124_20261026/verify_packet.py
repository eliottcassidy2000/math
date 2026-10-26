#!/usr/bin/env python3
"""Verify a closed mathematical packet and replay all supplied checkers."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def require(test, message):
    if not test:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = {str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file()}
    require(actual == set(manifest['files']) | {'MANIFEST.json'}, 'closed packet file set')
    require(not any(p.is_symlink() for p in HERE.rglob('*')), 'symlink substitution')
    for name, digest in manifest['files'].items():
        require(sha(HERE / name) == digest, 'packet hash: ' + name)
    for name, digest in manifest['sources'].items():
        require(sha(ROOT / name) == digest, 'source hash: ' + name)
    for name, digest in manifest['reviewed_files'].items():
        require(sha(HERE / name) == digest, 'independent review binding: ' + name)
    for folder, filename in manifest['child_manifests']:
        child = json.loads((HERE / folder / filename).read_text())
        for name, digest in child['files'].items():
            require(sha(HERE / folder / name) == digest, 'child manifest: ' + folder + '/' + name)
    field_source = json.loads((HERE / 'finite_field/SOURCE_PINS.json').read_text())
    incoming = field_source['incoming_input']
    blob = (ROOT / incoming['path']).read_bytes()
    git_blob = hashlib.sha1(b'blob ' + str(len(blob)).encode() + b'\0' + blob).hexdigest()
    require(git_blob == incoming['github_blob_sha'], 'incoming field-source blob')
    concurrent = field_source['concurrent_input']
    require(sha(ROOT / concurrent['future_path']) == concurrent['sha256'], 'field/tree shared convention')
    result = []
    for suite in manifest['suites']:
        script, receipt = suite['script'], suite['receipt']
        frozen = (HERE / receipt).read_bytes()
        decoded = json.loads(frozen)
        require(decoded['status'] == suite['status'], 'declared receipt status: ' + script)
        for optimized in (False, True):
            command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [str(HERE / script)]
            run = subprocess.run(command, cwd=HERE, capture_output=True, check=False)
            require(run.returncode == 0, script + ': ' + run.stderr.decode(errors='replace'))
            require(run.stdout == frozen, script + ': replay differs from frozen receipt')
        result.append({'script': script, 'normal_and_optimized': 'identical to frozen'})
    print(json.dumps({'status': 'PASS', 'files': len(manifest['files']),
        'sources': len(manifest['sources']), 'review_bindings': len(manifest['reviewed_files']),
        'producer_gates': manifest['producer_gates'],
        'independent_executable_gates': manifest['independent_executable_gates'],
        'replays': result}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
