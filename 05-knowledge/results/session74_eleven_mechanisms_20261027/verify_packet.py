#!/usr/bin/env python3
"""Check a frozen Session 74 packet and replay only its declared current scripts."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
SCHEMA = 'session74-frozen-packet-v1'

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def safe_path(name):
    need(isinstance(name, str), 'manifest path must be text')
    p = PurePosixPath(name)
    need(name == p.as_posix() and not p.is_absolute() and '..' not in p.parts
         and name not in ('', '.') and '\\' not in name, 'unsafe path: ' + repr(name))
    result = ROOT.joinpath(*p.parts)
    for parent in (result, *result.parents):
        if parent == ROOT:
            break
        need(not parent.is_symlink(), 'symlink path: ' + name)
    return result

def select_manifest(name=None):
    if name is None:
        found = sorted(ROOT.glob('FILE_MANIFEST.checkpoint*.json'))
        need(bool(found), 'no frozen checkpoint manifest')
        name = found[-1].name
    need('/' not in name and name.startswith('FILE_MANIFEST.checkpoint') and name.endswith('.json'), 'invalid manifest name')
    return name

def verify_files(manifest_name):
    manifest = json.loads(safe_path(manifest_name).read_text())
    need(manifest['schema'] == SCHEMA, 'manifest schema')
    need(manifest['source_pin'] == '4959b4933cba7e6561b9f19b1d68278ecf6cb2f5', 'research source pin')
    names = [r['path'] for r in manifest['files']]
    need(names == sorted(set(names)) and bool(names), 'unique sorted nonempty file inventory required')
    for row in manifest['files']:
        p = safe_path(row['path'])
        need(p.is_file(), 'missing/nonregular file: ' + row['path'])
        data = p.read_bytes()
        need(len(data) == row['bytes'], 'size mismatch: ' + row['path'])
        need(hashlib.sha256(data).hexdigest() == row['sha256'], 'SHA256 mismatch: ' + row['path'])
        need(stat.S_IMODE(p.stat().st_mode) == row['mode'], 'mode mismatch: ' + row['path'])
    actual = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                    if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts
                    and p != ROOT / manifest_name)
    need(actual == names, 'unlisted or absent packet file')
    for job in manifest.get('replays', []):
        for key in ('script', 'normal', 'optimized'):
            safe_path(job[key])
            need(job[key] in names, 'unlisted replay ' + key)
        need(job['script'].endswith('.py'), 'replay must name Python script')
        need(job['script'] not in manifest.get('inert_source_snapshots', []), 'cannot execute inert source snapshot')
        need(safe_path(job['normal']).read_bytes() == safe_path(job['optimized']).read_bytes(),
             'native normal/optimized mismatch: ' + job['script'])
    for job in manifest.get('negative_controls', []):
        safe_path(job['script'])
        need(job['script'] not in manifest.get('inert_source_snapshots', []), 'cannot execute inert source snapshot')
        need(job['script'] in names and job['exit_code'] != 0 and bool(job['stderr_contains']),
             'invalid negative-control declaration')
    # Native manifests remain byte-identical. Their selected payloads and the
    # exact files reviewed by auditors are additionally bound here.
    for binding in manifest.get('evidence_bindings', []):
        data = safe_path(binding['path']).read_bytes()
        need(binding['path'] in names, 'unlisted evidence binding')
        need(hashlib.sha256(data).hexdigest() == binding['sha256'], 'native/review evidence mismatch')
        if 'equal_to' in binding:
            need(data == safe_path(binding['equal_to']).read_bytes(), 'review/native replay bytes differ')
    if manifest.get('native_copy_ledger'):
        ledger = json.loads(safe_path(manifest['native_copy_ledger']).read_text())
        need(ledger['schema'] == 'session74-native-preserving-curation-v1', 'native copy ledger schema')
        native_names = [r['packet_path'] for r in ledger['files']]
        need(native_names == sorted(set(native_names)), 'native copy paths must be unique and sorted')
        for row in ledger['files']:
            path = row['packet_path']
            need(path == row['native_lane'] + '/' + row['native_path'] and path in names,
                 'native copy mapping differs')
            data = safe_path(path).read_bytes()
            need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                 'native copy ledger mismatch: ' + path)
            need(stat.S_IMODE(safe_path(path).stat().st_mode) == row['mode'], 'native copy mode mismatch')
    return manifest

def replay(manifest, manifest_name):
    results = []
    with tempfile.TemporaryDirectory(prefix='session74-portable-replay-') as temp:
        packet = Path(temp) / 'packet'
        shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))
        env = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'PYTHONHOME')}
        env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
        for job in manifest.get('replays', []):
            expected = safe_path(job['normal']).read_bytes()
            script = packet / job['script']
            for optimized in (False, True):
                output = Path(temp) / ('receipt-' + str(len(results)) + ('-optimized' if optimized else '-normal') + '.json')
                args = [x.replace('{packet}', str(packet)).replace('{output}', str(output)) for x in job.get('args', [])]
                flags = ['-B'] + (['-I'] if job.get('isolated', True) else []) + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(script), *args],
                                      cwd=script.parent, env=env, capture_output=True,
                                      timeout=job.get('timeout_seconds', 600))
                need(proc.returncode == 0, 'replay failed: ' + job['script'] + '\n' + proc.stderr.decode(errors='replace'))
                actual_receipt = output.read_bytes() if job.get('receipt_file', False) else proc.stdout
                need(actual_receipt == expected, 'replay receipt mismatch: ' + job['script'])
            results.append({'script': job['script'], 'normal': 'PASS', 'optimized': 'PASS',
                            'receipt_sha256': hashlib.sha256(expected).hexdigest()})
        for job in manifest.get('negative_controls', []):
            script = packet / job['script']
            for optimized in (False, True):
                flags = ['-B'] + (['-I'] if job.get('isolated', True) else []) + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(script), *job['args']],
                                      cwd=script.parent, env=env, capture_output=True,
                                      timeout=job.get('timeout_seconds', 600))
                need(proc.returncode == job['exit_code'] and not proc.stdout
                     and job['stderr_contains'].encode() in proc.stderr,
                     'negative control missed intended gate: ' + job['script'])
            results.append({'script': job['script'], 'negative_control': 'PASS', 'normal_and_optimized': True})
        for row in manifest['files']:
            p = packet / row['path']
            need(p.is_file() and not p.is_symlink() and hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256'],
                 'replay changed a frozen input/expectation: ' + row['path'])
            need(stat.S_IMODE(p.stat().st_mode) == row['mode'], 'replay changed file mode: ' + row['path'])
        actual = sorted(p.relative_to(packet).as_posix() for p in packet.rglob('*')
                        if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts
                        and p != packet / manifest_name)
        need(actual == [r['path'] for r in manifest['files']], 'replay added or removed packet files')
    return results

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--manifest')
    args = parser.parse_args()
    manifest_name = select_manifest(args.manifest)
    manifest = verify_files(manifest_name)
    print(json.dumps({'status': 'PASS', 'files_verified': len(manifest['files']),
                      'source_pin': manifest['source_pin'],
                      'evidence_bindings_verified': len(manifest.get('evidence_bindings', [])),
                      'manifest': manifest_name,
                      'replays': replay(manifest, manifest_name) if args.replay else 'NOT_RUN'}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
