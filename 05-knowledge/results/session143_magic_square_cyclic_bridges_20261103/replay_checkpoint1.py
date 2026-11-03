#!/usr/bin/env python3
"""Portable Session143 checkpoint replay and fail-closed corruption controls."""
import argparse, hashlib, json, shutil, subprocess, sys, tempfile
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def manifest_rows(root, number):
    return json.loads((root / ('SOURCE_FILES_CHECKPOINT' + number + '.json')).read_text())['files']


def integrity(root, number):
    seen = set()
    rows = manifest_rows(root, number)
    for row in rows:
        p = Path(row['path'])
        need(str(p) == row['path'] and not p.is_absolute() and '..' not in p.parts and row['path'], 'unsafe manifest path')
        need(row['path'] not in seen, 'duplicate manifest path')
        seen.add(row['path'])
        f = root / p
        need(f.is_file() and not f.is_symlink(), 'missing source ' + row['path'])
        data = f.read_bytes()
        need(len(data) == row['bytes'], 'source size ' + row['path'])
        need(sha(data) == row['sha256'], 'source pin ' + row['path'])
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-sources-only', type=Path)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    number = Path(__file__).stem.removeprefix('replay_checkpoint')
    config_name = 'REPLAY_CONFIG_CHECKPOINT' + number + '.json'
    if args.verify_sources_only:
        print(json.dumps({'status': 'PASS', 'source_files': integrity(args.verify_sources_only, number)}, sort_keys=True))
        return
    count = integrity(source, number)
    config = json.loads((source / config_name).read_text())
    records = []
    with tempfile.TemporaryDirectory(prefix='session143-checkpoint' + number + '-') as td:
        root = Path(td) / 'relocated'
        root.mkdir()
        cwd = Path(td) / 'unrelated-cwd'
        cwd.mkdir()
        for row in manifest_rows(source, number):
            dest = root / row['path']
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source / row['path'], dest)
        for name in [Path(__file__).name, config_name, 'SOURCE_FILES_CHECKPOINT' + number + '.json']:
            shutil.copyfile(source / name, root / name)
        def run(test, mode):
            command = [sys.executable, '-B'] + (['-O'] if mode == 'optimized' else [])
            command += [str(root / test['script']), *test.get('argv', [])]
            cp = subprocess.run(command, cwd=cwd, capture_output=True, timeout=360)
            expected = test.get('returncode', 0)
            need(cp.returncode == expected, test['label'] + ' returncode ' + str(cp.returncode) + ' ' + cp.stderr.decode()[-1600:])
            if test.get('output'):
                need(cp.stdout == (source / test['output']).read_bytes(), test['label'] + ' frozen output mismatch')
            if test.get('reason'):
                need(test['reason'].encode() in cp.stdout + cp.stderr, test['label'] + ' rejection reason absent')
            records.append({'label': test['label'], 'mode': mode, 'returncode': cp.returncode, 'stdout_sha256': sha(cp.stdout), 'expected_rejection': expected != 0})
        for mode in ('normal', 'optimized'):
            for test in config['tests']:
                run(test, mode)
            run({'label': 'source-integrity', 'script': Path(__file__).name, 'argv': ['--verify-sources-only', str(root)]}, mode)
        for mutation in config.get('mutations', []):
            f = root / mutation['target']
            original = f.read_text()
            need(original.count(mutation['before']) == 1, 'mutation target ' + mutation['label'])
            f.write_text(original.replace(mutation['before'], mutation['after']))
            for mode in ('normal', 'optimized'):
                run(mutation, mode)
            f.write_text(original)
        f = root / config['proof_corruption']
        original = f.read_bytes()
        need(bool(original), 'empty corruption target')
        f.write_bytes(bytes([original[0] ^ 1]) + original[1:])
        for mode in ('normal', 'optimized'):
            run({'label': 'same-size-proof-corruption', 'script': Path(__file__).name, 'argv': ['--verify-sources-only', str(root)], 'returncode': 1, 'reason': 'source pin ' + config['proof_corruption']}, mode)
        f.write_bytes(original)
        need(integrity(root, number) == count, 'restored source integrity')
    report = {'status': 'PASS', 'checkpoint': number, 'fresh_relocation': True, 'unrelated_cwd': True, 'normal_and_optimized': True, 'source_files_verified': count, 'total_executions': len(records), 'positive_executions': sum(not r['expected_rejection'] for r in records), 'expected_rejections': sum(r['expected_rejection'] for r in records), 'source_manifest_sha256': sha((source / ('SOURCE_FILES_CHECKPOINT' + number + '.json')).read_bytes()), 'scope': config['scope'], 'runs': records}
    result = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(result)
    print(result, end='')


if __name__ == '__main__':
    main()
