#!/usr/bin/env python3
"""Literal-hash verification and actual manifest-only isolated replays."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_checked(command, directory, expected_stdout):
    completed = subprocess.run(command, cwd=directory, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=False)
    require(completed.returncode == 0, completed.stderr.decode())
    require(completed.stdout == expected_stdout, 'isolated stdout differs: ' + str(command))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest-only', action='store_true', required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    manifest = json.loads((here / 'manifest.json').read_text())
    for name, expected in manifest['sha256'].items():
        require(not Path(name).is_absolute() and '..' not in Path(name).parts, 'safe manifest path')
        require(digest(here / name) == expected, 'manifest hash: ' + name)
    for record in json.loads((here / 'source-pins.json').read_text())['sources']:
        require(digest(here / record['archive_path']) == record['sha256'], 'literal archive pin')
    outcomes = []
    with tempfile.TemporaryDirectory(prefix='blind-refinement-manifest-') as directory:
        root = Path(directory)
        for name in list(manifest['sha256']) + ['manifest.json']:
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(here / name, target)
        for optimized in [False, True]:
            prefix = [sys.executable] + (['-O'] if optimized else [])
            fresh = root / 'fresh-receipt.json'
            run_checked(prefix + [str(root / 'predictor.py'), '--receipt', str(fresh)], root,
                        (root / 'normal.out').read_bytes())
            require(fresh.read_bytes() == (root / 'receipt.json').read_bytes(), 'isolated receipt differs')
            for script, expected in [('validate_inputs.py', 'input-validation.out'),
                                     ('hostile_controls.py', 'hostile-controls.out')]:
                run_checked(prefix + [str(root / script)], root, (root / expected).read_bytes())
            generated = root / 'fresh-atlas.json'
            run_checked(prefix + [str(root / 'freeze_inputs.py'), '--output', str(generated)], root,
                        (root / 'input-generation.out').read_bytes())
            require(generated.read_bytes() == (root / 'parents-through-seven.json').read_bytes(),
                    'isolated input regeneration differs')
            outcomes.append({'optimized': optimized, 'receipt_sha256': digest(fresh),
                             'regenerated_atlas_sha256': digest(generated)})
        bare = root / 'minimal-blind-inputs'
        bare.mkdir()
        for name in ['predictor.py', 'parents-through-seven.json']:
            shutil.copyfile(root / name, bare / name)
        for optimized in [False, True]:
            prefix = [sys.executable] + (['-O'] if optimized else [])
            fresh = bare / 'fresh.json'
            run_checked(prefix + [str(bare / 'predictor.py'), '--receipt', str(fresh)], bare,
                        (root / 'normal.out').read_bytes())
            require(fresh.read_bytes() == (root / 'receipt.json').read_bytes(), 'minimal blind receipt differs')
    print(json.dumps({'status': 'PASS', 'mode': 'actual manifest-only and two-file blind relocation',
                      'manifest_sha256': digest(here / 'manifest.json'),
                      'files_verified': len(manifest['sha256']), 'replays': outcomes,
                      'minimal_blind_files': ['predictor.py', 'parents-through-seven.json']}, sort_keys=True))


if __name__ == '__main__':
    main()
