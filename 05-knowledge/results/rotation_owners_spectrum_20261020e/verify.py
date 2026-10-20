"""Authenticate this research packet and optionally replay its exact controls."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = ('black_check.py', 'owners_check.py', 'regular_check.py',
           'spectrum_check.py', 'arithmetic_check.py', 'owners_independent.py')


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def normalized(data):
    return data.replace(b'\r\n', b'\n')


def digest(path):
    data = path.read_bytes()
    if path.suffix != '.gz':
        data = normalized(data)
    return hashlib.sha256(data).hexdigest()


def inputs():
    own = [p for p in HERE.iterdir() if p.is_file()
           and p.suffix in ('.py', '.md', '.out')
           and p.name != 'verify.out']
    deps = [ROOT / p for p in (
        '05-knowledge/results/deficiency_three_20261020/ownership_check.py',
        '05-knowledge/results/defect_boundary_20261020/sparse_count.py',
        '05-knowledge/results/metagraph_growth_20261020/enumeration.py',
        '05-knowledge/results/metagraph_growth_20261020/growth_check.py',
        '05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz',
    )]
    return sorted(set(own + deps))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = HERE / 'manifest.json'
    if args.record:
        manifest.write_text(json.dumps({
            'normalization': 'CRLF to LF for text; binary gzip unchanged',
            'files': {p.relative_to(ROOT).as_posix(): digest(p) for p in inputs()}
        }, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    data = json.loads(manifest.read_text(encoding='utf-8'))
    need(set(data['files']) == {p.relative_to(ROOT).as_posix() for p in inputs()},
         'manifest coverage changed')
    for name, expected in data['files'].items():
        need(digest(ROOT / name) == expected, f'input changed: {name}')
    print('SHA256 inputs verified:', len(data['files']), flush=True)
    for name in SCRIPTS:
        tree = ast.parse((HERE / name).read_text(encoding='utf-8'))
        need(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
             f'removable assertion in {name}')
        need((HERE / name.replace('.py', '.out')).is_file(), 'missing receipt')
    print('All six sources parse; no removable assertions.', flush=True)
    if args.replay:
        for name in SCRIPTS:
            expected = normalized((HERE / name.replace('.py', '.out')).read_bytes())
            for flags in ([], ['-O']):
                result = subprocess.run([sys.executable, *flags, str(HERE / name)],
                                        cwd=ROOT, capture_output=True, check=False)
                need(result.returncode == 0,
                     f'{name} {flags}: {result.stderr.decode(errors="replace")}')
                need(normalized(result.stdout) == expected,
                     f'transcript mismatch: {name} {flags}')
            print('Normal and optimized transcript equality:', name, flush=True)
    print('PASS', flush=True)


if __name__ == '__main__':
    main()
