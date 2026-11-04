"""Replay frozen receipts and check artifacts; never regenerate expectations."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINES = (
    ('profiles/check.py', 'profiles/receipt.json'),
    ('rank_audit/check.py', 'rank_audit/receipt.json'),
    ('continuation/check.py', 'continuation/normal.json'),
    ('continuation_audit/check.py', 'continuation_audit/receipt.json'),
    ('arithmetic/verify.py', 'arithmetic/receipt.json'),
    ('rank_audit/arithmetic_check.py', 'rank_audit/arithmetic_receipt.json'),
)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text(encoding='utf8'))
    for rel, expected in manifest['files'].items():
        require(sha(HERE / rel) == expected, 'artifact hash mismatch: ' + rel)
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p.name != 'MANIFEST.json'
              and p.suffix != '.pyc'}
    require(actual == set(manifest['files']), 'artifact inventory mismatch')
    sources = json.loads((HERE / 'SOURCE_PINS.json').read_text(encoding='utf8'))
    for rel, expected in sources['files'].items():
        require(sha(ROOT / rel) == expected, 'inherited source changed; review it: ' + rel)
    return len(manifest['files']), len(sources['files'])


def main():
    files, sources = verify_pins()
    require(sys.argv[1:] in ([], ['--pins-only']), 'usage: verify.py [--pins-only]')
    if sys.argv[1:] == ['--pins-only']:
        print(json.dumps({'status': 'PASS', 'artifact_files': files,
                          'inherited_source_files': sources, 'mode': 'pins-only'}, sort_keys=True))
        return
    # Freeze ALL expectations before starting ANY producer.
    expected = {script: (HERE / receipt).read_text(encoding='utf8') for script, receipt in ENGINES}
    for text in expected.values():
        json.loads(text)
    runs = []
    for script, receipt in ENGINES:
        path = HERE / script
        for optimized in (False, True):
            command = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(path)]
            result = subprocess.run(command, cwd=path.parent, capture_output=True,
                                    text=True, encoding='utf8', check=True)
            # Text mode normalizes platform line endings, not JSON content.
            require(result.stdout == expected[script], 'frozen receipt mismatch: ' + script)
            require(not result.stderr, 'unexpected stderr: ' + script)
            runs.append([script, 'optimized' if optimized else 'normal'])
            print('PASS ' + script + (' -O' if optimized else ''), flush=True)
    require(verify_pins() == (files, sources), 'pins changed during replay')
    print(json.dumps({'status': 'PASS', 'artifact_files': files, 'inherited_source_files': sources,
                      'engine_runs': len(runs), 'expected_receipts_modified': False}, sort_keys=True))


if __name__ == '__main__':
    main()
