"""Reproducible exact order12 scout; prescribed finite universes only."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    sources = {name: digest((HERE / name).read_bytes()) for name in ('scout.cpp', 'run.py')}
    bank = ROOT / '05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw = bank.read_bytes()
    if digest(raw) != 'b4e559ad6cb31980852034f9c7f61e366a7cb479fee4d9a641ad3356745200d3':
        raise RuntimeError('inherited bank hash')
    members = json.loads(gzip.decompress(raw))['members']
    expected = [1, 1, 2, 4, 12, 56, 456, 6880]
    lines = []
    for n, count in enumerate(expected, 1):
        if len(members[str(n)]) != count:
            raise RuntimeError(('bank size', n))
        lines.append(str(count))
        lines.extend(str(row['code']) for row in members[str(n)])
    compiler = shutil.which('g++')
    if not compiler:
        raise RuntimeError('g++ required')
    with tempfile.TemporaryDirectory(prefix='twelve-module-scout-') as tmp:
        tmp = Path(tmp)
        inp, exe = tmp / 'bank.txt', tmp / 'scout.exe'
        inp.write_text('\n'.join(lines) + '\n', encoding='ascii')
        subprocess.run([compiler, '-O3', '-std=c++17', str(HERE / 'scout.cpp'), '-o', str(exe)], check=True, capture_output=True)
        result = json.loads(subprocess.run([str(exe), str(inp)], check=True, stdout=subprocess.PIPE, text=True).stdout)
    if sources != {name: digest((HERE / name).read_bytes()) for name in sources}:
        raise RuntimeError('source changed')
    result.update(source_sha256=sources, bank_sha256=digest(raw), status='FINITE-EXACT scout; universal two-port OPEN')
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
