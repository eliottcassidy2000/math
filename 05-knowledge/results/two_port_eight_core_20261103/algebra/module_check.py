"""Exact finite module exclusions and controls of the accompanying identities."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    names = ('module_checks.cpp', 'module_check.py', 'scout.cpp')
    pins = {name: sha((HERE / name).read_bytes()) for name in names}
    raw = (ROOT / '05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz').read_bytes()
    if sha(raw) != 'b4e559ad6cb31980852034f9c7f61e366a7cb479fee4d9a641ad3356745200d3':
        raise RuntimeError('bank hash')
    members = json.loads(gzip.decompress(raw))['members']
    lines = []
    for n, count in enumerate([1, 1, 2, 4, 12, 56, 456, 6880], 1):
        if len(members[str(n)]) != count:
            raise RuntimeError('bank count')
        lines.append(str(count))
        lines.extend(str(r['code']) for r in members[str(n)])
    with tempfile.TemporaryDirectory(prefix='twelve-module-audit-') as tmp:
        tmp = Path(tmp)
        inp, exe = tmp / 'bank.txt', tmp / 'check.exe'
        inp.write_text('\n'.join(lines) + '\n', encoding='ascii')
        compiler = shutil.which('g++')
        if not compiler:
            raise RuntimeError('g++ required')
        subprocess.run([compiler, '-O3', '-std=c++17', str(HERE / 'module_checks.cpp'), '-o', str(exe)], check=True, capture_output=True)
        result = json.loads(subprocess.run([str(exe), str(inp)], check=True, capture_output=True, text=True).stdout)
    if pins != {name: sha((HERE / name).read_bytes()) for name in names}:
        raise RuntimeError('source changed')
    result.update(source_sha256=pins, bank_sha256=sha(raw))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
