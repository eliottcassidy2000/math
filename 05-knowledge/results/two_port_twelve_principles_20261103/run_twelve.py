"""Complete declared 12-vertex families, not a complete order-12 census."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import argparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    files = [HERE/'run_twelve.py', HERE/'twelve.cpp']
    sources = {f.name: sha(f.read_bytes()) for f in files}
    bank = ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw = bank.read_bytes()
    members = json.loads(gzip.decompress(raw))['members']['6']
    if len(members) != 56:
        raise RuntimeError('bank count')
    compiler = shutil.which('g++')
    if not compiler:
        raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='twelve-principles-') as temp:
        temp = Path(temp)
        inp, exe = temp/'quotients.txt', temp/'twelve.exe'
        inp.write_text('56\n'+'\n'.join(str(r['code']) for r in members)+'\n', encoding='ascii')
        subprocess.run([compiler, '-O3', '-std=c++17', str(HERE/'twelve.cpp'), '-o', str(exe)],
                       check=True, capture_output=True)
        result = json.loads(subprocess.run([str(exe), str(inp)], check=True,
                                          capture_output=True, text=True).stdout)
    if sources != {f.name: sha(f.read_bytes()) for f in files}:
        raise RuntimeError('source changed')
    result.update(status='FINITE-EXACT declared families only', source_sha256=sources,
                  bank_sha256=sha(raw), universe={
                      '0': 'Every 2048 named one-vertex extension of transitive11',
                      '1': 'Every 2048 named one-vertex extension of cyclic11, steps1..5',
                      '2': 'Every 2048 named one-vertex extension of Paley11, residues1,3,4,5,9',
                      '3': 'Every 56 inherited order6 quotient classes inflated by TT2 at each vertex'},
                  filters='None before counting named presentations. '
                          'Cover skew evaluated only on endpoint-parity survivors; '
                          'cubic tested on every presentation. Families may overlap up to isomorphism.')
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
