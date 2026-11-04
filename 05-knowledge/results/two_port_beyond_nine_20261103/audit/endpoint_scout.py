"""Exact inherited-bank scout for the even-order endpoint-balance hypothesis."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    sources = {name: (HERE/name).read_bytes() for name in ('endpoints.cpp','endpoint_scout.py')}
    bank = ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw = bank.read_bytes()
    members = json.loads(gzip.decompress(raw))['members']
    expected = [1,1,2,4,12,56,456,6880,191536]
    compiler = shutil.which('g++')
    if not compiler:
        raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='two-port-endpoint-audit-') as tmp:
        tmp = Path(tmp)
        lines = ['9']
        for n, count in enumerate(expected,1):
            if len(members[str(n)]) != count:
                raise RuntimeError(('bank cardinality',n))
            lines.append(f'{n} {count}')
            lines.extend(str(row['code']) for row in members[str(n)])
        inp = tmp/'bank.txt'
        inp.write_text('\n'.join(lines)+'\n',encoding='ascii')
        exe = tmp/'endpoints.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'endpoints.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result = subprocess.run([str(exe),str(inp)],check=True,capture_output=True,text=True,encoding='utf8')
        receipt = json.loads(result.stdout)
    for name, content in sources.items():
        if (HERE/name).read_bytes() != content:
            raise RuntimeError(('source changed',name))
    receipt.update(status='FINITE-EXACT scout; universal hypothesis unproved',
                   bank_sha256=hashlib.sha256(raw).hexdigest(),
                   source_sha256={name:hashlib.sha256(content).hexdigest() for name,content in sources.items()})
    text = json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text,encoding='utf8',newline='\n')
    else:
        print(text,end='')


if __name__ == '__main__':
    main()
