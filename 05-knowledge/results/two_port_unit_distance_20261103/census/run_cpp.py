"""Compile and replay the bounded census; expectations are never overwritten."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-n',type=int,default=9)
    args=parser.parse_args()
    if not 1<=args.max_n<=9:
        raise RuntimeError('orders1..9 only')
    sources={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_cpp.py','exact.cpp')}
    bank=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw=bank.read_bytes()
    data=json.loads(gzip.decompress(raw))['members']
    expected=(1,1,1,2,4,12,56,456,6880,191536)
    compiler=shutil.which('g++')
    if not compiler:
        raise RuntimeError('g++ is required; no software is installed automatically')
    with tempfile.TemporaryDirectory(prefix='two-port-census-') as tmp:
        tmp=Path(tmp)
        binary=tmp/'census.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'exact.cpp'),'-o',str(binary)],check=True,capture_output=True)
        inp=tmp/'bank.txt'
        lines=[str(args.max_n)]
        for n in range(1,args.max_n+1):
            members=data[str(n)]
            if len(members)!=expected[n]:
                raise RuntimeError(('bank count',n,len(members)))
            lines.append(f'{n} {len(members)}')
            lines.extend(f"{r['code']} {r['aut']}" for r in members)
        inp.write_text('\n'.join(lines)+'\n',encoding='ascii')
        result=subprocess.run([str(binary),str(inp)],check=True,capture_output=True,text=True,encoding='utf8')
        receipt=json.loads(result.stdout)
    for name,value in sources.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=value:
            raise RuntimeError('source changed during run: '+name)
    receipt.update(source_sha256=sources,bank_path=bank.relative_to(ROOT).as_posix(),
                   bank_sha256=hashlib.sha256(raw).hexdigest(),
                   status='PASS_INJECTIVE' if not receipt['collision_count'] else 'COLLISION')
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__':main()
