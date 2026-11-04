"""Read the complete pinned bank; build in a temporary directory; never install."""
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
BANK=ROOT/'05-knowledge/results/consecutive_flags_order10_20261028'


def sha(data):return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--limit',type=int)
    args=parser.parse_args()
    if args.limit is not None and not 1<=args.limit<=9733056:raise RuntimeError('bad limit')
    sources={name:sha((HERE/name).read_bytes()) for name in ('run.py','census.cpp')}
    manifest=json.loads((BANK/'bank_manifest.json').read_text(encoding='utf8'))
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('available C++17 compiler g++ required; nothing installed')
    pins={}
    with tempfile.TemporaryDirectory(prefix='two-port-ten-') as tmp:
        tmp=Path(tmp)
        for name in ('parents9.u64','flags10.varint'):
            raw=(BANK/(name+'.gz')).read_bytes();data=gzip.decompress(raw)
            expected=manifest['artifacts'][name+'.gz']
            if sha(raw)!=expected['sha256'] or sha(data)!=expected['uncompressed_sha256']:
                raise RuntimeError('bank hash mismatch: '+name)
            pins[name+'.gz']=sha(raw)
            (tmp/name).write_bytes(data)
        binary=tmp/'census.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'census.cpp'),'-o',str(binary)],check=True,capture_output=True)
        command=[str(binary),str(tmp/'parents9.u64'),str(tmp/'flags10.varint')]
        if args.limit is not None:command.append(str(args.limit))
        result=subprocess.run(command,check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
        receipt=json.loads(result.stdout)
    if sources!={name:sha((HERE/name).read_bytes()) for name in sources}:raise RuntimeError('source changed during run')
    if not args.limit:
        old=json.loads((BANK/'verification10.json').read_text(encoding='utf8'))
        if receipt['sidecar_class_collisions']==receipt['sidecar_automorphism_collisions']==0:
            if receipt['sidecar_automorphism_histogram']!=old['automorphism_order_histogram'] or receipt['sidecar_labelled_mass']!=2**45:
                raise RuntimeError('full bank symmetry/mass disagreement')
    receipt.update(source_sha256=sources,bank_sha256=pins,
                   status='SIDE_CAR_COLLISION' if receipt['collisions'] else 'PASS_INJECTIVE')
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__':main()
