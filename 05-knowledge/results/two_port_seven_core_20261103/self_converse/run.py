"""Reproduce the full self-converse family census from small complete cores."""
import argparse,hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(body):return hashlib.sha256(body).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--max-half',type=int,default=5);args=p.parse_args()
    if not 1<=args.max_half<=5:raise RuntimeError('half-orders1..5')
    sources={n:sha((HERE/n).read_bytes()) for n in ('run.py','census.cpp')}
    bank=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives.json'
    raw=bank.read_bytes();members=json.loads(raw)['members']
    expected=(0,1,1,2,4,12)
    lines=[str(args.max_half)]
    for h in range(1,args.max_half+1):
        if len(members[str(h)])!=expected[h]:raise RuntimeError('core bank count')
        lines.append(str(expected[h]));lines.extend(f"{r['code']} {r['aut']}" for r in members[str(h)])
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='two-port-sc11-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';exe=tmp/'census.exe'
        inp.write_text('\n'.join(lines)+'\n',encoding='ascii')
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'census.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result=subprocess.run([str(exe),str(inp)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
        receipt=json.loads(result.stdout)
    if sources!={n:sha((HERE/n).read_bytes()) for n in sources}:raise RuntimeError('sources changed')
    if [r['classes'] for r in receipt['orders']]!=[2,8,88,2752,279968][:args.max_half]:raise RuntimeError('recovered class-count control')
    receipt.update(source_sha256=sources,bank_sha256=sha(raw),status='COLLISION' if any(r['collisions'] for r in receipt['orders']) else 'PASS_WITHIN_SELF_CONVERSE_FAMILY')
    print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
