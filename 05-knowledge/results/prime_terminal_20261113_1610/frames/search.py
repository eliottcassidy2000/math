"""Pinned all-order test of two particular two-run affine profiles."""
from pathlib import Path
import gzip
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'multicore_metagraph_20261113/rooted_growth.json.gz'

def main():
    raw=json.loads(gzip.decompress(SOURCE.read_bytes()))
    bank={}
    for row in raw['records']:
        if row['n']==8:
            bank.setdefault(row['ordinary_code'],row['rows'])
    if len(bank)!=6880:
        raise RuntimeError('Incomplete ordinary order-eight bank')
    data='6880\n'+'\n'.join(' '.join(map(str,r)) for r in bank.values())+'\n'
    exe=HERE/'rooted9.exe'
    subprocess.run(['g++','-O3','-std=c++17',str(HERE/'rooted9.cpp'),'-o',str(exe)],check=True)
    run=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True)
    out=json.loads(run.stdout)
    out['scope']='All order-nine strong frames at H95/113/115/275/277/599/601 via all ordinary order-eight one-vertex extensions; matching<=1 roots only'
    out['source_sha256']=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    (HERE/'search_receipt.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('hits','audit_targets')}))
    print('First hits:',out['hits'][:2])

if __name__=='__main__':main()
