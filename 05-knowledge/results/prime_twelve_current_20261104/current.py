import argparse,gzip,json,subprocess,tempfile
from pathlib import Path
from hashlib import sha256
from collections import Counter
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'current.json');args=p.parse_args()
pins={x:sha256((HERE/x).read_bytes()).hexdigest() for x in ['current.cpp','current.py']}
bank=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
data=json.loads(gzip.decompress(bank.read_bytes()))
lines=[]
for n in range(2,10):
 for item in data['members'][str(n)]:lines.append(f"{n} {item['code']}")
with tempfile.TemporaryDirectory(prefix='prime-current-') as temp:
 exe=Path(temp)/'current.exe'
 subprocess.run(['g++','-O3','-std=c++17',str(HERE/'current.cpp'),'-o',str(exe)],check=True)
 result=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
records=[]
counts={}
for line in result.stdout.splitlines():
 n,code,strong,regular,prime,cards=map(int,line.split())
 records.append(dict(n=n,code=code,strong=bool(strong),regular=bool(regular),prime=bool(prime),all_deletion_currents_zero=bool(cards)))
for n in range(2,10):
 r=[x for x in records if x['n']==n]
 counts[n]=dict(classes=len(data['members'][str(n)]),blind=len(r),strong=sum(x['strong'] for x in r),regular=sum(x['regular'] for x in r),prime=sum(x['prime'] for x in r),all_cards=sum(x['all_deletion_currents_zero'] for x in r))
if sum(v['classes'] for v in counts.values())!=198947:raise RuntimeError('bank universe')
if pins!={x:sha256((HERE/x).read_bytes()).hexdigest() for x in pins}:raise RuntimeError('source changed')
receipt=dict(status='FINITE-EXACT complete inherited bank through9; not an order12 result',counts=counts,records=records,bank_sha256=sha256(bank.read_bytes()).hexdigest(),source_sha256=pins)
args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
print(json.dumps(counts,indent=2))

