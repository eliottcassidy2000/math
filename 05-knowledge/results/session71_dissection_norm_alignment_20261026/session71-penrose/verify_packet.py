#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def need(p,m):
    if not p:raise RuntimeError(m)

def run(args):
    return subprocess.run([sys.executable,*args],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)

records=[]
for flags,name in [([], 'normal'),(['-O'],'optimized')]:
    p=run([*flags,'check_carriers.py']);need(p.returncode==0,f'{name} producer failed')
    need(p.stdout==(ROOT/f'{name}.json').read_bytes(),f'{name} stored output differs')
    records.append({'mode':name,'exit_code':p.returncode,'sha256':hashlib.sha256(p.stdout).hexdigest()})
    h=run([*flags,'check_carriers.py','--inject-fault'])
    need(h.returncode==1 and b'injected hostile' in h.stderr,f'{name} inactive hostile gate')
    need(h.stdout==(ROOT/f'hostile.{name}.stdout').read_bytes(),f'{name} hostile stdout differs')
need(records[0]['sha256']==records[1]['sha256'],'optimized output differs')
pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
for p in pins:
    need(hashlib.sha256((ROOT/p['local']).read_bytes()).hexdigest()==p['sha256'],f"source hash: {p['local']}")
print(json.dumps({'status':'PASS','source_snapshots_verified':len(pins),'normal_optimized_match':True,
                  'injected_failure_exits':[1,1],'executions':records},indent=2,sort_keys=True))
