#!/usr/bin/env python3
"""Read-only replay of all98 weak-to-strong B states and exact duality."""
import hashlib,json,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent;W=P/'witnesses';checks=0;seen=set();hashes={}
def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
for f in sorted(W.glob('*.json')):
 raw=f.read_bytes();t=json.loads(raw);state=(t['parity'],t['residue'])
 check(t['family']=='B' and state in {(p,d)for p in(0,1)for d in range(24,73)} and state not in seen,'state family or duplicate');seen.add(state)
 if state[0]==1:
  source=json.loads((W/f'B_p0_d{96-state[1]}.json').read_text());words=[]
  for w in source['words'][::-1]:
   z=[]
   for atom in w:
    a=list(atom)
    if a[0]in('block','tail'):a[2]=1-a[2]
    z.append(a)
   words.append(z)
  check(t['words']==words and t.get('split_prefix_blocks',[])==source.get('split_prefix_blocks',[]),'exact dual')
 r=subprocess.run([sys.executable]+(['-O']if sys.flags.optimize else[])+[str(P/'verify_strong_transition.py'),str(f)],capture_output=True,text=True)
 check(r.returncode==0,'transition replay '+f.name+' '+r.stderr);report=json.loads(r.stdout);checks+=report['checks'];hashes[f.name]=hashlib.sha256(raw).hexdigest()
check(seen=={(p,d)for p in(0,1)for d in range(24,73)},'complete table')
print(json.dumps({'status':'passed','states':len(seen),'checks':checks,'complete':True,'certificate_sha256':hashes,'input_hypothesis':'Existing weak B pair; no stronger base or prefix hypothesis','output_prefix':[1,3,6]},indent=2,sort_keys=True))
