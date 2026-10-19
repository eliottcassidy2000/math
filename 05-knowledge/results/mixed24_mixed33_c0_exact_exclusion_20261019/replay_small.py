#!/usr/bin/env python3
"""Fresh standard-library catalogue generation and complete mathematical replay."""
from pathlib import Path
from itertools import combinations
import json,gzip,hashlib,time
import exact_replay_engine as R
ROOT=Path(__file__).resolve().parent

def catalogue(l,u):
 n=2+l+u;edges=list(combinations(range(2,n),2));out=set()
 for mask in range(1<<len(edges)):
  rel={(b,1) for b in range(2,2+l)}|{(0,a) for a in range(2+l,n)}|{e for j,e in enumerate(edges) if mask>>j&1}
  for k in range(n):rel|={(a,b) for a,z in rel if z==k for w,b in rel if w==k}
  out.add(tuple(sum(1<<a for a,b in rel if b==v) for v in range(n)))
 return out

def main():
 start=time.time();summary=[]
 for l,u in [(2,4),(3,3)]:
  label=f'{l}{u}0';path=ROOT/f'full_certificates_{label}.json.gz';records=json.loads(gzip.decompress(path.read_bytes()));cat=catalogue(l,u)
  R.check(len(cat)==4824,'catalogue count');R.check(len(records)==len(cat) and {tuple(r['pred']) for r in records}==cat,'full catalogue coverage')
  for i,r in enumerate(records):
   R.check(r['core_index']==i,'core index');R.replay(tuple(r['pred']),r)
  summary.append({'scope':[l,u,0],'cores':len(cat),'residual_orders':sum(r.get('residual_orders',0) for r in records),'record_sha256':hashlib.sha256(path.read_bytes()).hexdigest()});print(json.dumps(summary[-1]),flush=True)
 out={'status':'FULL_INDEPENDENT_REPLAY_PASSED','scopes':summary,'seconds':time.time()-start};(ROOT/'independent_small_replay.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
if __name__=='__main__':main()
