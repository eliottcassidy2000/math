"""Exact finite certificate; no rank or isogeny-oracle calls."""
from pathlib import Path
import json, hashlib
from itertools import combinations
ROOT=Path(__file__).resolve().parent
INPUT=ROOT.parent/'ranks/phase2/highrank_inputs.json'
def disc(a):
 a1,a2,a3,a4,a6=a
 b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
 b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 return -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
def trace(a,p):
 a1,a2,a3,a4,a6=[x%p for x in a]
 s=0
 for x in range(p):
  d=((a1*x+a3)**2+4*(x**3+a2*x*x+a4*x+a6))%p
  v=pow(d,(p-1)//2,p);s+=(-1 if v==p-1 else v)
 return -s
def main():
 raw=INPUT.read_bytes(); curves=sorted(json.loads(raw)['curves'],key=lambda c:c['id'])
 ids=[c['id'] for c in curves];aa=[list(map(int,c['ainvs'])) for c in curves]
 dd=[disc(a) for a in aa]
 if any(d!=int(c['discriminant']) for d,c in zip(dd,curves)):raise ValueError('discriminant mismatch')
 pending=set(combinations(range(len(ids)),2)); witnesses={};bank={}
 for p in range(3,1001,2):
  if any(p%d==0 for d in range(2,int(p**.5)+1)):continue
  row=[trace(a,p) if d%p else None for a,d in zip(aa,dd)]; bank[str(p)]=row
  for i,j in list(pending):
   if row[i] is not None and row[j] is not None and row[i]!=row[j]:
    witnesses[i,j]=p;pending.remove((i,j))
  if not pending:break
 out={'source_url':json.loads(raw)['url'],'input_sha256':hashlib.sha256(raw).hexdigest(),'created_utc':'2026-10-29T07:06:00Z','ids':ids,'traces':bank,'upper_triangle_first_witness_prime':[[witnesses.get((i,j)) for j in range(i+1,len(ids))] for i in range(len(ids))],'unseparated_pairs':[[ids[i],ids[j]] for i,j in sorted(pending)]}
 (ROOT/'certificate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(json.dumps({'curves':len(ids),'separated':len(witnesses),'unseparated':len(pending),'largest_prime':p,'witness_histogram':{str(q):list(witnesses.values()).count(q) for q in sorted(set(witnesses.values()))}},sort_keys=True))
if __name__=='__main__':main()
