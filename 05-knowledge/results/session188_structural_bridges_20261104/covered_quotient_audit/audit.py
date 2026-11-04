"""Literal-permutation independent audit; imports no producer helpers."""
from itertools import combinations,permutations,product
from functools import lru_cache
from pathlib import Path
import hashlib,json,subprocess,sys,math
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'covered_quotient'
def need(ok,message):
 if not ok:raise RuntimeError(message)
def rows(n,code):
 R=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  R[u if code>>k&1 else v]|=1<<(v if code>>k&1 else u)
 return tuple(R)
def ext(R,a,c):
 n=len(R);U=list(R)+[1<<(n+1),0]
 for u in range(n):
  for bit,new in ((a,n),(c,n+1)):
   if bit>>u&1:U[u]|=1<<new
   else:U[new]|=1<<u
 return tuple(U)
@lru_cache(None)
def literal(R):
 if not R:return (1,())
 ss=[0]*len(R)
 for p in permutations(range(len(R))):
  if all(R[u]>>v&1 for u,v in zip(p,p[1:])):ss[p[0]]+=1
 return sum(ss),tuple(ss)
def induced(R,S):
 ns=[i for i in range(len(R)) if S>>i&1]
 return tuple(sum(1<<j for j,v in enumerate(ns) if R[u]>>v&1) for u in ns)
def profile(R):
 n=len(R);q=[0]*(n+1)
 # Assign each name an ordered path slot, then multiply literal Hamilton counts.
 for k in range(1,n+1):
  for assignment in product(range(k),repeat=n):
   masks=[sum(1<<u for u in range(n) if assignment[u]==j) for j in range(k)]
   if all(masks):q[k]+=math.prod(literal(induced(R,M))[0] for M in masks)
 return q
def cycle(blocks,back=False):
 offsets=[0]
 for B in blocks:offsets.append(offsets[-1]+len(B))
 R=[0]*offsets[-1]
 for i,B in enumerate(blocks):
  j=(i+(-1 if back else 1))%3
  for u,row in enumerate(B):R[offsets[i]+u]=(row<<offsets[i])|(((1<<len(blocks[j]))-1)<<offsets[j])
 return tuple(R),offsets

def main():
 pins={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in ('PROOF.md','verify.py')}
 tomography=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   R=rows(n,code);H,ss=literal(R)
   for u,s in enumerate(ss):
    need(literal(ext(R,1<<u,0))[0]-2*literal(ext(R,0,1<<u))[0]+H==2*s,'literal start tomography')
    tomography+=1
 children=[(0,),(2,0),(2,4,1),(6,4,0)];qs={R:profile(R) for R in children};run_checks=0
 for blocks in product(children,repeat=3):
  if sum(map(len,blocks))>7:continue
  R,offsets=cycle(blocks);U,_=cycle(blocks,True);ss=literal(R)[1];uu=literal(U)[1];ps=[qs[B]+[0]*4 for B in blocks]
  for i,B in enumerate(blocks):
   j=(i+1)%3;l=(i+2)%3
   expected=sum(ps[i][k]*(ps[j][k]*ps[l][k-1]-ps[j][k-1]*ps[l][k]) for k in range(1,4))
   actual=sum(ss[offsets[i]:offsets[i+1]])-sum(uu[offsets[i]:offsets[i+1]])
   need(actual==expected,'literal cyclic-run formula');run_checks+=1
 # Literal unequal minimum case and arbitrary-child hostile.
 blocks=((2,0),(2,0),(6,4,0));R,offsets=cycle(blocks);U,_=cycle(blocks,True)
 need(sum(literal(R)[1][:2])!=sum(literal(U)[1][:2]),'literal strict unequal skew')
 bad=(8,1,3,6);p=profile(bad);d=profile(induced(bad,14));need(p[2]*d[1]-p[1]*d[2]==-8,'literal deletion hostile')
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized,'producer normal optimized agreement')
 for name,data in [('producer.normal.stdout',normal),('producer.optimized.stdout',optimized)]: (HERE/name).write_bytes(data)
 need(pins=={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in pins},'source stayed frozen')
 producer=json.loads(normal)
 summary={k:producer[k] for k in ('status','tomography_checks','stirling_minors_through_order_30','arbitrary_child_run_checks')};summary['literal_unequal_batteries']=len(producer['literal_unequal_batteries'])
 result={'status':'PASS','independent_tomography_checks':tomography,'independent_cyclic_run_checks':run_checks,'independent_method':'Literal permutations for all Hamilton counts; explicit vertex-to-slot assignments for child covers; no producer imports.','producer_summary':summary,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
