#!/usr/bin/env python3
"""New ordinal reservoir operation only; no inherited theorem census."""
import itertools,json,pathlib,random,subprocess
ROOT=pathlib.Path(__file__).resolve().parent

def h(rows):
 n=len(rows)
 if not n:return 1
 dp=[{} for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,c in dp[mask].items():
   rest=rows[v]&~mask
   while rest:
    bit=rest&-rest;rest-=bit;w=bit.bit_length()-1
    dest=dp[mask|bit];dest[w]=dest.get(w,0)+c
 return sum(dp[-1].values())
def induced(rows,vs):
 return [sum(1<<j for j,w in enumerate(vs) if rows[v]>>w&1) for v in vs]
def delete(rows,z):return induced(rows,[v for v in range(len(rows)) if v!=z])
def cut(rows,z):
 vs=[v for v in range(len(rows)) if v!=z];g=[0]*len(vs)
 for i,a in enumerate(vs):
  for j,b in enumerate(vs):
   if rows[z]>>a&1 and rows[b]>>z&1 and rows[a]>>b&1:g[i]|=1<<j;g[j]|=1<<i
 return g
def dom(g):
 out=[]
 for mask in range(1<<len(g)):
  covered=mask
  for v in range(len(g)):
   if mask>>v&1:covered|=g[v]
  if covered==(1<<len(g))-1:out.append(mask)
 return out
def shell(rows,k):
 out=rows[:]
 for _ in range(k):
  n=len(out);out=[r|(7<<n) for r in out]+[1<<(n+1),1<<(n+2),1<<n]
 return out

def check(rows,z,k):
 a=h(delete(rows,z));b=h(rows);d=len(dom(cut(rows,z)));t=shell(rows,k);n=len(rows)
 ht=h(t);ha=h(delete(t,z));gt=cut(t,z)
 assert (ht,ha)==(3**k*b,3**k*a)
 assert gt[:n-1]==cut(rows,z) and not any(gt[n-1:])
 assert len(dom(gt))==d
 assert ht-ha-d+1==3**k*(b-a-d+1)+(3**k-1)*(d-1)
 # Exact Boolean coefficient for every subset of marked shell vertices.
 marks=[n+3*i for i in range(k)]
 for s in range(1<<k):
  coeff=0
  for q in range(1<<k):
   if q&~s:continue
   keep=[v for v in range(len(t)) if v not in [marks[i] for i in range(k) if not(q>>i&1)]]
   u=induced(t,keep);delta=h(u)-h(delete(u,keep.index(z)))
   coeff+=(-1)**((s.bit_count()-q.bit_count()))*delta
  assert coeff==2**s.bit_count()*(b-a)
 return {'n':n,'z':z,'k':k,'H':ht,'deletion_H':ha,'d':d,'slack':ht-ha-d+1}
checks=0;initial=[]
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  rows=[0]*n
  for i,(a,b) in enumerate(pairs):
   if code>>i&1:rows[a]|=1<<b
   else:rows[b]|=1<<a
  for z in range(n):
   for k in range(3):check(rows,z,k);checks+=1
for k in range(5):
 # Direct DP at k=4 is deliberately avoided; finite initial controls k<=2.
 if k<=2:initial.append(check([2,4,1],0,k))
SOURCE_HEAD='f01a7113774301573ebdacc6a49918b3ae54d242'
paths=subprocess.check_output(['git','ls-tree','-r','--name-only',SOURCE_HEAD],text=True).splitlines()
pool=[p for p in paths if p.endswith('.md') and (p.startswith('01-canon/') or p.startswith('05-knowledge/results/'))]
receipt={'status':'FINITE-EXACT','universe':'all labelled rooted tournaments of orders 1..4; k=0,1,2; every shell Boolean subset','operation_checks':checks,'initials':initial,'hostile_empty_cut':check([0],0,2),'source_sampling':{'seed':26211031,'head':SOURCE_HEAD,'pool_size':len(pool),'paths':random.Random(26211031).sample(pool,3)},'claims_not_tested':['Chan-Pak coverage','arbitrary target correction','strong tournament preservation']}
(ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
