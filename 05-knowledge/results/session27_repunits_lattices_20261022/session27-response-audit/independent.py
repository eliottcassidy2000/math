"""Independent permutation controls; no imports of response or DP implementations."""
import itertools, json, hashlib
from pathlib import Path
from fractions import Fraction
OUT=Path(__file__).parent
ROOT=Path(__file__).resolve().parent.parent
def H(q):
 n=len(q)
 return sum(all(q[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
def extend(q,w):
 k=len(q); n=k+len(w); r=[a[:]+[0]*len(w) for a in q]+[[0]*n for _ in w]
 for i,a in enumerate(w):
  for u in range(k):r[u][k+i]=(a>>u)&1;r[k+i][u]=1-r[u][k+i]
  for j in range(i+1,len(w)):r[k+i][k+j]=1
 return r
def core(k,m):
 q=[[0]*k for _ in range(k)]
 for i,(u,v) in enumerate(itertools.combinations(range(k),2)):q[u][v]=(m>>i)&1;q[v][u]=1-q[u][v]
 return q
def direct(q):
 k=len(q);s=[0]*k;e=[0]*k;c=[[0]*k for _ in q]
 # Different cut construction: choose first path vertex subset, independently enumerate each path.
 for p in itertools.permutations(range(k)):
  if all(q[p[i]][p[i+1]] for i in range(k-1)):s[p[0]]+=1;e[p[-1]]+=1
 for mask in range(1,(1<<k)-1):
  a=[i for i in range(k) if mask>>i&1]; b=[i for i in range(k) if not mask>>i&1]
  left=[p for p in itertools.permutations(a) if all(q[p[i]][p[i+1]] for i in range(len(p)-1))]
  right=[p for p in itertools.permutations(b) if all(q[p[i]][p[i+1]] for i in range(len(p)-1))]
  for p in left:
   for r in right:c[p[-1]][r[0]]+=1
 return sum(s),s,e,c
checks=0; batteries=0;flipchecks=0
for k in range(1,5):
 for m in range(1<<(k*(k-1)//2)):
  q=core(k,m);h,s,e,c=direct(q);l=[e[u]-s[u]+sum(c[u]) for u in range(k)];p=[s[u]+e[u]+sum(c[u]) for u in range(k)]
  def f(w):return H(extend(q,w))
  single=[f((1<<u,)) for u in range(k)]
  for a in range(1<<k):
   expected=h+sum(l[u] for u in range(k) if a>>u&1)-sum((c[u][v]+c[v][u]) for u,v in itertools.combinations(range(k),2) if a>>u&1 and a>>v&1)
   if f((a,))!=expected:raise RuntimeError('one port')
   checks+=1
  pr=[f((1<<u,0))-single[u] for u in range(k)]
  if pr!=p:raise RuntimeError('potential')
  for u,v in itertools.combinations(range(k),2):
   d=f((1<<u,1<<v))-f((1<<v,1<<u));j=pr[u]-pr[v]-d
   b=h+(single[u]-h)+(single[v]-h)-f(((1<<u)|(1<<v),))
   if Fraction(b+j,2)!=c[u][v]:raise RuntimeError('tomography')
   r=[row[:] for row in q];r[u][v]^=1;r[v][u]^=1
   if H(r)-h != (c[v][u]-c[u][v])*(2*q[u][v]-1):raise RuntimeError('edge derivative')
   flipchecks+=1
  batteries+=1
records=json.loads((ROOT/'score_root_hostile.json').read_text());hostiles=[]
for record in records:
 masks=record['out_masks'];q=[[int(masks[u]>>v&1) for v in range(8)] for u in range(8)]
 if any(q[u][u] or any(q[u][v]+q[v][u]!=1 for v in range(u)) for u in range(8)):raise RuntimeError('not tournament')
 h=H(q); ch=H([r[:7] for r in q[:7]]);scores=list(map(sum,q));delta=[Fraction(s,1)-Fraction(7,2) for s in scores]
 if h!=record['H'] or ch!=record['core_H'] or scores!=record['scores']:raise RuntimeError('hostile mismatch')
 if sum(delta)!=0 or sum(x*x for x in delta)!=2 or any(x.denominator!=2 for x in delta):raise RuntimeError('E8 root')
 if [sum(r[:7]) for r in q[:7]]!=[3]*7:raise RuntimeError('not regular core')
 hostiles.append(dict(H=h,core_H=ch,scores=scores,deviation=[str(x) for x in delta],norm_squared=2,permutations=40320))
if hostiles[0]['scores']!=hostiles[1]['scores'] or hostiles[0]['H']==hostiles[1]['H']:raise RuntimeError('not same point hostile')
res=dict(independent_permutation_one_port_checks=checks,tomography_batteries=batteries,single_edge_flip_checks=flipchecks,hostiles=hostiles)
(OUT/'independent.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
# Verify the claimed construction rather than just trusting serialized masks.
a=[[int((v-u)%7 in {1,2,4}) for v in range(7)] for u in range(7)]
b=[r[:] for r in a]
for u,v in [(0,1),(1,3),(3,0)]:
 if not b[u][v]:raise RuntimeError('triangle orientation')
 b[u][v]=0;b[v][u]=1
for q,record in zip([a,b],records):
 built=extend(q,(7,))
 if [sum(x<<v for v,x in enumerate(row)) for row in built]!=record['out_masks']:raise RuntimeError('construction mismatch')
# Common realization checked against a separate permutation enumerator, not endpoint DP.
import sys,random
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'session27-response'))
import equivalence as eq
rng=random.Random(2702);ac=0
for k,maxm in [(0,7),(1,5),(2,4),(3,3),(4,3)]:
 ps,size,matrices,alpha=eq.realization(k)
 for _ in range(20):
  m=rng.randrange(maxm+1);q=core(k,rng.randrange(1<<(k*(k-1)//2)));w=tuple(rng.randrange(1<<k) for _ in range(m))
  v=alpha
  for mask in w:v=eq.advance(v,matrices[mask])
  f=eq.final(q,ps,size)
  if sum(x*f[s] for s,x in v.items())!=H(extend(q,w)):raise RuntimeError('automaton vs permutation')
  ac+=1
res['automaton_vs_permutation_checks']=ac
(OUT/'independent.json').write_text(json.dumps(res,indent=2)+'\n')
print('Common realization vs permutation:',ac,'PASS; hostile construction verified.')
