from itertools import permutations
from collections import Counter
import json

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
def orbit(p):return {tuple((p[(i-a)%5]+b)%5 for i in range(5)) for a in range(5) for b in range(5)}
def adj(p):
 A=[[False]*10 for _ in range(10)]
 for i in range(10):
  for j in range(10):
   if i//5==j//5:A[i][j]=(j-i)%5 in (1,2)
   elif i<5:A[i][j]=(j-5)==p[i]
   else:A[i][j]=i-5!=p[j]
 return A

def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S in range(1,1<<n):
  for i in range(n):
   if D[S][i]:
    for j in range(n):
     if not(S>>j&1) and A[i][j]:D[S|1<<j][j]+=D[S][i]
 return sum(D[-1])
seen=set();rows=[];checks=0
# Direct automorphism classification of the internal regular5 tournament.
R=adj(tuple(range(5)));auts=[]
for q in permutations(range(5)):
 if all(R[i][j]==R[q[i]][q[j]] for i in range(5) for j in range(5)):auts.append(q)
if set(auts)!={tuple((i+a)%5 for i in range(5)) for a in range(5)}:raise ValueError('internal automorphism group')
for p in permutations(range(5)):
 o=orbit(p)
 if min(o) in seen:continue
 seen.add(min(o));p=min(o);A=adj(p)
 if sorted(map(sum,A))!=[3]*5+[6]*5:raise ValueError('intrinsic cells')
 if len({sign(q) for q in o})!=1:raise ValueError('parity gauge')
 for q in o:
  if min(orbit(q))!=p:raise ValueError('orbit closure')
  checks+=1
 inv=[p.index(i) for i in range(5)]
 cp=min(orbit(tuple(-inv[-i%5]%5 for i in range(5))))
 J=sum(1 if (p[j]-p[i])%5 in (1,2) else -1 for i in range(5) for j in range(5) if (j-i)%5 in (1,2))
 K=sum(1 if (p[(i+1)%5]-p[i])%5 in (1,2) else -1 for i in range(5))
 rows.append(dict(permutation=p,orbit_size=len(o),parity=sign(p),H=h(A),converse=cp,J=J,K=K))
if len(rows)!=8 or Counter(r['orbit_size'] for r in rows)!={5:4,25:4}:raise ValueError('orbit census')
if Counter(r['parity'] for r in rows)!={-1:4,1:4}:raise ValueError('parity collapse')
merged={min(tuple(r['permutation']),tuple(r['converse'])) for r in rows}
if len(merged)!=7:raise ValueError('merged count')
for modulus in (8,16):
 buckets={}
 for r in rows:buckets.setdefault(r['J']%modulus,set()).add(min(tuple(r['permutation']),tuple(r['converse'])))
 if (all(len(b)==1 for b in buckets.values())) != (modulus==16):raise ValueError('J lift')
if len({(r['J']%8,r['K']%4) for r in rows})!=8:raise ValueError('two moment mod8 decoder')
for a in range(8):
 for b in range(8):
  if len({(a*r['J']+b*r['K'])%8 for r in rows})>4:raise ValueError('one linear mod8 obstruction')
print(json.dumps(dict(status='PASS',universe='all120permutations of regular5-cell matching block',checks=checks,internal_automorphisms=auts,classes=rows,merged_classes=7,mod8_pair_decoder='(J mod8,K mod4) injective on all8ordinarytypes',linear_mod8_obstruction='every aJ+bK has atmost4residues'),indent=2))
