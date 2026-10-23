from itertools import product, combinations
from fractions import Fraction as Q
from collections import Counter
from math import gcd,isqrt
from pathlib import Path
import json
out={}
# E8 independent rational coordinate representation.
R=set()
for ix in combinations(range(8),2):
 for a,b in product([-1,1],repeat=2):
  v=[Q(0)]*8; v[ix[0]]=Q(a);v[ix[1]]=Q(b);R.add(tuple(v))
for s in product([-1,1],repeat=8):
 if sum(x<0 for x in s)%2==0:R.add(tuple(Q(x,2) for x in s))
out['grades']=dict(Counter(int(r[0]+r[1]) for r in R))
assert out['grades']=={-2:1,-1:56,0:126,1:56,2:1}
plus={r for r in R if r[0]+r[1]==1};theta=(Q(1),Q(1))+ (Q(0),)*6
assert all(tuple(t-r for t,r in zip(theta,v)) in plus for v in plus)
assert {tuple(a+b for a,b in zip(v,w)) for v in plus for w in plus}&R=={theta}
models=[]
for d in range(1,isqrt(114)+1):
 if 114%d:continue
 p,q=d+1,114//d+1
 if gcd(p,q)!=1:continue
 c=1-Q(6*(p-q)**2,p*q)
 cb=1-Q(3*113**2,115)
 A=Q(57*115,8)*(c-cb)
 Af=Q(3*(p-1)*(p-2)*(q-1)*(q-2)*(p*p*q*q-(p-q)**2),16*p*q)
 assert A==Af
 models.append([p,q,str(c),str(A)])
out['models']=models
# All mod57 root fibres independently grouped by valuation support.
fibres={a:[r for r in range(57) if r*r%57==a] for a in range(57)}
for a,rs in fibres.items():
 if rs:assert len(rs)==2**sum(a%p!=0 for p in [3,19])
out['root_histogram57']=dict(Counter(map(len,fibres.values())))
out['norm57']={'a':5,'b':4,'norm':5**2+2*4**2,'sqrt_minus2_mod57':[r for r in range(57) if (r*r+2)%57==0]}
out['mod18_base2']={y:(2**y+y*y)%9 for y in range(1,19,2)}
# Prime-divisor splitting for primitive odd Leyland numbers, direct factorization independently.
def factor(n):
 f={};p=2
 while p*p<=n:
  while n%p==0:f[p]=f.get(p,0)+1;n//=p
  p+=1
 if n>1:f[n]=1
 return f
small=[]
for x in range(2,9):
 for y in range(x+1,15):
  n=x**y+y**x
  if n>10**8:continue
  f=factor(n)
  if (x+y)%2 and gcd(x,y)==1:
   e=x if x%2==0 else y
   for p in f:assert pow(-e%p,(p-1)//2,p)==1
  if len(f)==2 and all(v==1 and p%4==3 for p,v in f.items()):small.append([x,y,n,f])
out['independent_small_blums']=small
Path('lie-arithmetic-audit').mkdir(exist_ok=True)
Path('lie-arithmetic-audit/independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
# Full Leyland census with trial division independent of author's prime sieve.
full_rows=[]; full_blums=[]
for x in range(2,11):
 for y in range(x+1,40):
  n=x**y+y**x
  if n<=10**12:
   full_rows.append((x,y,n)); f=factor(n)
   if len(f)==2 and all(k==1 and p%4==3 for p,k in f.items()):
    full_blums.append([x,y,n,{str(p):k for p,k in f.items()}])
author=json.loads(Path('arithmetic/receipt.json').read_text())
assert len(full_rows)==author['leyland_rows']==111
assert full_blums==author['blum_rows']
out['full_census_independent_trial_division']={'rows':len(full_rows),'blum_rows':full_blums}
Path('lie-arithmetic-audit/independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print('Full census matches: 111 rows and 12 Blum rows.')
# Literal decoded path lifts, physical central deletion and CRT sign equivariance.
def states(k):return [frozenset(i+1 for i in range(k-1) if mask>>i&1) for mask in range(2**(k-1))]+[None]
def decode(k,I,K):
 A=set(range(k)); B=set(range(k+1,2*k+1))
 I=A if I is None else set(I); K=A if K is None else set(K)
 J={2*k-v for v in A-K}
 return tuple(sorted(I)+sorted(J)+[k]+sorted(A-I)+sorted(B-J))
def valid(k,P):return len(set(P))==2*k+1 and all((a>b if {a,b} in [{0,k},{k,2*k}] else a<b) for a,b in zip(P,P[1:]))
def lifts(S,k):return [(None,0)] if S is None else [(S,0),(S|{k},1)]
checks=0
for k in range(2,6):
 for I,K in product(states(k),repeat=2):
  parent=decode(k,I,K); assert valid(k,parent)
  for a,rs in fibres.items():
   if not rs or ((a%3==0)!=(I is None)) or ((a%19==0)!=(K is None)):continue
   mapped=[]
   for (J,e),(L,f) in product(lifts(I,k),lifts(K,k)):
    child=decode(k+1,J,L); assert valid(k+1,child)
    deletion=tuple(v if v<k else v-1 if v==k+1 else v-2 for v in child if v not in [k,k+2])
    assert deletion==parent
    u=((-1)**e*pow(a%3,1,3))%3;v=((-1)**f*pow(a%19,5,19))%19
    r=(19*u+39*v)%57;assert r*r%57==a
    mapped.append(r)
   assert sorted(mapped)==rs
   checks+=1
out['literal_deletion_to_root_fibres']=checks
Path('lie-arithmetic-audit/independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print('Literal deletion/CRT fibre comparisons:',checks)
