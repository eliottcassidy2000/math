#!/usr/bin/env python3
"""Standalone exact controls; no research-tree imports or assert statements."""
import itertools as it
import json, math
from collections import Counter

def require(ok,label):
 if not ok: raise RuntimeError(label)

def seed(bit):
 mask=(8,10)[bit]; a=[0]*5
 for j,(x,y) in enumerate(it.combinations(range(5),2)):
  tail,head=(x,y) if mask>>j&1 else (y,x)
  a[tail]|=1<<head
 return tuple(a)

def hp(a):
 n=len(a); dp={(1<<i,i):1 for i in range(n)}
 for mask in range(1,1<<n):
  for last in range(n):
   count=dp.get((mask,last),0)
   choices=a[last]&~mask
   while count and choices:
    b=choices&-choices; choices-=b; j=b.bit_length()-1
    dp[mask|b,j]=dp.get((mask|b,j),0)+count
 return sum(dp.get(((1<<n)-1,i),0) for i in range(n))

def covers(a):
 out=set(); n=len(a)
 for p in it.permutations(range(n)):
  for c in range(1<<(n-1)):
   paths=[]; start=0
   for j in range(n-1):
    if c>>j&1: paths.append(p[start:j+1]); start=j+1
    elif not a[p[j]]>>p[j+1]&1: break
   else:
    paths.append(p[start:]); out.add(tuple(sorted(paths)))
 return dict(sorted(Counter(map(len,out)).items()))

def extend(a,p=0):
 n=len(a)
 return tuple(v|((1<<n) if i!=p else 0) for i,v in enumerate(a))+(1<<p,)

def ordinal(parts):
 n=sum(map(len,parts)); a=[]; start=0
 for part in parts:
  end=start+len(part); later=((1<<n)-1)^((1<<end)-1)
  a.extend((v<<start)|later for v in part); start=end
 return tuple(a)

def arc_distance(a,b):
 return sum((x^y).bit_count() for x,y in zip(a,b))//2

def rank(rows):
 piv={}
 for x in rows:
  while x:
   p=x.bit_length()-1
   if p in piv:x^=piv[p]
   else:piv[p]=x;break
 return len(piv)

def basis(d,r):
 return [sum(1<<y for y in range(1<<d) if y&m==m) for m in range(1<<d) if m.bit_count()<=d-r]

def code(b):
 c=[0]
 for x in b:c += [y^x for y in c]
 return c

def eval_rank(b,S):
 return rank([sum(((v>>y)&1)<<j for j,y in enumerate(S)) for v in b])

report={}; seeds=[seed(0),seed(1)]
require(arc_distance(*seeds)==1,'seed distance')
require(not any(all(bool(seeds[0][i]>>j&1)==bool(seeds[1][p[i]]>>p[j]&1) for i,j in it.combinations(range(5),2)) for p in it.permutations(range(5))),'nonisomorphism')
for a in seeds:
 require(hp(a)==9,'seed Hamilton count')
 require(covers(a)=={1:9,2:33,3:31,4:10,5:1},'cover profile')
 for start in range(5):
  seen=1<<start
  for _ in range(5):
   for i in range(5):
    if seen>>i&1:seen|=a[i]
  require(seen==31,'strong seed')
require([hp(extend(a)) for a in seeds]==[17,19],'local observer')
for bits in it.product(range(2),repeat=2):
 parts=[seeds[b] for b in bits]
 require(hp(ordinal(parts))==81,'literal two-block product')
 for y in range(2):
  p=parts.copy();p[y]=extend(p[y])
  require(hp(ordinal(p))==9*(17+2*bits[y]),'literal contextual observer')
report['literal']={'seed_rows':seeds,'seed_H':9,'query_H':[17,19],'two_block_graphs':12}
records=[]
for d in range(2,6):
 for r in range(2,d+1):
  b=basis(d,r); M=1<<d; K=len(b)
  bank=[y for y in range(M) if y.bit_count()<=d-r]
  require(len(bank)==K and eval_rank(b,bank)==K,'optimal bank')
  item={'d':d,'r':r,'K':K,'M':M,'distance':1<<r}
  if K<=11:
   words=code(b); require(min(x.bit_count() for x in words if x)==1<<r,'distance')
   for q in range(min(M,6)+1):
    S=list(range(q)); counts=Counter(sum(((w>>y)&1)<<j for j,y in enumerate(S)) for w in words)
    rk=eval_rank(b,S)
    require(len(counts)==1<<rk and set(counts.values())=={1<<(K-rk)},'exact posterior fibers')
   item['enumerated_words']=len(words)
  # Direct labelled arc distance on constructed 5M-vertex objects.
  samples=code(b[:min(K,5)])
  origin=ordinal([seeds[0]]*M)
  for w in samples:
   graph=ordinal([seeds[w>>y&1] for y in range(M)])
   require(arc_distance(origin,graph)==w.bit_count(),'literal code isometry')
  records.append(item)
report['codes']=records
# Exhaustive raw-arc single-error decoder on the d=3,r=2 affine code:
# 40 vertices, 780 possible corruptions per codeword, even non-code arcs.
d=3;r=2; words=code(basis(d,r)); graphs=[ordinal([seeds[w>>y&1] for y in range(8)]) for w in words]
checks=0
for idx,a in enumerate(graphs):
 for i,j in it.combinations(range(40),2):
  corrupt=list(a);corrupt[i]^=1<<j;corrupt[j]^=1<<i
  distances=[arc_distance(corrupt,g) for g in graphs]
  require(distances[idx]==1 and sum(x==min(distances) for x in distances)==1,'arbitrary single arc recovery')
  checks+=1
report['arbitrary_single_arc_decoder_checks']=checks
# Restricted local probe promise matters: aggregate query returns subset weight.
for bits in it.product(range(2),repeat=2):
 require(hp(ordinal([extend(seeds[b]) for b in bits]))==17**(2-sum(bits))*19**sum(bits),'aggregate observer')
report['status']='PASS'
print(json.dumps(report,indent=2,sort_keys=True))
