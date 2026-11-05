"""Independent graph checks for incremental reachability and boundary factorization."""
import itertools,json,random,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ribbon'))
from transfer import graph,count,held_karp

def floyd(a):
 r=[sum(1<<j for j in range(len(a))if i==j or a[i][j])for i in range(len(a))]
 for k in range(len(a)):
  for i in range(len(a)):
   if r[i]>>k&1:r[i]|=r[k]
 return r

def incremental(word):
 r=[]
 for n,a in enumerate(word):
  out=sum(1<<i for i in range(max(0,n-3),n)if a>>(i-(n-3))&1)
  ins=((1<<n)-1)^out;new=1<<n
  for i in range(n):
   if out>>i&1:new|=r[i]
  r=[x|new if x&ins else x for x in r]+[new]
 return r

def split(word):
 a=graph(word,3);n=len(a);r=floyd(a)
 if any(r[i]&7 for i in range(n-3,n)):return None
 todo=set(range(n));cs=[]
 while todo:
  i=min(todo);c={j for j in todo if r[i]>>j&1 and r[j]>>i&1};cs.append(c);todo-=c
 cs=sorted(cs,key=lambda c:-sum(a[min(c)][min(d)]for d in cs if c!=d))
 left=max(i for i,c in enumerate(cs)if c&{0,1,2});right=min(i for i,c in enumerate(cs)if c&set(range(n-3,n)))
 assert left<right
 I=sorted(set().union(*cs[:left+1]));M=sorted(set().union(*cs[left+1:right]));J=sorted(set().union(*cs[right:]))
 assert all(a[i][j]for i in I for j in M+J) and all(a[i][j]for i in M for j in J)
 return I,M,J

def induced(a,vs):return [[a[i][j]for j in vs]for i in vs]
def H(a):return held_karp(a)if a else 1

canonical=0
alphas=[range(0,8,1<<max(0,3-j))for j in range(6)]
for word in itertools.product(*alphas):
 assert incremental(word)==floyd(graph(word,3));canonical+=1
rng=random.Random(197);random_reach=0;factor_checks=0;ratios=set()
for n in range(7,13):
 for _ in range(100):
  word=[0 if rng.random()<.6 else rng.randrange(8)for _ in range(n)]
  assert incremental(word)==floyd(graph(word,3));random_reach+=1
  parts=split(word)
  if parts is None:continue
  I,M,J=parts;a=graph(word,3);two=graph(word*2,3)
  q=H(induced(a,M))*H(induced(two,J+[n+i for i in I]));ratios.add(q)
  prefix=[rng.randrange(8)for _ in range(5)]
  hs=[count(prefix+word*k,3)for k in range(1,5)]
  assert all(hs[k+1]==q*hs[k]for k in range(3)),(word,parts,q,hs)
  factor_checks+=1
result={'canonical_order6_reach_checks':canonical,'random_reach_checks':random_reach,
 'literal_SCC_factorization_checks':factor_checks,'distinct_verified_multipliers':sorted(ratios)}
(HERE/'PHYSICAL_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
