"""Exact finite controls for the proof packet; no assert statements."""
from itertools import combinations, permutations, product
from fractions import Fraction as Q
from math import gcd, isqrt, factorial, prod
from collections import Counter, defaultdict
from pathlib import Path
import json
P=Path(__file__).resolve().parent
checks=0
def require(v,msg=''):
 global checks
 checks+=1
 if not v: raise RuntimeError(msg)
def rank(vs):
 b={}
 for x in vs:
  while x:
   p=x.bit_length()-1
   if p in b:x^=b[p]
   else:b[p]=x;break
 return len(b)
def syn(w,cols):
 s=0
 for i,c in enumerate(cols):
  if w>>i&1:s^=c
 return s
def sigma(n):
 return sum(d for d in range(1,isqrt(n)+1) if n%d==0 for d in ([d] if d*d==n else [d,n//d]))
def sf(n):
 r=1;p=2
 while p*p<=n:
  e=0
  while n%p==0:n//=p;e^=1
  if e:r*=p
  p+=1
 return r*n
def prime(n):
 return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def code(cols,s=0):return [w for w in range(1<<len(cols)) if syn(w,cols)==s]
# User's four erased positions in inherited composite bank.
atoms=[4,2,5,3,7,13,11]; cols=[4,2,3,0,1,5,2]
E=[0,1,2,5];K=[3,4,6]
legal=lambda w:not(w&1 and w&2)
erasures=[]
for bits in range(8):
 fixed=sum(((bits>>j)&1)<<i for j,i in enumerate(K))
 raw=[fixed|sum(((u>>j)&1)<<i for j,i in enumerate(E)) for u in range(16)]
 raw=[w for w in raw if syn(w,cols)==0]
 good=[w for w in raw if legal(w)]
 erasures.append({'fixed': [atoms[i] for i in K if fixed>>i&1], 'parity_products':[prod(atoms[i] for i in range(7) if w>>i&1) for w in raw], 'legal_products':[prod(atoms[i] for i in range(7) if w>>i&1) for w in good]})
 require(len(raw)==2)
require(any(r['legal_products']==[70,364] or set(r['legal_products'])=={70,364} for r in erasures))
# Explicit Hamming and extended affine banks; scan d*t^2-1 deterministically.
def bank(dim,affine=False):
 base=[2,3,5,7,11]
 labels=list(range(1,1<<dim)) if not affine else [1|(v<<1) for v in range(1<<dim)]
 rows=[];used=set()
 for v in labels:
  d=prod(base[i] for i in range(dim+int(affine)) if v>>i&1)
  t=1
  while not prime(d*t*t-1) or d*t*t-1 in used:t+=1
  p=d*t*t-1;used.add(p)
  require(sf(p+1)==d)
  rows.append({'column':v,'prime':p,'squarefree_part':d,'multiplier':t})
 return rows
banks={}
for name,d,a in [('hamming7',3,False),('hamming15',4,False),('affine16',4,True)]:
 rows=bank(d,a);cs=[r['column'] for r in rows];C=code(cs)
 banks[name]={'rows':rows,'rank':rank(cs),'dimension':len(cs)-rank(cs),'enumerator':dict(sorted(Counter(w.bit_count() for w in C).items()))}
 require(len(C)==1<<(len(cs)-rank(cs)))
# Seven-prime Fano chart, cuts: exact minimal number of future languages.
cols=list(range(1,8));C=code(cols);D=code(cols,1)
cuts=Counter()
for mask in range(128):
 L=[i for i in range(7) if mask>>i&1];R=[i for i in range(7) if not(mask>>i&1)]
 left={sum(((w>>i)&1)<<j for j,i in enumerate(L)) for w in C}
 rights=defaultdict(set)
 for w in C:
  a=sum(((w>>i)&1)<<j for j,i in enumerate(L));b=sum(((w>>i)&1)<<j for j,i in enumerate(R));rights[a].add(b)
 lam=rank([cols[i] for i in L])+rank([cols[i] for i in R])-3
 require(len({frozenset(x) for x in rights.values()})==1<<lam)
 cuts[(len(L),lam)]+=1
# Orthant-code exact moments through degree four, compared with dual formula.
def moment(words,m):
 def compositions(total,n):
  if n==1:yield(total,);return
  for v in range(total+1):
   for rest in compositions(total-v,n-1):yield(v,)+rest
 ans=Q(0)
 for es in compositions(m,7):
  odd=sum((e%2)<<i for i,e in enumerate(es))
  corr=Q(sum((-1)**((w&odd).bit_count()) for w in words),len(words))
  ans+=Q(factorial(m),prod(factorial(e) for e in es))*corr/prod(e+1 for e in es)
 return ans
moments={}
for s,words in [(0,C),(1,D)]:
 for size in range(4):
  for inds in combinations(range(7),size):
   marg=Counter(tuple(w>>i&1 for i in inds) for w in words)
   require(len(marg)==1<<size and len(set(marg.values()))==1)
 moments[str(s)]=[str(moment(words,m)) for m in range(5)]
require(moments['0']==['1','0','7/3','0','259/10'])
require(moments['1']==['1','0','7/3','0','139/10'])
def box_cdf(x,widths):
 if x<=0:return Q(0)
 ans=Q(0);n=len(widths)
 for mask in range(1<<n):
  t=x-sum(widths[i] for i in range(n) if mask>>i&1)
  if t>0:ans+=(-1)**mask.bit_count()*t**n
 return ans/(factorial(n)*prod(widths))
def tail(words,t,widths):
 return sum((box_cdf(sum(a for i,a in enumerate(widths) if not(w>>i&1))-t,widths) for w in words),Q(0))/len(words)
for words,h in [(C,7),(D,6)]:
 for delta in [Q(1,2),Q(1,3),Q(1,10)]:
  require(tail(words,Q(h)-delta,[Q(1)]*7)==delta**7/(16*factorial(7)))
require(tail(D,Q(13,2),[Q(1)]*7)==0)
require(tail(C,Q(13,2),[Q(1)]*7)==Q(1,10321920))
# Affine geometry: Steiner S(3,4,8), its 56 pointed planes.
planes=code([8|v for v in range(8)])
planes=[w for w in planes if w.bit_count()==4]
triples=list(combinations(range(8),3))
for t in triples:
 mask=sum(1<<i for i in t)
 contain=[p for p in planes if p&mask==mask]
 require(len(contain)==1)
 missing=(contain[0]^mask).bit_length()-1
 require(missing==t[0]^t[1]^t[2])
require(len(planes)==14 and len(triples)==56)
# The 56 order-three matrices of GL(3,2), each with a unique 1+2 split.
matrices=[]
def apply(M,x):return syn(x,M)
for M in product(range(1,8),repeat=3):
 if rank(M)<3:continue
 if all(apply(M,apply(M,apply(M,x)))==x for x in range(8)) and M!=(1,2,4): matrices.append(M)
require(len(matrices)==56)
# Exact equivariant bridge: oriented Fano triangles -> GL3 order3 -> 56 matchings.
def matching(M):
 K={apply(M,x)^x for x in range(8)}
 fixed=[x for x in range(1,8) if apply(M,x)==x]
 require(len(K)==4 and len(fixed)==1)
 a=fixed[0]
 return tuple(sorted(tuple(sorted((x,a^apply(M,x)))) for x in K))
def all_matchings(xs):
 if not xs:yield ();return
 x=xs[0]
 for i in range(1,len(xs)):
  y=xs[i]
  for rest in all_matchings(xs[1:i]+xs[i+1:]):yield tuple(sorted(((x,y),)+rest))
MM={matching(M):M for M in matrices}
allM=list(all_matchings(list(range(8))))
directions=Counter(len({x^y for x,y in M}) for M in allM)
require(directions=={1:7,2:42,4:56})
require(set(MM)=={M for M in allM if len({x^y for x,y in M})==4})
triangles={}
for a,b,c in permutations(range(1,8),3):
 if rank([a,b,c])<3:continue
 tri=min((a,b,c),(b,c,a),(c,a,b))
 ma=tuple(sorted(tuple(sorted(e)) for e in [(0,a^b^c),(a,a^b),(b,b^c),(c,c^a)]))
 require(ma in MM)
 T=MM[ma];require(apply(T,a)==b and apply(T,b)==c and apply(T,c)==a)
 triangles[tri]=ma
require(len(triangles)==56 and len(set(triangles.values()))==56)
# Equivariance under every GL3 matrix, for every one of the 56 objects.
GL=[M for M in product(range(1,8),repeat=3) if rank(M)==3]
for G in GL:
 for tri,ma in triangles.items():
  a,b,c=[apply(G,x) for x in tri]
  newtri=min((a,b,c),(b,c,a),(c,a,b))
  mapped=tuple(sorted(tuple(sorted((apply(G,x),apply(G,y)))) for x,y in ma))
  require(triangles[newtri]==mapped)
# Non-equivalence control: affine triples with fixed origin split into 21,7,28 GL orbits.
triple_types=Counter('contains_zero' if 0 in T else ('line' if T[0]^T[1]^T[2]==0 else 'triangle') for T in combinations(range(8),3))
require(triple_types=={'contains_zero':21,'line':7,'triangle':28})
# Native unit-distance boundary control for the orientation bit.
tri=(1,2,4);rev=(1,4,2);M=triangles[tri];Mi=triangles[rev]
coords={i:(10*i+10,0) for i in range(8)};coords[1]=(0,0);coords[3]=(1,0)
edges=[(i,j) for i,j in combinations(range(8),2) if sum((a-b)**2 for a,b in zip(coords[i],coords[j]))==1]
require(edges==[(1,3)])
def contraction(M,q):
 count=0
 for colors in product(range(q),repeat=4):
  assigned={v:colors[j] for j,e in enumerate(M) for v in e}
  count+=all(assigned[a]!=assigned[b] for a,b in edges)
 return count
geometric_counts=[contraction(N,6) for N in [M,Mi]]
require(geometric_counts==[0,1080])
# Actual Collatz parity cylinders at an unchanged seven-step integer clock.
parity_sources={}
for n in range(1,129):
 x=n;w=0;orbit=[n]
 for j in range(7):
  b=x%2;w|=b<<j;x=(3*x+1)//2 if b else x//2;orbit.append(x)
 require(w not in parity_sources)
 parity_sources[w]={'source':n,'target':x,'first_descent':next((j for j,v in enumerate(orbit[1:],1) if v<n),None)}
collatz_cosets=[]
for sy in range(8):
 words=code(cols,sy);rows=[parity_sources[w] for w in words]
 collatz_cosets.append({'syndrome':sy,'sources':[r['source'] for r in rows], 'descend_by7':sum(r['first_descent'] is not None for r in rows),'endpoint_below_source':sum(r['target']<r['source'] for r in rows),'target_sum':sum(r['target'] for r in rows)})
require(len({r['endpoint_below_source'] for r in collatz_cosets})>1)
# Square-sum Hamilton paths: enumerate bounded witnesses for balanced double decoding.
ham_receipts=[]
for n in [15,16,17,23,25,26,27,28,29,30,31,32]:
 adj={i:[j for j in range(1,n+1) if i!=j and isqrt(i+j)**2==i+j] for i in range(1,n+1)}
 found=[];nodes=[0]
 def dfs(path,used):
  nodes[0]+=1
  if len(found)>=100 or nodes[0]>1000000:return
  if len(path)==n:found.append(path);return
  for v in sorted((j for j in adj[path[-1]] if j not in used),key=lambda j:sum(k not in used for k in adj[j])):
   dfs(path+[v],used|{v})
 starts=sorted(adj,key=lambda j:len(adj[j]))
 for start in starts:
  dfs([start],{start})
  if len(found)>=100 or nodes[0]>1000000:break
 amb=[]
 if n%2==0:
  for x in found:
   m=n//2;c=Q(sum(x[1::2])-sum(x[::2]),m)
   if c and c.denominator==1:
    y=[a+(int(c) if i%2==0 else -int(c)) for i,a in enumerate(x)]
    if sorted(y)==list(range(1,n+1)):
     require(all(x[i]+x[i+1]==y[i]+y[i+1] for i in range(n-1)))
     amb.append({'x':x,'y':y,'shift':int(c)})
 ham_receipts.append({'n':n,'paths_sampled':len(found),'dfs_nodes':nodes[0], 'capped':nodes[0]>1000000 or len(found)>=100,'ambiguous':amb[:1]})
# Exhaustive permutation edge-sum decoding for fixed connected path through n=8.
permutation_summary=[]
for n in range(2,9):
 fibres=defaultdict(list)
 for x in permutations(range(1,n+1)):
  fibres[tuple(x[i]+x[i+1] for i in range(n-1))].append(x)
 mx=max(map(len,fibres.values()));require(mx<=(1 if n%2 else 2))
 permutation_summary.append({'n':n,'fibres':len(fibres),'max_fibre':mx,'two_element_fibres':sum(len(v)==2 for v in fibres.values())})
# Borwein exact first corner at denominators 15 and 113 (not numerical integration).
borwein=[]
for t,m in [(1,7),(2,56)]:
 widths=[Q(1,2*j+1) for j in range(1,m+1)];eps=sum(widths)-t
 require(eps>0 and eps-widths[-1]<=0 and eps<2*min(widths))
 rational=eps**m/(2**m*factorial(m)*prod(widths))
 if m==7:require(rational==Q(6879714958723010531,935615849440640907310521750000))
 borwein.append({'threshold':t,'first_m':m,'last_m':m-1,'first_denominator':2*m+1,'epsilon':str(eps),'deficit_over_pi':str(rational)})
result={'checks':checks,'guarded_erasures':erasures,'banks':banks,'fano_cut_census':{str(k):v for k,v in sorted(cuts.items())},'orthant_moments':moments,'orthant_tail_at_6.5':{'C0':'1/10321920','C1':'0'},'affine_planes':len(planes),'pointed_planes':56,'GL3_order3':len(matrices),'matching_direction_counts':dict(directions),'equivariant_triangles':len(triangles),'affine_triple_orbits':dict(triple_types),'unit_distance_orientation_counts':geometric_counts,'collatz_cosets':collatz_cosets,'square_sum_sample':ham_receipts,'permutation_decoding':permutation_summary,'borwein':borwein}
(P/'receipts.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['banks','borwein']},indent=2))
