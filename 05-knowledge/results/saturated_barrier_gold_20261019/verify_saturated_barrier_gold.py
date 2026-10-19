#!/usr/bin/env python3
"""Exact saturated-barrier and two-free-port Gold certificates.
Standard library only. No imported discovery arrays or earlier verifier.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product,permutations
import hashlib,json

Q=(0,0,0,0,7,13,14,12,206,5,685,535)
K=7;N=len(Q);SAT=(1,4,5,6,7,9);TREE=((7,9),(4,5))
checks=0

def require(ok,why):
 global checks
 checks+=1
 if not ok:raise ValueError(why)

def inc(p,a,b):return a!=b and not((p[a]>>b&1)or(p[b]>>a&1))
def closed(p):
 n=len(p)
 return all(not(p[a]>>a&1)and all(not(p[a]>>b&1)or not(p[b]&~p[a])for b in range(n))for a in range(n))
def connected(p,vertices):
 seen={vertices[0]}
 while True:
  nxt=seen|{v for v in vertices if any(inc(p,v,w)for w in seen)}
  if nxt==seen:return len(seen)==len(vertices)
  seen=nxt

def model(alpha,beta,deleted=()):
 kept=[q for q in range(N)if q not in deleted]
 h=max((alpha[q]for q in kept),default=0);j=max((beta[q]for q in kept),default=0)
 labels=kept+list(range(N,N+h+j));ix={q:i for i,q in enumerate(labels)};p=[0]*len(labels)
 for q in kept:
  p[ix[q]]=sum(1<<ix[r]for r in kept if Q[q]>>r&1)
  p[ix[q]]|=sum(1<<ix[N+i]for i in range(h-alpha[q]))
 for i in range(h):p[ix[N+i]]=sum(1<<ix[N+z]for z in range(i))
 for i in range(j):
  p[ix[N+h+i]]=sum(1<<ix[N+z]for z in range(h+i))
  p[ix[N+h+i]]|=sum(1<<ix[q]for q in kept if i>=beta[q])
 require(closed(p),'model is transitively closed and acyclic')
 require(all(((p[ix[q]]>>ix[r])&1)==((Q[q]>>r)&1)for q in kept for r in kept),'induced core unchanged')
 return labels,tuple(p)

def count(labels,p,events=()):
 ix={q:i for i,q in enumerate(labels)};pp=list(p)
 for a,b in events:pp[ix[b]]|=1<<ix[a]
 full=(1<<len(pp))-1
 @lru_cache(None)
 def dp(s):
  if s==full:return 1
  return sum(dp(s|1<<q)for q,d in enumerate(pp)if not(s>>q&1)and not(d&~s))
 return dp(0)

def cells(alpha,beta,deleted=()):
 labels,p=model(alpha,beta,deleted)
 E=count(labels,p);F=count(labels,p,(TREE[0],));G=count(labels,p,(TREE[1],));J=count(labels,p,TREE)
 cc=(J,F-J,G-J,E-F-G+J)
 require(E>0 and min(cc)>=0 and sum(cc)==E,'nonnegative four-cell partition')
 return cc,labels,p

def margins(cc):
 E=sum(cc)
 return tuple(E-cc[2*(r//2)]-cc[2*(r//2)+1]-cc[r]for r in range(4))

def profiles(a,b,c,f,g,h):
 return ((a,0,b,c,0,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,f,0,g,h))

require(closed(Q),'core is a poset')
DEG=tuple(sum(inc(Q,a,b)for b in range(N))for a in range(N))
DELTA=tuple(K-d for d in DEG)
require(DEG==(6,7,3,6,7,7,7,7,6,7,5,6),'core degrees')
require(tuple(q for q in range(N)if DEG[q]==K)==SAT,'exact saturated set')
require(connected(Q,SAT),'saturated incomparability graph connected')
require(all(q in SAT or any(inc(Q,q,s)for s in SAT)for q in range(N)),'saturated set dominates entire core')
require(Q[2]==0 and all(not(Q[q]>>10&1)for q in range(N)),'free ports are minimal and maximal')
require(all(q in SAT or q in (0,2,3,8,10,11)for q in range(N)),'core port partition')
require(all(any(Q[q]>>s&1 for s in SAT)for q in (8,10,11)),'upper ports above barrier vertices')
require(all(any(Q[s]>>q&1 for s in SAT)for q in (0,2,3)),'lower ports below barrier vertices')

# Generic bounded monotone profiles are checked against the explicit rectangle.
def general_profiles(side):
 eligible=[q for q in range(N)if q not in SAT and (not any(Q[q]>>s&1 for s in SAT)if side==0 else not any(Q[s]>>q&1 for s in SAT))]
 result=[]
 for values in product(*(range(DELTA[q]+1)for q in eligible)):
  t=[0]*N
  for q,v in zip(eligible,values):t[q]=v
  if sum(v>0 for v in t)>K:continue
  if any((t[a]<t[b]if side==0 else t[a]>t[b])for a in range(N)for b in range(N)if Q[b]>>a&1):continue
  result.append(tuple(t))
 return sorted(result)
A=sorted(profiles(a,b,c,0,0,0)[0]for a,b,c in product(range(2),range(5),range(2)))
B=sorted(profiles(0,0,0,f,g,h)[1]for f,g,h in product(range(2),range(3),range(2)))
require(general_profiles(0)==A and general_profiles(1)==B,'complete profile rectangle')
require((len(A),len(B))==(20,12),'profile counts')

# An independent core-word/gap enumeration provides every bounded atom cell count.
@lru_cache(None)
def gap(t,last=1000000):
 if not t:return 1
 return sum(gap(t[1:],z)for z in range(min(t[0],last)+1))
FIBERS=Counter();words=0

def visit(word,placed):
 global words
 if len(word)==N:
  words+=1;pos={q:i for i,q in enumerate(word)};left=min(pos[s]for s in SAT);right=max(pos[s]for s in SAT)
  pre=tuple(word[:left]);suf=tuple(reversed(word[right+1:]));cell=2*(not(pos[7]<pos[9]))+(not(pos[4]<pos[5]))
  require(set(pre)<={0,2,3} and set(suf)<={8,10,11},'only lower/upper ports enter the boundary words')
  FIBERS[pre,suf,cell]+=1
  return
 for q in range(N):
  if not(placed>>q&1)and not(Q[q]&~placed):visit(word+[q],placed|1<<q)
visit([],0)
require(words==61180,'core extension total')

def gap_cells(a,b):
 cc=[0]*4
 for (pre,suf,r),mult in FIBERS.items():cc[r]+=mult*gap(tuple(a[q]for q in pre))*gap(tuple(b[q]for q in suf))
 return tuple(cc)

bounded_min=[None]*4;bounded_witness=[None]*4;atom_digest=hashlib.sha256();orders=[];balanced=[None,None]
for ia,a in enumerate(A):
 for ib,b in enumerate(B):
  require(all(x+y<=d for x,y,d in zip(a,b,DELTA)),'shared degree allocations compatible')
  cc,labels,p=cells(a,b);E=sum(cc);mm=margins(cc);orders.append(len(p))
  require(max(sum(inc(p,q,r)for r in range(len(p)))for q in range(len(p)))<=K,'every bounded atom seven-thin')
  require(cc==gap_cells(a,b),'independent core-word gap counts agree with ideal DP')
  for r,m in enumerate(mm):
   require(9*m>E,'bounded margin strictly exceeds E/9')
   v=Fraction(m,E)
   if bounded_min[r]is None or v<bounded_min[r]:bounded_min[r]=v;bounded_witness[r]=[ia,ib,E,list(cc)]
  first=Fraction(cc[0]+cc[1],E)
  balanced[0]=first if balanced[0]is None else min(balanced[0],first)
  balanced[1]=first if balanced[1]is None else max(balanced[1],first)
  atom_digest.update((json.dumps([list(a),list(b),E,list(cc)],separators=(',',':'))+'\n').encode())
require(list(map(str,bounded_min))==['78113/453252','25363/180473','50695/453252','4768/24065'],'exact bounded minima')

# The free-port coefficients are actual deletion models, independently equal
# to finite differences of full-core models at gap caps 0,1,2.
free_min=[None]*4;free_witness=[None]*4;free_digest=hashlib.sha256();free_orders=[];generator_count=0
for a,c,f,h in product(range(2),repeat=4):
 grid={}
 for b,g in product(range(3),repeat=2):
  aa,bb=profiles(a,b,c,f,g,h);grid[b,g]=cells(aa,bb)[0]
 for lower,upper in product(range(3),repeat=2):
  # States 0 and 1 retain the free port at that cap. State 2 deletes it.
  deleted=tuple(([2]if lower==2 else[])+([10]if upper==2 else[]))
  aa,bb=profiles(a,lower if lower<2 else 0,c,f,upper if upper<2 else 0,h)
  cc,labels,p=cells(aa,bb,deleted);E=sum(cc);mm=margins(cc);free_orders.append(len(p));generator_count+=1
  dl=[(lower,1)]if lower<2 else[(2,1),(1,-1)]
  du=[(upper,1)]if upper<2 else[(2,1),(1,-1)]
  diff=tuple(sum(x*y*grid[i,j][r]for i,x in dl for j,y in du)for r in range(4))
  require(cc==diff,'deletion generator equals the full-core finite difference')
  for r,m in enumerate(mm):
   require(50*m>3*E,'free-port margin strictly exceeds 3E/50')
   v=Fraction(m,E)
   if free_min[r]is None or v<free_min[r]:free_min[r]=v;free_witness[r]=[a,c,f,h,lower,upper,E,list(cc)]
  free_digest.update((json.dumps([a,c,f,h,lower,upper,E,list(cc)],separators=(',',':'))+'\n').encode())
require(generator_count==144,'entire two-free-port generator universe')
require(list(map(str,free_min))==['771/6245','3231/24464','386/6245','3116/16323'],'exact free-port minima')

# Verify the elementary affine-tail identity on every ordered subset of a
# three-port alphabet, with general bounded threshold r, not just this core.
affine_controls=0
for threshold in range(4):
 for x,y in product(range(threshold+1),repeat=2):
  caps={0:x,1:threshold,2:y}
  for length in range(4):
   for w in permutations(range(3),length):
    for b in range(threshold,threshold+8):
     actual=dict(caps);actual[1]=b
     lhs=gap(tuple(actual[q]for q in w));base=gap(tuple(caps[q]for q in w))
     derivative=gap(tuple(caps[q]for q in w[1:]))if w and w[0]==1 else 0
     require(lhs==base+(b-threshold)*derivative,'general affine-tail gap identity')
     affine_controls+=1

# A genuine unbounded family approaches the best recession margin.
sharpness=[]
for n in (1,2,4,7):
 aa,bb=profiles(1,n,0,1,n,0)
 cc,labels,p=cells(aa,bb);E=sum(cc);m=margins(cc)[2]
 require(E==12490*n*n+78784*n+122888,'sharpness extension polynomial')
 require(m==772*n*n+8213*n+18243,'sharpness Gold-margin polynomial')
 require(cc==gap_cells(aa,bb),'unbounded sharpness family independent gap counts')
 require(6245*m>386*E,'finite family lies strictly above its limiting margin')
 sharpness.append([n,E,m])
require(Fraction(772,12490)==Fraction(386,6245),'sharp limiting constant')

# Hostile controls: neither connectedness nor domination may be discarded.
# Two two-antichains separated by an inserted singleton.
P=(0,0,3,7,7);Q0=(0,0,3,3);S0=(0,1,2,3)
require(max(sum(inc(P,a,b)for b in range(5))for a in range(5))==1,'disconnected hostile is one-thin')
require(not connected(Q0,S0),'hostile saturated set disconnected')
require(all(P[2]>>s&1 for s in (0,1))and all(P[s]>>2&1 for s in (3,4)),'outside vertex lies between saturated components')
# r < z < {s,t}; the retained core {r,s,t} has non-dominating S={s,t}.
P=(0,1,3,3);Q0=(0,1,1)
require(connected(Q0,(1,2))and not any(inc(Q0,0,s)for s in (1,2)),'connected hostile fails domination')
require(P[1]>>0&1 and all(P[s]>>1&1 for s in (2,3)),'lower exterior can lie above a nondominated core vertex')

result={'status':'PASS','guards':checks,'core':Q,'core_extensions':words,'degrees':DEG,'saturated_barrier':SAT,'tree':TREE,
'bounded_certificate':{'lower_profiles':len(A),'upper_profiles':len(B),'atoms':len(A)*len(B),'orders':[min(orders),max(orders)],'margin':'1/9, strict','exact_cell_minima':list(map(str,bounded_min)),'minimizers':bounded_witness,'first_probability_interval':list(map(str,balanced)),'count_table_sha256':atom_digest.hexdigest()},
'unbounded_certificate':{'lower_free_port':2,'upper_free_port':10,'other_port_cap':1,'generators':generator_count,'orders':[min(free_orders),max(free_orders)],'margin':'3/50, strict','exact_cell_minima':list(map(str,free_min)),'minimizers':free_witness,'count_table_sha256':free_digest.hexdigest()},
'independent_controls':{'boundary_word_fibers':len(FIBERS),'bounded_atoms_recounted_by_gap_words':len(A)*len(B),'deletion_generators_checked_by_full_core_differences':generator_count,'general_affine_gap_controls':affine_controls,'structural_hostiles':2},'sharpness_family':{'extension_polynomial_coefficients_descending':[12490,78784,122888],'FT_margin_polynomial_coefficients_descending':[772,8213,18243],'checked_values':sharpness}}
print(json.dumps(result,indent=2))
