#!/usr/bin/env python3
import itertools as it, math, json, functools, random, collections
from pathlib import Path
P=Path(__file__).resolve().parent.parent
counts=collections.Counter()
def ck(v,s):
 counts[s]+=1
 if not v: raise RuntimeError(s)
def tournament(n,b):
 a=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(it.combinations(range(n),2)):a[u][v]=(b>>k)&1;a[v][u]=1-a[u][v]
 return a
def brute_h(a,vertices=None):
 vs=list(range(len(a))) if vertices is None else vertices
 return sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in it.permutations(vs))
def forest_profile(a):
 # Subset HP by literal permutations, then unordered partitions via least vertex.
 n=len(a); h={s:brute_h(a,[i for i in range(n) if s>>i&1]) for s in range(1,1<<n)}
 @functools.lru_cache(None)
 def f(s):
  if not s:return (1,)+(0,)*n
  out=[0]*(n+1); low=s&-s;t=s
  while t:
   if t&low:
    z=f(s^t)
    for k in range(1,n+1):out[k]+=h[t]*z[k-1]
   t=(t-1)&s
  return tuple(out)
 return [math.factorial(k)*v for k,v in enumerate(f((1<<n)-1))]
def coeff(q):
 # Finite automaton: append unvisited exterior or star after an exterior.
 m=len(q)-1
 @functools.lru_cache(None)
 def walk(used,last):
  out=[0]*(m+2)
  if used==(1<<m)-1:out[0]=1
  if last!=0 and (last==-1 or q[last][0]):
   z=walk(used,0)
   for k in range(m+1):out[k+1]+=z[k]
  for j in range(1,m+1):
   if not used>>(j-1)&1 and (last==-1 or q[last][j]):
    z=walk(used|1<<(j-1),j)
    for k in range(m+2):out[k]+=z[k]
  return tuple(out)
 return walk(0,-1)
def subst(q,a):
 n=len(a);N=n+len(q)-1
 return [[(a[i][j] if i<n and j<n else q[0 if i<n else i-n+1][0 if j<n else j-n+1]) for j in range(N)] for i in range(N)]
# Fano independently constructs each orientation then counts extendable simple paths.
norm=[4,1,3,7,6,5,2]; H=[]
for mask in range(128):
 a=[[int(i!=j and (j-i)%7 in [1,2,4]) for j in range(7)] for i in range(7)]
 for k in range(7):
  if mask>>k&1:
   for i,j in it.combinations([k,(k+1)%7,(k+3)%7],2):a[i][j]^=1;a[j][i]^=1
 def extend(last,remaining):
  if not remaining:return 1
  return sum(extend(j,remaining-{j}) for j in remaining if a[last][j])
 h=sum(extend(i,set(range(7))-{i}) for i in range(7));H.append(h)
 s=0
 for i in range(7):
  if mask>>i&1:s^=norm[i]
 ck(h==171+18*(s==0),'fano_base')
 ck(all(sum(r)==3 for r in a),'fano_regular')
# Independent full Walsh with integer-scaled coefficients, conditional averages and derivatives.
W=[sum(h*(1-2*((mask&b).bit_count()%2)) for mask,h in enumerate(H)) for b in range(128)]
ck(W[0]==22176 and sorted(x for x in W[1:] if x)==[288]*7,'walsh_amplitudes')
ck(all(not W[b] or b.bit_count()==4 for b in range(1,128)),'walsh_degree')
ck(H[0]-H[1]==18,'pointwise_not_blind')
# Weighted one-step choice destroys lumpability unless all direction weights agree.
ck(sum((i+1)*H[1^(1<<i)] for i in range(7))!=sum((i+1)*H[2^(1<<i)] for i in range(7)),'nonuniform_average_hostile')
# All-r code rank and triple relation generating rank, tested through r=7 without cubes.
def rank(rows):
 basis={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k not in basis:basis[k]=x;break
   x^=basis[k]
 return len(basis)
for r in range(2,8):
 n=(1<<r)-1
 triples=[(1<<(a-1))^(1<<(b-1))^(1<<((a^b)-1)) for a in range(1,n+1) for b in range(a+1,n+1)]
 ck(rank(triples)==n-r,'triple_span_all_r_instances')
 ck({sum((a&b).bit_count()%2 for a in range(1,n+1)) for b in range(1,n+1)}=={1<<(r-1)},'dual_weight_all_r_instances')
# Profile extraction is independent of producer permutation-cut method.
factors=[]
known_order={}
for n in range(1,6):
 for mask in range(1<<(n*(n-1)//2)):
  a=tournament(n,mask);f=forest_profile(a)
  ck(f[1]==brute_h(a),'partition_profile_base')
  h=max(1,(n-1)//2); key=(n,tuple(f[1:h+1]))
  ck(key not in known_order or known_order[key]==f,'fixed_order_sufficiency_census')
  known_order[key]=f
  if n==5 and mask in [40,69,76]:ck(f==([0,15,90,210,240,120] if mask==76 else [0,15,78,198,240,120]),'exact_same_H_witness_profiles')
  if n<=3 or (n==5 and mask in [40,76]):factors.append((a,f))
# Context coefficients do not rely on C3 kernel or combinations-of-star-slots.
for m in range(5):
 for b in range(1<<(m*(m+1)//2)):
  q=tournament(m+1,b);c=coeff(q)
  ck(not any(c[m//2+2:]),'context_budget')
  if b<8:
   for a,f in factors:
    if len(a)+m<=7:ck(sum(x*y for x,y in zip(c[1:],f[1:]))==brute_h(subst(q,a)),'literal_substitution')
# C3 experimental coefficients computed via independent automaton, up to j=4.
rows=[]
for j in range(1,5):
 q=[[0]*(2*j+1) for _ in range(2*j+1)]
 for u,v in it.combinations(range(2*j+1),2):
  bu=0 if u==0 else (1 if u<=j else 2);bv=0 if v==0 else (1 if v<=j else 2)
  q[u][v]=int((bu==bv and u<v) or ((bv-bu)%3==1));q[v][u]=1-q[u][v]
 c=coeff(q);rows.append(list(c[1:j+2]));ck(c[j+1]==math.factorial(j)**2,'triangular_diagonal')
 canonical=json.loads((P/'tournament/exact.json').read_text())['triangle_rows_1_to_10'][j-1][1:]
 ck(list(c[1:j+2])==canonical,'triangular_rows')
# All 512 balanced 3x3 top interfaces. Permanent computed by Ryser, not matching enumeration.
for bits in range(512):
 q=tournament(7,0)
 for i in range(1,4):q[0][i]=1;q[i][0]=0
 for i in range(4,7):q[0][i]=0;q[i][0]=1
 B=[[(bits>>(3*i+j))&1 for j in range(3)] for i in range(3)]
 for i in range(3):
  for j in range(3):q[i+1][j+4]=B[i][j];q[j+4][i+1]=1-B[i][j]
 perm=sum((-1)**(3-s.bit_count())*math.prod(sum(B[i][j] for j in range(3) if s>>j&1) for i in range(3)) for s in range(8))
 ck(coeff(q)[4]==6*perm,'r3_permanent_ryser')
# Hostile nonuniform port directly, no substitution implementation.
a=[[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]]
b=[list(x) for x in zip(*a)]
def ext(a):return [row+[int(i==0)] for i,row in enumerate(a)]+[[0,1,1,1,0]]
ck(forest_profile(a)==forest_profile(b),'converse_profile')
port=[brute_h(ext(a)),brute_h(ext(b))];ck(port==[3,15],'port_hostile')
print(json.dumps({'status':'PASS','methods':'No producer code imports. Path-extension DFS; subset Hamiltonian permutations and least-vertex partitions; quotient finite automaton; Ryser permanent.','counts':dict(counts),'fano_histogram':dict(collections.Counter(H)),'triangle_rows':rows,'nonuniform_counts':port},indent=2,sort_keys=True))
