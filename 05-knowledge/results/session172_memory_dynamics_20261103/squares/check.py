#!/usr/bin/env python3
"""Exact checks: undirected oriented Hamilton counts and one-vertex contexts."""
import itertools,json,math
from functools import lru_cache

def require(x,msg):
 if not x: raise RuntimeError(msg)

def graph(n, edges):
 a=[set() for _ in range(n)]
 for i,j in edges:a[i].add(j);a[j].add(i)
 return a

def count(a):
 n=len(a)
 @lru_cache(None)
 def f(mask,last):
  if mask==(1<<n)-1:return 1
  return sum(f(mask|1<<v,v) for v in a[last] if not mask>>v&1)
 return sum(f(1<<v,v) for v in range(n))

def response(a):
 n=len(a);E=[0]*n;K=[[0]*n for _ in range(n)];H=0
 for p in itertools.permutations(range(n)):
  bad=[t for t in range(n-1) if p[t+1] not in a[p[t]]]
  if not bad:H+=1;E[p[-1]]+=1
  cuts=range(n-1) if not bad else bad if len(bad)==1 else []
  for t in cuts:K[p[t]][p[t+1]]+=1
 return H,E,K

def extend(a,S):
 b=[set(x) for x in a]+[set(S)];n=len(a)
 for v in S:b[v].add(n)
 return b

def square(n):return graph(n,((i,j) for i in range(n) for j in range(i+1,n) if math.isqrt(i+j+2)**2==i+j+2))

def paths(a):
 n=len(a);out=[]
 def rec(p,mask):
  if len(p)==n:out.append(p);return
  for v in sorted(a[p[-1]]):
   if not mask>>v&1:rec(p+(v,),mask|1<<v)
 for v in range(n):rec((v,),1<<v)
 return out

def join_count(a,i,j):
 # Independent cut-cover enumeration by subset path DP, not by extended paths.
 n=len(a)
 @lru_cache(None)
 def ends(mask,last):
  if mask==1<<last:return 1
  rest=mask^(1<<last)
  return sum(ends(rest,v) for v in a[last] if rest>>v&1)
 # Enumerate all simple paths ending i and j, aggregate by their exact sets.
 def banks(root):
  bank={}
  def rec(v,mask):
   bank[mask]=bank.get(mask,0)+1
   for w in a[v]:
    if not mask>>w&1:rec(w,mask|1<<w)
  rec(root,1<<root)
  return bank
 bi,bj=banks(i),banks(j);full=(1<<n)-1
 return sum(v*bj.get(full^s,0) for s,v in bi.items())

rank_rows={}
exhaust=[]
for n in range(1,6):
 edges=list(itertools.combinations(range(n),2));contexts=0
 for bits in range(1<<len(edges)):
  a=graph(n,[e for k,e in enumerate(edges) if bits>>k&1]);H,E,K=response(a)
  require(H==count(a),'base')
  rank_rows.setdefault(n,[]).append(E+[K[i][j] for i,j in itertools.combinations(range(n),2)])
  require(all(K[i][j]==K[j][i] for i in range(n) for j in range(n)),'symmetric')
  for sm in range(1<<n):
   S=[i for i in range(n) if sm>>i&1]
   predicted=2*sum(E[i] for i in S)+sum(K[i][j] for i in S for j in S)
   require(predicted==count(extend(a,S)),f'context {n} {bits} {sm}')
   contexts+=1
 exhaust.append({'n':n,'graphs':1<<len(edges),'contexts':contexts})
records=[]
for n in (15,16,17,18,22,23,24,25):
 a=square(n);p=paths(a);require(len(p)==count(a),'square independent DP')
 records.append({'n':n,'oriented_H':len(p),'sample_paths':[[v+1 for v in q] for q in p[:1]]})
resurrection=[]
for n,i,j in [(22,2,13),(24,11,24)]:
 a=square(n);b=square(n+1)
 require(b[-1]=={i-1,j-1},'actual arithmetic neighbor set')
 require(count(a)==0,'zero base')
 K=join_count(a,i-1,j-1)
 require(2*K==count(b),'join certificate')
 p=paths(b);covers=[]
 for q in p:
  t=q.index(n)
  require(0<t<n,'new vertex must be internal')
  if q[t-1]==i-1:
   covers.append([[v+1 for v in q[:t]],[v+1 for v in q[t+1:]]])
 require(len(covers)==K,'cover list')
 resurrection.append({'base_n':n,'new_ports':[i,j],'K':K,'oriented_next_H':2*K,'covers':covers})
from fractions import Fraction

def rank(rows):
 pivots={}
 for row in rows:
  row=list(map(Fraction,row))
  for i in sorted(pivots):
   if row[i]:
    q=row[i];row=[x-q*y for x,y in zip(row,pivots[i])]
  if any(row):
   i=next(i for i,x in enumerate(row) if x);q=row[i]
   pivots[i]=[x/q for x in row]
 return len(pivots)
ranks={str(n):rank(rows) for n,rows in rank_rows.items()}
for n in range(3,6):require(ranks[str(n)]==n*(n+1)//2,'full response rank')

def profile(a):
 n=len(a);hp={}
 for mask in range(1,1<<n):
  vs=[i for i in range(n) if mask>>i&1]
  hp[mask]=sum(all(q[i+1] in a[q[i]] for i in range(len(q)-1)) for q in itertools.permutations(vs))
 @lru_cache(None)
 def f(mask):
  if not mask:return (1,)+(0,)*n
  out=[0]*(n+1);low=mask&-mask;s=mask
  while s:
   if s&low:
    for k,x in enumerate(f(mask^s)[:-1]):out[k+1]+=hp[s]*x
   s=(s-1)&mask
  return tuple(out)
 return f((1<<n)-1)
seed_edges=[[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4)],[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(2,3)]]
seeds=[graph(5,e) for e in seed_edges]
profiles=[profile(a) for a in seeds]
require(profiles[0]==profiles[1]==(0,12,60,54,14,1),'full profile collision')
require(sorted(map(len,seeds[0]))!=sorted(map(len,seeds[1])),'nonisomorphism')
reads=[count(extend(a,[1])) for a in seeds]
require(reads==[0,4],'read head')
empty5=graph(5,[]);three_pieces=graph(5,[(0,1),(2,3)])
require(response(empty5)==response(three_pieces),'equal one-context memory')
def two_new(a):
 b=extend(a,range(len(a)))
 return extend(b,range(len(a)))
two_reads=[count(two_new(a)) for a in (empty5,three_pieces)]
require(two_reads==[0,48],'two-context separation')
print(json.dumps({'convention':'oriented Hamilton paths; reverse orders counted separately','exhaustive':exhaust,'response_ranks':ranks,'one_context_equal_two_context_different':two_reads,'equal_profile_seeds':{'edges':seed_edges,'profile':profiles[0],'leaf_at_label_1_responses':reads},'square_prefixes':records,'resurrections':resurrection},sort_keys=True,indent=2))
