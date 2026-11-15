from functools import lru_cache
def append(r,p):
 n=len(r);return tuple(a | ((i!=p)<<n) for i,a in enumerate(r))+(1<<p,)
@lru_cache(None)
def F(r):
 n=len(r);b=1<<n;dp=[[0]*n for _ in range(b)];h=[1]+[0]*(b-1)
 for s in range(1,b):
  w=s
  while w:
   z=w&-w;i=z.bit_length()-1;w^=z;t=s^z
   if not t:dp[s][i]=1;continue
   u=t&r[i]
   while u:
    v=u&-u;u^=v;dp[s][i]+=dp[t][v.bit_length()-1]
  h[s]=sum(dp[s])
 d=[-h[-1]]*n
 for s in range(1,b):
  if s.bit_count()%2:continue
  t=2*h[-1^s]
  for i,c in enumerate(dp[s]):d[i]+=t*c
 return h[-1],d
def generator(n, chronological):
 if chronological:
  rr=[((0,),())]
  for i in range(1,n):
   rr=[(append(r,q),w+(q,)) for r,w in rr for q in range(i)]
  yield from rr
 else:
  from itertools import combinations
  edges=list(combinations(range(n),2))
  for code in range(1<<len(edges)):
   r=[0]*n
   for i,(a,b) in enumerate(edges):
    if (code>>i)&1:r[a]|=1<<b
    else:r[b]|=1<<a
   yield tuple(r),code

def mixed(r):
 _,d=F(r)
 return -1 in d and 1 in d

def census(max_n, chronological):
 from collections import Counter
 from itertools import permutations
 from json import dumps
 def sccs(r):
  m=[a|(1<<i) for i,a in enumerate(r)]
  for i in range(len(r)):
   for j in range(len(r)):
    if (m[j]>>i)&1:m[j]|=m[i]
  c={}
  for q,u in enumerate(m):c.setdefault(u,[]).append(q)
  return [tuple(c[x]) for x in sorted(c,key=int.bit_count,reverse=True)]
 def canon(r,c):
  return min(tuple(sum(((r[i]>>j)&1)<<k for k,j in enumerate(b)) for i in b) for b in permutations(c))
 result={}
 for n in range(1,max_n+1):
  counts=Counter()
  cores=Counter()
  signs=Counter()
  tested=0
  for r,w in generator(n,chronological):
   tested+=1
   h,d=F(r)
   if -1 not in d or 1 not in d: continue
   counts[h]+=1
   cs=sccs(r)
   if max(map(len,cs))<=8:
    key=str([canon(r,c) for c in cs])
    cores[key]+=1
   a,b=(any(mixed(append(r,q)) for q in range(n) if d[q]==e) for e in (-1,1))
   signs[str((a,b))]+=1
  result[n]={'sources':tested,'mixed_by_H':dict(counts),'mixed_by_SCC':dict(cores),'child_mixed':dict(signs)}
  F.cache_clear()
 return result

if __name__=='__main__':
 from argparse import ArgumentParser
 from json import dumps
 parser=ArgumentParser()
 parser.add_argument('--chronological',type=int,default=8)
 parser.add_argument('--generic',type=int,default=6)
 a=parser.parse_args()
 print(dumps({'chronological':census(a.chronological,True),'generic':census(a.generic,False)},indent=2))
