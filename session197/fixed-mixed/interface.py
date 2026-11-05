from functools import lru_cache
from itertools import permutations,product
from collections import Counter
@lru_cache(None)
def shapes(k):
 types=[(i,j)for i in range(k+1)for j in range(k+1)if(i,j)!=(k,k)];ss=set()
 def rec(at,C,used,P):
  for q in range(max(0,1-len(C)),P+2-len(C)):ss.add(tuple(sorted(C+[(k,k)]*q)))
  for a in range(at,len(types)):
   x,y=types[a];ports={z for z in(x,y)if z<k}
   if ports&used:continue
   rec(a+1,C+[(x,y)],used|ports,P+(x<k)+(y<k))
 rec(0,[],set(),0);return tuple(sorted(ss,key=lambda s:(len(s),s)))
@lru_cache(None)
def trans(k,side,S):
 ss=shapes(k);ix={s:i for i,s in enumerate(ss)};ans=[]
 for C in ss:
  if side=='R':
   pred=[None]+[i for i,(s,e)in enumerate(C)if e==k or not S>>e&1]
   succ=[None]+[i for i,(s,e)in enumerate(C)if s<k and S>>s&1]
  else:
   pred=[None]+[i for i,(s,e)in enumerate(C)if e<k and S>>e&1]
   succ=[None]+[i for i,(s,e)in enumerate(C)if s==k or not S>>s&1]
  counts=Counter()
  for a,b in product(pred,succ):
   if a is not None and a==b:continue
   D=[x for i,x in enumerate(C)if i!=a and i!=b]
   D.append((C[a][0]if a is not None else k,C[b][1]if b is not None else k));D=tuple(sorted(D))
   P=sum(x<k for c in D for x in c)
   if len(D)-P>1:continue
   if D not in ix:raise RuntimeError((k,side,S,C,a,b,D))
   counts[ix[D]]+=1
  ans.append(tuple(counts.items()))
 return tuple(ans)
def step(v,k,side,S):
 out=[0]*len(v)
 for x,row in zip(v,trans(k,side,S)):
  if x:
   for j,m in row:out[j]+=x*m
 return tuple(out)
def H(v,k):return sum(x for x,s in zip(v,shapes(k))if len(s)==1)
def brute(rows,ports):
 k=len(ports);n=len(rows);mp={v:i for i,v in enumerate(ports)};counts=Counter()
 for p in permutations(range(n)):
  for mask in range(1<<(n-1)):
   comps=[];last=0;ok=True
   for i in range(n):
    if i==n-1 or mask>>i&1:comps.append(p[last:i+1]);last=i+1
    elif not rows[p[i]]>>p[i+1]&1:ok=False;break
   if not ok or comps!=sorted(comps):continue
   s=tuple(sorted((mp.get(c[0],k),mp.get(c[-1],k))for c in comps))
   P=sum(x<k for c in s for x in c)
   if len(s)-P<=1:counts[s]+=1
 return tuple(counts[s]for s in shapes(k))
def physical(rows,ports,side,S):
 n=len(rows);special={p for i,p in enumerate(ports)if S>>i&1};out=rows[:]+[0]
 for v in range(n):
  win=(v in special)if side=='R'else(v not in special)
  if win:out[n]|=1<<v
  else:out[v]|=1<<n
 return out
