from itertools import combinations, product, permutations
from collections import Counter
from functools import lru_cache

def shapes(k):
 out=[]
 for r in range(1,k+2):
  for starts in combinations(range(k+1),r):
   for ends in product(range(k+1),repeat=r):
    be=[x for x in ends if x<k]
    if len(set(be))!=len(be):continue
    if any(b in starts and starts[ends.index(b)]!=b for b in be):continue
    out.append(tuple(zip(starts,ends)))
 return tuple(out)

@lru_cache(None)
def transitions(k,S,drop):
 ss=shapes(k); ix={s:i for i,s in enumerate(ss)};ans=[]
 for src in ss:
  C=[(s if s<k else k+1,e if e<k else k+1) for s,e in src]
  pred=[None]+[i for i,(s,e) in enumerate(C) if e==k+1 or not(S>>e&1)]
  succ=[None]+[i for i,(s,e) in enumerate(C) if s<k and S>>s&1]
  row=Counter()
  for a,b in product(pred,succ):
   if a is not None and a==b:continue
   D=[se for i,se in enumerate(C) if i!=a and i!=b]
   D.append((C[a][0] if a is not None else k,C[b][1] if b is not None else k))
   def remap(z):
    return k if z in (drop,k+1) else z-(z>drop)
   D=tuple(sorted((remap(s),remap(e)) for s,e in D))
   if sum(s==k for s,e in D)>1:continue
   if D not in ix:raise RuntimeError((src,S,drop,D))
   row[ix[D]]+=1
  ans.append(tuple(row.items()))
 return tuple(ans)

def step(state,k,S,drop):
 out=[0]*len(state)
 for count,tr in zip(state,transitions(k,S,drop)):
  if count:
   for i,m in tr:out[i]+=m*count
 return tuple(out)

def H(state,k):return sum(c for c,s in zip(state,shapes(k)) if len(s)==1)

def brute(rows,ports):
 n=len(rows);k=len(ports);mp={v:i for i,v in enumerate(ports)};d=Counter()
 for p in permutations(range(n)):
  for cuts in range(1<<(n-1)):
   comp=[];last=0;legal=True
   for i in range(n):
    if i==n-1 or cuts>>i&1:
     comp.append((mp.get(p[last],k),mp.get(p[i],k)));last=i+1
    elif not(rows[p[i]]>>p[i+1]&1):legal=False;break
   if not legal or sum(s==k for s,e in comp)>1:continue
   # canonical component ordering by start label, not path vertices
   if [s for s,e in comp]!=sorted(s for s,e in comp):continue
   d[tuple(comp)]+=1
 return tuple(d[s] for s in shapes(k))

def physical(rows,ports,S,drop):
 n=len(rows);out=list(rows)+[0]
 for u in range(n):
  if u in ports and S>>ports.index(u)&1:out[n]|=1<<u
  else:out[u]|=1<<n
 q=ports+[n];q.pop(drop)
 return out,q

if __name__=='__main__':
 for k in (1,2,3):print(k,len(shapes(k)),sum(len(r) for S in range(1<<k) for dr in range(k+1) for r in transitions(k,S,dr)))
 seed=brute([2,0],[0,1]);print('seed',H(seed,2))
 for S in range(4):
  for drop in range(3):
   o,p=physical([2,0],[0,1],S,drop)
   if step(seed,2,S,drop)!=brute(o,p):raise RuntimeError('Mismatch')
 print('initial checks passed')
