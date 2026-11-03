from itertools import combinations, permutations
from random import Random
R=Random(3)
n=7
pairs=list(combinations(range(n),2)); E=list(combinations(range(n),2)); cycs=[]
for k in [3,5,7]:
 for s in combinations(range(n),k):
  for p in permutations(s[1:]): cycs.append((s[0],)+p)
def orient(mask,u,v):
 i=E.index((min(u,v),max(u,v)))
 return ((mask>>i)&1)==(u<v)
def cycles(mask):
 return [c for c in cycs if all(orient(mask,c[i],c[(i+1)%len(c)]) for i in range(len(c)))]
def h_sub(cs,mask):
 tris=[set(c) for c in cs if len(c)==3 and all((mask>>v)&1 for v in c)]
 one=sum(all((mask>>v)&1 for v in c) for c in cs)
 two=sum(not a&b for a,b in combinations(tris,2))
 return 1+2*one+4*two
for z in range(1000):
 mask=R.randrange(1<<len(E)); old=cycles(mask)
 for ix,e in enumerate(E):
  new=cycles(mask^(1<<ix)); so=set(old); sn=set(new); add=sn-so; rem=so-sn
  if add and rem:
   terms=[]
   for sign,coll in [(1,add),(-1,rem)]:
    for c in coll:
     comp=((1<<n)-1)^sum(1<<i for i in c)
     terms.append((sign,c,h_sub(old,comp)))
   raw=len(add)-len(rem); weighted=sum(s*w for s,c,w in terms)
   if raw!=weighted and raw*weighted<=0:
    print('mask',mask,'edge',e,'old',h_sub(old,(1<<n)-1),'new',h_sub(new,(1<<n)-1),'raw',raw,'weighted',weighted,'terms',sorted(terms)); print('arcs',[p if mask>>i&1 else p[::-1] for i,p in enumerate(E)]); raise SystemExit
