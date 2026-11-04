import itertools,json,math
from query import hp,extend
from verify import SEEDS
A=(4,1,2,7,15);B=(10,0,3,6,15)
def covers(a,k):
 n=len(a);h={0:1}
 for mask in range(1,1<<n):
  ids=[i for i in range(n) if mask>>i&1]
  h[mask]=hp([sum(1<<j for j,v in enumerate(ids) if a[u]>>v&1) for u in ids])
 out=0
 for colors in itertools.product(range(k),repeat=n):
  sets=[sum(1<<i for i,c in enumerate(colors) if c==j) for j in range(k)]
  if all(sets):out+=math.prod(h[s] for s in sets)
 return out
def pairs(a):return [list((u,v)) for u,v in itertools.combinations(range(len(a)),2) if (a[u]^a[v])&~((1<<u)|(1<<v))==0]
expected=[3,42,162,240,120]
profiles=[[covers(t,k) for k in range(1,6)] for t in (A,B)]
if profiles!=[expected,expected]:raise RuntimeError('profile')
if pairs(A)!=[[3,4]] or pairs(B)!=[]:raise RuntimeError('twins')
values=[hp(extend(t,0,1)) for t in (A,B)]
if values!=[15,9]:raise RuntimeError('full query')
def substitute(seed,module):
 n=len(module);out=list(module)+[0]*5
 for u in range(5):
  for v in range(5):
   if seed[u]>>v&1:out[n+u]|=1<<(n+v)
  for v in range(n):
   if seed[u]>>5&1:out[n+u]|=1<<v
   else:out[v]|=1<<(n+u)
 return out
ambient=[]
for j,seed in enumerate(SEEDS):
 base=[hp(substitute(seed,t)) for t in (A,B)]
 queries=[hp(substitute(seed,extend(t,0,1))) for t in (A,B)]
 if base!=[1611+6*j]*2 or queries!=[21015+30*j,17793+18*j]:raise RuntimeError('ambient')
 ambient.append({'state':j,'H':base,'port_values':queries})
print(json.dumps({'modules':[A,B],'profiles':profiles,'homogeneous_pairs':[pairs(A),pairs(B)],'query':{'a':0,'c':1,'values':values},'strong_contexts':ambient,'status':'PASS'},indent=2,sort_keys=True))
