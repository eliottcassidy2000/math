"""Independent W9 subset recurrence and backward equitable signatures.
No producer imports; checks survive python -O.
"""
from collections import Counter
import json

def require(x,msg):
    if not x: raise RuntimeError(msg)

def arc(a,b):
    i,j=divmod(a,3); k,l=divmod(b,3)
    return ((k-i)%3==1) if i!=k else (l-j)%3==1

def tau(a):
    i,j=divmod(a,3); return 3*((-i)%3)+(-j)%3
pairs=[]
for v in range(1,9):
    if v<tau(v): pairs.append((v,tau(v)))
owner={v:i for i,p in enumerate(pairs) for v in p}
states=[(m,v) for m in range(1,16) for v in range(1,9) if m>>owner[v]&1]
out={s:[(s[0]|(1<<i),v) for i,p in enumerate(pairs) if not s[0]>>i&1 for v in p if arc(s[1],v)] for s in states}
require(len(states)==64,'64 states')
def quotient(named=False):
    ids={}; signatures={}; labels={}
    for level in range(4,0,-1):
        for s in states:
            if s[0].bit_count()!=level:continue
            sig=(level,int(s[0]==15 and arc(s[1],0)),int(s[1]==3) if named else None,tuple(sorted(Counter(ids[t] for t in out[s]).items())))
            if sig not in labels: labels[sig]=len(labels)
            ids[s]=labels[sig];signatures[s]=sig
    # Since acyclic graded, these unique recursive signatures equal coarsest
    # partition retaining level, output, and outgoing multiplicities.
    for a in states:
        for b in states:
            if ids[a]==ids[b]:
                require(Counter(ids[t] for t in out[a])==Counter(ids[t] for t in out[b]),'equitable')
    return ids
q=quotient();r=quotient(True)
require(len(set(q.values()))==16,'16 quotient');require(len(set(r.values()))==35,'35 quotient')
# Unfolded ordinary Hamilton DP independently gives H(W9).
dp={(1<<v,v):1 for v in range(9)}
for m in range(1,512):
    for v in range(9):
        c=dp.get((m,v),0)
        if not c:continue
        for w in range(9):
            if not m>>w&1 and arc(v,w):dp[(m|(1<<w),w)]=dp.get((m|(1<<w),w),0)+c
h=sum(dp.get((511,v),0) for v in range(9));require(h==3159,'H W9')
f={s:int(s[0]==15 and arc(s[1],0)) for s in states}
for level in (3,2,1):
    for s in states:
        if s[0].bit_count()==level:f[s]=sum(f[t] for t in out[s])
hfold=sum(f[(1<<owner[v],v)] for v in range(1,9));require(hfold==33,'folded')
print(json.dumps({'states':64,'arcs':sum(map(len,out.values())),'quotient':len(set(q.values())),'named_vertex_3':len(set(r.values())),'H_W9':h,'folded':hfold,'levels16':dict(Counter(m.bit_count() for m,v in [next(s for s in states if q[s]==c) for c in set(q.values())]))},sort_keys=True))
# Separate exact linear-algebra implementation (SymPy rational elimination).
from sympy import Matrix
base={s:(s[0].bit_count(),int(s[0]==15 and arc(s[1],0))) for s in states}
vectors=[]
for c in sorted(set(base.values())):
    f={s:int(base[s]==c) for s in states}
    for k in range(4):
        vectors.append([f[s] for s in states]);f={s:sum(f[t] for t in out[s]) for s in states}
mat=Matrix(vectors).T
rk=mat.rank();require(rk==12,'linear rank12')
basis=[list(mat[:,j]) for j in mat.rref()[1]]
products=[[x*y for x,y in zip(a,b)] for i,a in enumerate(basis) for b in basis[i:]]
product_rank=Matrix(basis+products).rank()
print(json.dumps({'linear_rank':rk,'quadratic_product_rank':product_rank}))
# Named observer linear signatures do not even induce the stable partition.
colours={s:(base[s],s[1]==3) for s in states};cols=[]
for c in sorted(set(colours.values())):
 f={s:int(colours[s]==c) for s in states}
 for _ in range(4):
  cols.append([f[s] for s in states]);f={s:sum(f[t] for t in out[s]) for s in states}
require(Matrix(cols).rank()==21,'named rank21')
sigs=list(zip(*cols));require(len(set(sigs))==34,'named signatures34')
x,y=(8,7),(4,8)
require(sigs[states.index(x)]==sigs[states.index(y)],'linear blindness')
f={s:int(s[0].bit_count()==4 and not arc(s[1],0) and s[1]!=3) for s in states}
for _ in range(2):f={s:sum(f[t] for t in out[s]) for s in states}
z={s:sum(f[t]**2 for t in out[s]) for s in states}
require((z[x],z[y])==(3,9),'nonlinear separator')
print(json.dumps({'named_linear_rank':21,'named_linear_signature_classes':34,'nonlinear_separator':[z[x],z[y]]}))
menus={}
for mask in range(256):
 ids={};labels={}
 for level in range(4,0,-1):
  for s in states:
   if s[0].bit_count()!=level:continue
   observation=s[1] if mask>>(s[1]-1)&1 else 0
   sig=(base[s],observation,tuple(sorted(Counter(ids[t] for t in out[s]).items())))
   ids[s]=labels.setdefault(sig,len(labels))
 menus[mask]=len(labels)
minimal=[m for m in menus if menus[m]==64 and all(menus[m^(1<<i)]<64 for i in range(8) if m>>i&1)]
require(minimal==[126,190,222,237,245,249],'minimal menus')
require(menus[1]+menus[4]<menus[0]+menus[5],'non submodular')
require(menus[1]+menus[2]>menus[0]+menus[3],'non supermodular')
print(json.dumps({'menus':256,'minimal_full_menus':minimal}))
