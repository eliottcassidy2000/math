import json
from collections import deque
from math import gcd
M={'O':((1,0,0,0),(0,2,0,0),(0,1,1,0),(1,2,0,0)), 'Z':((0,0,0,1),(0,1,1,0),(0,0,2,0),(0,0,2,1)), 'J':((0,0,0,1),(1,1,0,0),(0,1,0,0),(1,2,0,0))}
def apply(A,x): return tuple(sum(a*b for a,b in zip(row,x)) for row in A)
def supp(x): return sum((a>0)<<i for i,a in enumerate(x))
def bitvec(s): return tuple((s>>i)&1 for i in range(4))
def step(s,a): return supp(apply(M[a],bitvec(s)))
# Exact DFA equivalence, all sixteen starting support states.
blocks=[{s for s in range(16) if not s&1},{s for s in range(16) if s&1}]
while True:
    where={s:i for i,b in enumerate(blocks) for s in b}; new=[]
    for b in blocks:
        bins={}
        for s in b: bins.setdefault(tuple(where[step(s,a)] for a in M),set()).add(s)
        new.extend(bins.values())
    if len(new)==len(blocks): break
    blocks=new
# Quantitative vectors independently replay support through length 8.
front=[(s,bitvec(s)) for s in range(16)]; cases=0
for depth in range(9):
    for s,x in front: assert supp(x)==s; cases+=1
    if depth<8: front=[(step(s,a),apply(M[a],x)) for s,x in front for a in M]
# product compiler demonstration: actual 401 repair difference (0,12,12,12).
# For H-neutral words only O survives; require cheaper state U=573*2^k ==1 mod m.
def shortest(m):
    start=(14,(401%m,573%m,532%m,933%m)); todo=deque([(start,'')]); seen={start}
    while todo:
        (s,x),w=todo.popleft()
        if not s&1 and x[1]==1%m:return w,len(seen),x
        for a,A in M.items():
            y=(step(s,a),tuple(v%m for v in apply(A,x)))
            if y not in seen:seen.add(y);todo.append((y,w+a))
    return None,len(seen),None
# Only m<=11 to keep product compact. Nonunit odd moduli hostile included.
products={str(m):shortest(m) for m in [3,5,7,9,11]}
# Fixed signed-affine words: exact image and fibre counts, no orbit claim.
fibres=[]
for m in range(3,50,2):
 for r in range(5):
  for K in [r,2*r+1]:
   a=pow(3,r,m)*pow(pow(2,K,m),-1,m)%m;b=(r*r+1)%m
   sizes={}
   for x in range(m):sizes[(a*x+b)%m]=sizes.get((a*x+b)%m,0)+1
   d=gcd(3**r,m)
   assert len(sizes)==m//d and set(sizes.values())=={d}
   # Kernel-index sidecar via x = j*(m/d)+base.
   seen={( (a*x+b)%m, x//(m//d)) for x in range(m)}
   assert len(seen)==m
   fibres.append((m,r,K,d))
print(json.dumps({'support_cases':cases,'minimal_dfa_classes':[sorted(b) for b in blocks], 'support_transitions':{str(s):{a:step(s,a) for a in M} for s in range(16)},'product_examples':products,'affine_fibre_cases':len(fibres)},indent=2))
