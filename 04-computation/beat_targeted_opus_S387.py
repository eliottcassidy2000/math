# opus-2026-07-17-S387 -- HYP-7730: CAN THE BEAT CERTIFICATE BE DEFEATED ON PURPOSE?
#
# S384 proved the pairing lemma and left the certificate as a conjecture backed
# by blind hill-climbing -- weak evidence, given MISTAKE-152/154/156/157.  Here
# is a sharper handle and a TARGETED attack.
#
# CHARACTERISATION.  W_q is symmetric, so killset(v,q) depends only on the class
# {+-v mod q}.  Hence
#     #distinct blockers at q  =  |{ +-v mod q : v in V }|
# and blocking q requires  |{+-v mod q}| * k_q  >=  phi(q).
# The certificate fires at q iff  |{+-v mod q}| * k_q < phi(q),  i.e. iff q has a
# high ratio phi(q)/k_q AND the speeds collide enough mod q.
#
# THE ATTACK.  phi(q)/k_q is maximised at q = 90 (=12) and is SMALL for highly
# composite q (phi small) -- so a family whose beats all land on low-ratio q
# would defeat the certificate.  Build families with highly composite beats.
from math import gcd
from itertools import combinations
import random
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def killset(v,q,Wq): return frozenset(p for p in range(1,q) if gcd(p,q)==1 and (v*p)%q in Wq)
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def kq(q,Wq): return sum(1 for w in Wq if w!=0 and gcd(w,q)==1)
def pmclasses(V,q): return len({min(v%q,(-v)%q) for v in V})

print("(1) VERIFY: #distinct kill-sets == |{+-v mod q}| ?")
random.seed(387)
bad=0; n=0
for _ in range(120):
    q=random.randint(15,140); Wq=W(q)
    V=sorted(random.sample(range(1,300),13))
    a=len({killset(v,q,Wq) for v in V}); b=pmclasses(V,q)
    n+=1
    if a!=b: bad+=1
print(f"    {n} (V,q) pairs; mismatches = {bad}")
print("    (0 => blocking capacity is exactly the number of +- classes mod q)")

print()
print("(2) THE RATIO LANDSCAPE phi(q)/k_q -- which q can certify at all?")
print("    certificate at q needs |{+-v}| < phi(q)/k_q, and |{+-v}| >= 1")
rows=[]
for q in range(15,200):
    Wq=W(q); k=kq(q,Wq)
    if k==0: continue
    rows.append((phi(q)/k, q, phi(q), k))
rows.sort(reverse=True)
print("    best q:", [(q,f"{r:.1f}") for r,q,_,_ in rows[:6]])
print("    worst q:", [(q,f"{r:.1f}") for r,q,_,_ in rows[-6:]])
lowratio={q for r,q,_,_ in rows if r<=4}
print(f"    q with ratio <= 4 (hard to certify): {len(lowratio)} of {len(rows)}")

print()
print("(3) TARGETED ATTACK: build families whose beats are all LOW-ratio q")
print("    (highly composite beats => small phi => small phi/k => no certificate)")
def beats(V):
    B=set()
    for a,b in combinations(V,2):
        B.add(a+b); d=abs(a-b)
        if d>0: B.add(d)
    for v in V: B.add(2*v)
    return {q for q in B if q>=2}
def certifies(V,cap=400):
    best=None
    for q in sorted(beats(V)):
        if q>cap: continue
        if q<=14:
            if not any(v%q==0 for v in V): return (q,'sieve',phi(q))
            continue
        if any(v%q==0 for v in V): continue
        Wq=W(q); k=kq(q,Wq)
        if k==0: continue
        m=phi(q)-pmclasses(V,q)*k
        if best is None or m>best[2]: best=(q,'count',m)
    return best
# candidates: all speeds multiples of a highly composite base, plus offsets that
# keep sums/differences highly composite
worst=(10**9,None,None)
random.seed(3871)
for base in [210, 180, 120, 90, 60]:
    for trial in range(40):
        V=sorted({base*random.randint(1,6)+random.choice([0,0,0,base//2]) for _ in range(20)})[:13]
        if len(V)!=13: continue
        c=certifies(V); m=c[2] if c else -10**9
        if m<worst[0]: worst=(m,list(V),c)
print(f"    lowest certificate margin from composite-beat construction: {worst[0]}")
print(f"      V = {worst[1]}")
print(f"      best beat = {worst[2]}")
