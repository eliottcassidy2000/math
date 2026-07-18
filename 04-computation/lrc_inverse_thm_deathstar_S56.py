from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
random.seed(5)
def covers(fam): return all(any(v%q==0 for v in fam) for q in range(2,15))
def M_lt_13(fam):   # True iff M(fam) < 1/13  (bail as soon as a time reaches >=1/13)
    Q=2*max(fam)+2
    for q in range(2,Q+1):
        thr=q  # need 13*m >= q for M>=1/13
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 13*m>=q: return False   # M>=1/13
    return True
def is_dilated_AP(W):
    W=sorted(W); d=W[0]
    return W==[d*i for i in range(1,13)]
found=[]; tested=0; mlt=0
for trial in range(200000):
    if random.random()<0.6:
        W=list(range(1,13))
        for _ in range(random.randint(1,3)):
            i=random.randrange(12); W[i]=random.randint(1,45)
        W=sorted(set(W))
        if len(W)!=12: continue
        V=sorted(set(W+[random.choice([182,169,170,26,39,52,24,25,84,183])]))
    else:
        V=sorted(random.sample(range(1,55),13))
    if len(V)!=13 or not covers(V) or reduce(gcd,V)!=1: continue
    tested+=1
    if M_lt_13(V):
        mlt+=1
        W=V[:-1]
        if not is_dilated_AP(W):
            found.append((V,W))
            if len(found)<=8: print(f"  NON-AP CORE with M<1/13: V={V} W={W}")
print(f"tested {tested} covering, {mlt} had M<1/13, non-AP-core among them: {len(found)}")
if not found and mlt>0: print("  => inverse theorem HOLDS on sample (every M<1/13 covering family has dilated-AP core)")
