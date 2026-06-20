import sys
sys.path.insert(0, 'C:/Users/Eliott/Documents/GitHub/math/04-computation')
from fractions import Fraction
import itertools, random
from functools import reduce
from math import gcd
import importlib.util
spec=importlib.util.spec_from_file_location("m","C:/Users/Eliott/Documents/GitHub/math/04-computation/lrc14_tower_deletion_proof_kps.py")
# just re-define measure here to avoid running the print module
def lonely_measure(C, theta=Fraction(1,14)):
    arcs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0, d+1):
            arcs.append((Fraction(m,d)-w, Fraction(m,d)+w))
    segs=[]
    for lo,hi in arcs:
        for shift in (-1,0,1):
            a2=max(lo+shift,Fraction(0)); b2=min(hi+shift,Fraction(1))
            if a2<b2: segs.append((a2,b2))
    segs.sort()
    cur=Fraction(-1); union=[]
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in union)

thr2=Fraction(426,35035)
# Direct random check of the FULL claim on actual AP-tail 12-cores missing 2^a, with tails:
random.seed(1)
bad=0; tested=0
for a in [0,1,2,3]:
    bit=2**a
    for trial in range(4000):
        nh=random.choice([1,2,3,4])
        others=[d for d in range(1,14) if d!=bit]
        extra=random.sample(others, nh-1)
        holes=set(extra)|{bit}
        ntails=nh-1
        tails=tuple(sorted(random.sample(range(14,60), ntails))) if ntails>0 else ()
        C=tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))
        if len(C)!=12: continue
        if reduce(gcd,C)!=1: continue
        tested+=1
        L=lonely_measure(C)
        if L<thr2:
            bad+=1
            print("COUNTEREX", a, holes, tails, L, float(L))
print(f"random full-core test: tested={tested}, below_thr2={bad}")
