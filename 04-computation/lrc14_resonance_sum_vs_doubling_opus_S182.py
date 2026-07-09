"""
lrc14_resonance_sum_vs_doubling_opus_S182.py  (opus-2026-07-09-S182)

EXTENDING THE FORWARD LEAD (BSG -> Freiman 3k-4 -> density floor). The lonely measure decomposes as
L(S) = (6/7)^13 + R(S), R(S) = Sum_{n: n.v=0} What(n) = a theta-sum over the resonance lattice Lambda.
The density-floor inequality the whole covering case needs is |R| < (6/7)^13 (strict, with equality only
at the tight AP). S181 says looseness is governed by the resonance-lattice GEOMETRY, not the scalar E.
This asks the sharp a-priori question: is |R| (equivalently the tightness deficit (6/7)^13 - L) governed
MONOTONICALLY by the DOUBLING K=|S+S|/13 -- so that K > K_min => |R| < (6/7)^13 with a margin that is a
function of the doubling? If yes, the density floor follows from a doubling bound (Freiman/Plunnecke),
a-priori and uniform, sidestepping the empirical census.
"""
import numpy as np, random
from collections import Counter

NG = 300007
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14
MAIN = (6.0/7.0)**13

def sumset(S): return len(set(a+b for a in S for b in S))
def Edd(S):
    r = Counter(a+b for a in S for b in S); return sum(c*c for c in r.values())
def longest_ap(S):
    Sset=set(S); S=sorted(S); best=1
    for i,a in enumerate(S):
        for b in S[i+1:]:
            d=b-a
            if a-d in Sset: continue
            L=2; x=b+d
            while x in Sset: L+=1; x+=d
            best=max(best,L)
    return best
def lonely(S):
    M=np.zeros(NG)
    for v in S:
        d=np.abs(((v*TAU+0.5)%1.0)-0.5); M+=(d<=H)
    return float((M==0).mean())

def gap(dims, steps):
    import itertools
    return sorted(set(1+sum(x*s for x,s in zip(c,steps)) for c in itertools.product(*[range(d) for d in dims])))

fams=[]
fams.append(("AP {1..13}", list(range(1,14))))
fams.append(("2*AP {2..26}", [2*i for i in range(1,14)]))
fams.append(("near-AP {1..12}+20", list(range(1,13))+[20]))
fams.append(("near-AP {1..11}+{20,30}", list(range(1,12))+[20,30]))
fams.append(("near-AP {1..10}+{20,30,40}", list(range(1,11))+[20,30,40]))
fams.append(("GAP 4x4 P=5", gap([4,4],[1,5])[:13]))
fams.append(("GAP 4x4 P=7", gap([4,4],[1,7])[:13]))
fams.append(("GAP 7x2 P=9", gap([7,2],[1,9])[:13]))
fams.append(("GAP 3x5 P=11", gap([5,3],[1,11])[:13]))
fams.append(("GAP3D 2x2x4", gap([4,2,2],[1,5,20])[:13]))
fams.append(("Sidon-ish", [1,2,5,11,22,33,45,60,78,95,110,130,150]))
fams.append(("Fib-ish", [1,2,3,5,8,13,21,34,55,89,144,233,377]))
# a few random dissociated
rng=random.Random(182)
for s in range(3):
    r=random.Random(500+s); sp=r.randint(40,120)
    S=sorted(set([1]+r.sample(range(2,sp),11)+[sp]))
    while len(S)!=13: S=sorted(set(S+[r.randint(1,sp)]))[:13]
    fams.append((f"random sp{sp}", S))

print("="*104)
print(f"RESONANCE SUM |R|=|L-(6/7)^13| vs DOUBLING K.  (6/7)^13={MAIN:.5f}  [tight: R=-{MAIN:.3f}, L=0]")
print("="*104)
print(f"  {'set':>28} {'K=|S+S|/13':>10} {'lAP':>4} {'E':>5} {'L':>8} {'R=L-main':>9} {'|R|/main':>8} {'deficit':>8}")
rows=[]
for name,S in fams:
    S=sorted(set(S))
    if len(S)!=13: continue
    K=sumset(S)/13; L=lonely(S); R=L-MAIN; deficit=MAIN-L
    rows.append((K,name,longest_ap(S),Edd(S),L,R,abs(R)/MAIN,deficit))
rows.sort()
for K,name,lap,E,L,R,absRr,deficit in rows:
    print(f"  {name:>28} {K:>10.3f} {lap:>4} {E:>5} {L:>8.4f} {R:>9.4f} {absRr:>8.3f} {deficit:>8.4f}")
print("-"*104)
print("  If |R| < main STRICTLY for all K>K_min (only AP at K_min has |R|=main, L=0), and the deficit")
print("  (main-L) grows with K, the density floor |R|<(6/7)^13 is a DOUBLING bound: K>K_min => loose.")
print("="*104)
