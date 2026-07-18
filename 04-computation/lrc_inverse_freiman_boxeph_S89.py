#!/usr/bin/env python3
"""
Working the inverse-additive statement creatively  (boxeph-2026-07-18-S89)
=========================================================================
INV: M<1/13 covering => V minus v_max is a dilated AP.  Reduce/relate to known ideas:

 (F) FREIMAN 3k-4: a 12-set A with |A-A| <= 3*12-4 = 32 lies in an AP of length
     <= |A-A|-12+1; and |A-A|>=2*12-1=23 with EQUALITY iff A is a 12-term AP.
     So INV <=> 'M<1/13 => the 12-core has |core-core| = 23 (minimal)'.
 (S) THE M-SPECTRUM in (1/14,1/13): is it discrete? exactly {14m/(182m+1)}?
 (E) ADDITIVE ENERGY of V vs M (Balog-Szemeredi-Gowers / THM-730 AP-maximizes-triples).
 (O) OSTROWSKI/continued fraction of the maximizer.
We compute these across M<1/13 families and near-tight families and look for the
clean threshold that would let an inverse theorem fire.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def exact_M(V):
    best=F(0); qs=set([13,14])
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for dd in range(1,s+1):
                if s%dd==0: qs.add(dd)
    for q in qs:
        for a in range(1,q):
            if gcd(a,q)==1:
                m=min(min((v*a)%q,q-(v*a)%q) for v in V); c=F(m,q)
                if c>best: best=c
    return best
def cov(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1
def diffset(A):
    return set(a-b for a in A for b in A)
def add_energy(A):  # #{(a,b,c,d): a+b=c+d}
    from collections import Counter
    s=Counter(a+b for a in A for b in A)
    return sum(c*c for c in s.values())
def is_dilated_ap(A):
    As=sorted(A); d=As[1]-As[0]
    return all(As[i+1]-As[i]==d for i in range(len(As)-1))

# ---- collect M<1/13 covering families ----
random.seed(3)
fams=[]
for m in range(1,9): fams.append(sorted(list(range(1,13))+[182*m]))
for w in range(150,900):
    V=sorted(list(range(1,13))+[w])
    if len(set(V))==13: fams.append(V)
for _ in range(5000):
    V=sorted(random.sample(range(1,240),13))
    fams.append(V)
for w1 in [169,182,195,338,364,507]:
    for w2 in range(160,420,3):
        fams.append(sorted(list(range(1,12))+[w1,w2]))
# dedup
seen=set(); F0=[]
for V in fams:
    t=tuple(sorted(set(V)))
    if len(t)==13 and t not in seen: seen.add(t); F0.append(sorted(t))

below=[]
for V in F0:
    if not prim(V) or not cov(V): continue
    M=exact_M(V)
    if M<F(1,13):
        core=sorted(V)[:-1]  # V minus v_max
        below.append((M,V,core))

print("="*90)
print(f"M<1/13 covering families collected: {len(below)}")
print("="*90)
print("\n(F) FREIMAN: |core-core| for M<1/13 families (AP minimum for 12 elts = 23):")
sizes=set()
notap=[]
for M,V,core in below:
    dc=len(diffset(core))//1  # symmetric; |A-A| counts 0 and negatives => 2*11+1=23 for AP
    sizes.add(dc)
    if not is_dilated_ap(core): notap.append((M,V,core,dc))
print(f"   distinct |core-core| values over M<1/13 families: {sorted(sizes)}")
print(f"   (AP {{1..12}} has |A-A| = {len(diffset(list(range(1,13))))})")
print(f"   cores that are NOT dilated APs: {len(notap)}"
      + ("" if not notap else f"  <-- would REFUTE INV"))
for M,V,core,dc in notap[:5]:
    print(f"      M={M} |core-core|={dc} core={core} V={V}")

print("\n(S) THE M-SPECTRUM in (1/14, 1/13):  distinct M values (sorted):")
spec=sorted(set(M for M,_,_ in below))
for M in spec[:20]:
    # is it 14m/(182m+1)?
    num,den=M.numerator,M.denominator
    form = (den==13*num+1)
    print(f"   M={str(M):10s}={float(M):.6f}  num={num} den={den}  is 14m/(182m+1)-form (den=13num+1): {form}")
print(f"   total distinct M<1/13 covering values: {len(spec)}")
allform=all(M.denominator==13*M.numerator+1 for M in spec)
print(f"   ALL of the form m/(13m+1) [i.e. den=13*num+1]? {allform}  (=> discrete spectrum, CF [0;13,num])")

print("\n(E) ADDITIVE ENERGY: M<1/13 core energy vs a generic 12-set:")
if below:
    e_ap=add_energy(list(range(1,13)))
    print(f"   AP {{1..12}} additive energy = {e_ap}")
    gens=[add_energy(random.sample(range(1,60),12)) for _ in range(5)]
    print(f"   generic 12-set energies (5 samples): {gens}")
    core_es=set(add_energy(core) for _,_,core in below)
    print(f"   M<1/13 core energies (distinct): {sorted(core_es)}  (AP is the MAX for 12 elts)")

print("\nREADING: if every M<1/13 core is a dilated AP with |core-core|=23 (minimal) and")
print("the M<1/13 spectrum is exactly the discrete {14m/(182m+1)}, INV is 'M<1/13 => core")
print("has minimal difference set' = Freiman 3k-4 with the SHARP constant. The remaining")
print("content: M<1/13 (a Diophantine/resonance condition) => minimal |core-core| (additive).")
