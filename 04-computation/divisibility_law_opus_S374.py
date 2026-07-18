# opus-2026-07-17-S374 -- THE ARITHMETIC POSITION LAW.
#
# S373 asked where the uncovered gaps SIT.  Answer, from the refutation above:
# their position is governed by the DIVISIBILITY structure of V.  Define
#     q0(V) = the smallest modulus dividing NO speed
# (so D(V) = q0 - 1 is the largest Q with every q <= Q dividing some speed).
# Blocking q requires a speed divisible by q, so no lonely rational can have
# denominator < ... well, that is not quite right either -- q dividing some
# speed does not by itself block q, since a DIFFERENT p may still work.  So
# there are two questions:
#   (A) how often is min-denominator EXACTLY q0?  (the extended sieve lemma:
#       "q divides no speed => some p/q is lonely", classically proved only
#       for q <= 14)
#   (B) when it fails, by how much?
from math import gcd
from functools import reduce
import random
from collections import Counter
def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
def works(V,q):
    return any(lonely_at(V,p,q) for p in range(1,q) if gcd(p,q)==1)
def min_den(V,Qmax=400):
    for q in range(1,Qmax+1):
        if works(V,q): return q
    return None
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q

print("(4) THE LAW: is min-denominator = q0 (the first modulus dividing no speed)?")
random.seed(3741)
C=Counter(); big=[]
n=0
for _ in range(500):
    V=sorted(random.sample(range(2,300),13))
    if reduce(gcd,V)!=1: continue
    z=q0(V); m=min_den(V)
    if m is None: continue
    n+=1; C[m-z]+=1
    if m!=z: big.append((z,m,V))
print(f"    {n} primitive families")
for d in sorted(C): print(f"      min-den - q0 = {d:3d}:  {C[d]:4d}  ({100*C[d]/n:5.1f}%)")
print(f"    EXACT (min-den = q0): {100*C[0]/n:.1f}%")

print()
print("(5) WHEN THE LAW FAILS -- is it always q0 > 14?")
print("    (the classical sieve lemma proves q divides no speed => 1/q lonely")
print("     ONLY for q <= 14, since ||v/q|| >= 1/q needs 1/q >= 1/14)")
le14=sum(1 for z,m,V in big if z<=14); gt14=sum(1 for z,m,V in big if z>14)
tot14=sum(1 for _ in range(0))
allz=[]
random.seed(3741)
cnt_le=cnt_gt=fail_le=fail_gt=0
for _ in range(500):
    V=sorted(random.sample(range(2,300),13))
    if reduce(gcd,V)!=1: continue
    z=q0(V); m=min_den(V)
    if m is None: continue
    if z<=14:
        cnt_le+=1
        if m!=z: fail_le+=1
    else:
        cnt_gt+=1
        if m!=z: fail_gt+=1
print(f"    q0 <= 14: {cnt_le:4d} families, law fails {fail_le:3d} ({100*fail_le/max(cnt_le,1):.1f}%)")
print(f"    q0 >  14: {cnt_gt:4d} families, law fails {fail_gt:3d} ({100*fail_gt/max(cnt_gt,1):.1f}%)")
print()
print("    => if the law NEVER fails for q0 <= 14, that is the classical lemma")
print("       confirmed; the open content is exactly the q0 > 14 regime.")
if big:
    print()
    print("(6) EXAMPLES WHERE IT FAILS")
    for z,m,V in big[:5]:
        print(f"      q0 = {z:3d}, min-den = {m:3d}   V = {V}")
