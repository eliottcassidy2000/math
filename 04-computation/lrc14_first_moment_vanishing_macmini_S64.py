"""
mac-mini-2026-07-09-S64 -- DECISIVE TEST of the first-moment-vanishing path for the
7-structured dissociated good period.

Reframe (HYP-5600 + LEM-011): W(x) = sum_gaps (g - 1/7)_+ = uncovered measure; W(x)>0 <=> x good.
THM-664 Weyl:  E_grid[W] over j=0..V-1  =  sum_{n : V | n.e} What(n),  with What(0)=(6/7)^k.
LEM-011: What(n) VANISHES when 7|n_i (any i) or 7|N (N=sum n_i) -- "7-commensurate".

CLAIM to test: for 7-structured E (many e_i = 0 mod 7), the surviving corrections are small/zero,
so E_grid[W] stays ~ (6/7)^k > 0  ==> SOME j is good (first moment positive => a positive summand).
If E_grid[W] over j=1..V-1 is bounded below by a positive constant, the good period is FORCED with
NO moment certificate (which MISTAKE-128 killed). Split by gcd(7,V).

We compute E_grid[W] directly (exact rationals via integer arithmetic on the V-grid) and compare
to (6/7)^k. Decisive: is avg_{j>=1} W(j/V) >= c > 0 uniformly?
"""
from math import gcd
from functools import reduce
from fractions import Fraction
import random
random.seed(6464)

def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best

def W_at(E, j, V):
    """W(j/V) = sum over circular gaps of (gap - 1/7)_+, exact Fraction. phases = {e*j mod V}/V."""
    ph = sorted({(e*j) % V for e in E})
    m = len(ph)
    if m == 0: return Fraction(0)
    if m == 1:
        # single phase: one gap of length 1 (whole circle minus the point)
        g = Fraction(1)
        return max(g - Fraction(1,7), Fraction(0))
    tot = Fraction(0)
    for i in range(m):
        gap = Fraction((ph[(i+1) % m] - ph[i]) % V, V)
        d = gap - Fraction(1,7)
        if d > 0: tot += d
    return tot

def avg_W(E, V, include0):
    lo = 0 if include0 else 1
    s = Fraction(0); cnt = 0
    for j in range(lo, V):
        s += W_at(E, j, V); cnt += 1
    return s / cnt if cnt else Fraction(0)

def make_7struct(k, s, min_s7):
    sevens=[x for x in range(0,s+1,7)]
    if len(sevens)<min_s7: return None
    S7=sorted(random.sample(sevens, random.randint(min_s7, min(len(sevens),k-1))))
    rest_pool=[x for x in range(1,s) if x%7!=0]
    need=k-len(S7)
    if need<0 or len(rest_pool)<need: return None
    E=sorted(set(S7+[0,s]+random.sample(rest_pool,max(0,need))))
    return E if len(E)==k else None

k = 13
target = Fraction(6,7)**k
print(f"k={k}, (6/7)^k = {float(target):.6f}\n")

for label, sevenDivV in (("7 | Vmax", True), ("gcd(7,Vmax)=1", False)):
    print(f"=== {label} ===")
    n=0; minavg1=None; minavg0=None; worst=None; below_target0=0
    for _ in range(3000):
        s=random.randint(30,120)
        E=make_7struct(k,s,4)
        if E is None or not prim(E) or longest_ap(E)>k-6: continue
        mx=max(E)
        if sevenDivV:
            V = 7*((mx)//7 + random.randint(1,6))
        else:
            V = mx + random.randint(1, max(1,mx//4))
            while V % 7 == 0: V += 1
        if V <= mx or V < 14: continue
        n += 1
        a1 = avg_W(E, V, include0=False)   # j=1..V-1  (real periods)
        a0 = avg_W(E, V, include0=True)    # j=0..V-1  (full Weyl average)
        if minavg1 is None or a1 < minavg1: minavg1 = a1; worst = (list(E), V, float(a1), float(a0))
        if minavg0 is None or a0 < minavg0: minavg0 = a0
        if a0 < target: below_target0 += 1
        if n >= 250: break
    if n:
        print(f"  sampled {n} sets")
        print(f"  min over sets of  E_grid[W] (j=1..V-1) = {float(minavg1):.6f}   (target (6/7)^k={float(target):.6f})")
        print(f"  min over sets of  E_grid[W] (j=0..V-1) = {float(minavg0):.6f}")
        print(f"  sets with full-grid avg < (6/7)^k: {below_target0}/{n}  (0 => corrections are >=0 i.e. help)")
        print(f"  worst (min j>=1 avg): E={worst[0][:8]}... V={worst[1]}  avg1={worst[2]:.5f} avg0={worst[3]:.5f}")
    print()

print("INTERPRETATION:")
print(" - If min E_grid[W](j=1..V-1) > 0 uniformly  => first moment forces a good period (DONE, no moments).")
print(" - If full-grid avg >= (6/7)^k (below_target0=0) => the surviving corrections are NONNEGATIVE")
print("   for 7-structured sets: the 7-commensurate vanishing (LEM-011) leaves only helpful terms.")
