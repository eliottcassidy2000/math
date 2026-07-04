#!/usr/bin/env python3
"""
THE WITNESS-OVERLAP KERNEL as a RESONANCE-COUNTING argument (mac-mini-2026-07-03-S28).
Creative reason (circle method): the number of lonely a at denominator q is
  N(q) = sum_{a=1}^{q-1} prod_i 1[ ||v_i a / q|| >= 1/14 ]
       = q * (6/7)^13  +  q * E(q),   E(q) = sum_{m != 0, q | sum m_i v_i} prod_i shat(m_i),
where shat(m) = -sin(pi m / 7)/(pi m) (m!=0), shat(0)=6/7 (safe-arc Fourier).
CLAIM: N(q)=0 (no witness) REQUIRES a SMALL-coefficient resonance q | (sum m_i v_i), |m_i| <= K, because the
large-|m| tail of E(q) is bounded below -(6/7)^13 in magnitude (singular-series convergence). Small combos
sum m_i v_i (|m|<=K) are integers <= K*max|v|; their divisors number O(log max|v|). So NO-WITNESS moduli
sit on O(log M) resonant moduli => the first witness q* = O(log M).

VERIFY: for compressed lcm-blockers, is every no-witness modulus a small-combo divisor? count resonant moduli.
"""
from math import gcd, log
from functools import reduce
from itertools import product
import random

def gcd_all(xs): return reduce(gcd, xs)
def lcm(a,b): return a*b//gcd(a,b)
def nd(x): x=x%1; return min(x,1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))

def has_witness(sp, q):
    for a in range(1, q):
        if gcd(a,q)!=1: continue
        if all(nd(v*a/q)>=1/14 for v in sp): return True
    return False

def small_combo_divisors(sp, K, qlo, qhi):
    """moduli q in [qlo,qhi] dividing some nonzero sum m_i v_i with |m_i|<=K and at most 2 nonzero m_i."""
    res = set()
    n = len(sp)
    # 1-term: q | v_i
    for v in sp:
        for q in range(qlo, qhi+1):
            if v % q == 0: res.add(q)
    # 2-term: q | (a*v_i + b*v_j), |a|,|b|<=K
    for i in range(n):
        for j in range(i+1, n):
            for a in range(-K, K+1):
                for b in range(-K, K+1):
                    if a==0 and b==0: continue
                    val = a*sp[i] + b*sp[j]
                    if val==0: continue
                    av = abs(val)
                    for q in range(qlo, qhi+1):
                        if av % q == 0: res.add(q)
    return res

def build_lcm_blocker(N, rng):
    mods=list(range(15,80)); rng.shuffle(mods); runners=[]; i=0
    while i<len(mods) and len(runners)<13:
        grp=[]; L=1
        while i<len(mods):
            L2=lcm(L,mods[i])
            if L2<=N and len(grp)<8: L=L2; grp.append(mods[i]); i+=1
            else: break
        if L>1: runners.append(L*max(1,round(N/L)))
        else: i+=1
    for q in [8,9,5,7,11,13,2,3,4,6]:
        if len(runners)>=13: break
        if not any(s%q==0 for s in runners): runners.append(q)
    while len(runners)<13: runners.append(rng.randint(2,22))
    return sorted(set(runners))[:13]

if __name__ == "__main__":
    rng = random.Random(28)
    print("Resonance-kernel verification on compressed lcm-blockers.")
    print("=" * 88)
    print(f"{'N':>12} {'witness q*':>11} {'#no-wit q in[15,q*]':>20} {'#(1&2-term reso)':>17} {'no-wit SUBSET reso?':>19}")
    for N in [10**4, 10**6, 10**9, 10**12]:
        # pick the worst (highest q*) family at this N
        worst=None; wq=0
        for _ in range(2000):
            sp=build_lcm_blocker(N, rng)
            if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
            q=2
            while q<=400 and not has_witness(sp,q): q+=1
            if q<=400 and q>wq: wq=q; worst=sp
        if worst is None: continue
        sp=worst; qstar=wq
        nowit = [q for q in range(15, qstar) if not has_witness(sp, q)]
        reso = small_combo_divisors(sp, K=7, qlo=15, qhi=qstar-1)
        subset = all(q in reso for q in nowit)
        print(f"{N:>12} {qstar:>11} {len(nowit):>20} {len(reso):>17} {str(subset):>19}")
        if not subset:
            missing=[q for q in nowit if q not in reso]
            print(f"       no-wit moduli NOT explained by 1&2-term small resonance: {missing[:10]}")
    print("\n=> if no-wit SUBSET reso (True) and #reso = O(log N): 'no witness at q' => small resonance,")
    print("   and small resonances number O(log M) (divisors of ~K^13 combos <= K*M). So witness q*=O(log M):")
    print("   the CRT/divisor capacity of 13 speeds<=M bounds the blockable moduli. This is the creative reason.")
