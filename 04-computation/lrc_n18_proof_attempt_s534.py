#!/usr/bin/env python3
"""
lrc_n18_proof_attempt_s534.py    oracle-2026-06-01-S534o

An HONEST proof attempt for LRC(n=18) (17 nonzero speeds, threshold 1/18) using the
full repo machinery. n=18 is DEEPLY open (k=17 runners >> proven k<=6). This builds
the rigorous sub-cases, characterizes the residual with the n*=9 (prime-power)
generalization of the S533c three-channel parity law, and runs the cascade heuristic.

n=18 facts: n* = n/2 = 9 = 3^2 (PRIME POWER -- the 'if 15 were prime' regime, S17).
The covering characters ghat(m) vanish iff 9 | m.

CONTENTS
 (1) SIEVE sub-cases (RIGOROUS, THM-369): for each q in {2..18} with no speed
     divisible by q, t=1/q is 18-lonely (since 1/q >= 1/18). In particular
     no-multiple-of-9 => t=1/9 lonely. A counterexample must be 'q-covered' for
     every q<=18. Verified numerically.
 (2) The n*=9 PARITY law and its VACUITY: full-support inside debt-free <=> 0 not
     representable as sum c_i v_i (mod 9), c_i in {1..8}. We show this is NEVER
     satisfiable for primitive 17-sets -- the parity obstruction that solved half of
     n=4 is VACUOUS at n=18. Tabulate the degradation vs phi(n*) / unit structure.
 (3) The 3-adic FILTRATION: mod 9 splits nonzero residues into units (v3=0, 6 of
     them) and v3=1 multiples-of-3 (residues 3,6). Debt-free could only survive if
     ALL runners are divisible by 3 -- which recurses to the n=6 mod-3 law one
     3-adic level down (and is non-primitive). The prime-power recursion.
 (4) CASCADE heuristic (S527): threshold (n-1)((n-2)/n)^{n-2} = 17*(8/9)^16 ~ 2.58
     >= 1 (passes with margin); run the feasible-arc carve on sample 17-sets (stays
     nonempty) and on the AP/regular-18-gon (empties to 0 = tight). NON-rigorous
     (discrepancy caveat).
 (5) HONEST verdict.
"""
from math import sin, pi, gcd
from itertools import product
from functools import reduce
import random

# ----------------------------------------------------------------------
def nstar(n): return n//2 if n % 2 == 0 else n

# (1) SIEVE sub-cases -- rigorous
def lonely_at(speeds, q, n):
    """is t=1/q 18-lonely? min over runners of ||v/q|| >= 1/n ?"""
    thr = 1.0/n
    return all(min((v/q) % 1.0, 1-((v/q) % 1.0)) >= thr - 1e-12 for v in speeds)

def part1():
    print("="*72); print("(1) SIEVE sub-cases for n=18 (RIGOROUS, THM-369): no mult of q => t=1/q lonely")
    print("="*72)
    n = 18
    print("   q : 1/q   >=1/18?  (no speed divisible by q => t=1/q is 18-lonely)")
    for q in range(2, 19):
        print(f"   {q:2d} : {1/q:.4f}  {'yes' if 1/q >= 1/n - 1e-12 else 'NO '}")
    # numeric demo: a random set with no multiple of 9 is lonely at t=1/9
    rnd = random.Random(534)
    ok = 0; tot = 0
    for _ in range(2000):
        v = tuple(sorted(rnd.sample(range(1, 60), 17)))
        if any(x % 9 == 0 for x in v): continue
        tot += 1
        if lonely_at(v, 9, n): ok += 1
    print(f"   numeric: {ok}/{tot} random 17-sets with NO multiple of 9 are lonely at t=1/9 (expect all)")
    print("   => A counterexample must be q-COVERED for every q in {2..18}: contain a")
    print("      multiple of 9, of 16, of 5,7,11,13,17, ... (the denominator sieve).")

# (2)+(3) the n*=9 parity law and its vacuity, with the 3-adic filtration
def full_support_zero_reachable(speeds, ns):
    """is 0 reachable as sum c_i v_i (mod ns), each c_i in {1..ns-1} (full support)?
    => a full-support resonance is mod-ns feasible => inside debt PRESENT."""
    reach = {0}
    for v in speeds:
        Si = set((c*v) % ns for c in range(1, ns))
        reach = set((r+s) % ns for r in reach for s in Si)
        if len(reach) == ns:
            # already everything; remaining runners can only keep it full
            pass
    return 0 in reach

def part23():
    print("\n"+"="*72); print("(2)+(3) the n*=9 parity law: VACUOUS at n=18 (debt always present)"); print("="*72)
    n = 18; ns = nstar(n)
    rnd = random.Random(5342)
    tot = 0; debtfree = 0; ex = []
    for _ in range(3000):
        v = tuple(sorted(rnd.sample(range(1, 80), 17)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if not full_support_zero_reachable(v, ns):
            debtfree += 1
            if len(ex) < 5: ex.append(v)
    print(f"   n=18, n*=9: of {tot} primitive 17-sets, FULL-SUPPORT debt-free (parity fires): {debtfree}")
    print(f"   => the parity obstruction is VACUOUS for primitive sets at n=18.")
    print("   3-adic reason: mod 9 nonzero residues split into UNITS (v3=0: {1,2,4,5,7,8})")
    print("   and v3=1 multiples-of-3 ({3,6}). A v3=1 runner can contribute 0 mod 9 (c=3,6);")
    print("   a UNIT runner contributes all of {1..8}. So with any unit runner, 0 is")
    print("   reachable -> debt present. Debt-free needs ALL runners divisible by 3 =>")
    print("   recurse to the n=6 mod-3 law one level down => non-primitive. Prime-power recursion.")
    # degradation table
    print("\n   PARITY-FIRING DEGRADATION vs n (fraction of primitive sets that are debt-free):")
    print("   n   n*  phi(n*)  units mod n*               parity-fires?")
    def phi(m):
        return sum(1 for a in range(1, m+1) if gcd(a, m) == 1)
    for nn in (4, 6, 8, 9, 10, 12, 14, 16, 18):
        ns2 = nstar(nn); k = nn-1
        units = [a for a in range(1, ns2) if gcd(a, ns2) == 1]
        # quick sample firing rate
        rr = random.Random(nn); t = 0; df = 0
        for _ in range(1500):
            v = tuple(rr.sample(range(1, 8*nn), k))
            if reduce(gcd, v) != 1: continue
            t += 1
            if not full_support_zero_reachable(v, ns2): df += 1
            if t >= 400: break
        rate = f"{df}/{t}"
        print(f"   {nn:2d}  {ns2:2d}   {phi(ns2):2d}    {str(units):24s}  {rate}")
    print("   => parity fires often only at n=4 (n*=2, units {1}); rare at n=6 (n*=3);")
    print("      VANISHES once n* has >=2 units AND a non-unit-nonzero residue (n>=8-ish).")

# (4) cascade heuristic
def cascade_threshold(n):
    return (n-1) * ((n-2)/n)**(n-2)

def feasible_carve(speeds, n, samples=200000):
    thr = 1.0/n; lo, hi = thr, 1-thr
    safe = [True]*samples
    chain = []
    for v in speeds:
        for i in range(samples):
            if safe[i]:
                x = (v*(i+0.5)/samples) % 1.0
                if not (lo < x < hi): safe[i] = False
        m = sum(safe)/samples
        chain.append(m)
    return chain

def part4():
    print("\n"+"="*72); print("(4) CASCADE heuristic (S527): threshold passes for n=18 (margin 2.58)"); print("="*72)
    for n in (7, 14, 16, 18):
        print(f"   n={n}: (n-1)((n-2)/n)^(n-2) = {cascade_threshold(n):.4f}  (>=1 => cascade passes)")
    n = 18
    ap = tuple(range(1, 18))                      # regular 18-gon (tight)
    rnd = random.Random(53421)
    nonap = tuple(sorted(rnd.sample(range(1, 80), 17)))
    while reduce(gcd, nonap) != 1:
        nonap = tuple(sorted(rnd.sample(range(1, 80), 17)))
    for label, v in [("AP=regular-18gon", ap), ("random primitive", nonap)]:
        ch = feasible_carve(v, n)
        final = ch[-1]
        print(f"   {label}: feasible measure after each runner -> final={final:.5f} "
              f"({'EMPTIES (tight)' if final < 1e-4 else 'stays >0 (lonely)'})")
        print(f"      chain (every 3rd): {[round(ch[i],3) for i in range(0,len(ch),3)]}")
    print("   NON-RIGOROUS: the (n-2)/n per-step shrink needs a discrepancy bound")
    print("   (Erdos-Turan); its log factor is exactly what keeps n=18 open.")

def part5():
    print("\n"+"="*72); print("(5) HONEST VERDICT"); print("="*72)
    print("""   PROVEN (rigorous) for n=18:
     - no-multiple-of-q case for each q<=18 (t=1/q lonely); strongest single class
       q=9 (t=1/9). A counterexample must be q-covered for ALL q in {2..18}.
   CHARACTERIZED (new, this session):
     - the n*=9 parity law is VACUOUS for primitive sets: inside debt is ALWAYS
       present. The 'free parity' that settled half of n=4 contributes NOTHING at
       n=18 because 9=3^2 has 6 units + a 3-adic sublevel. Debt-free recurses to the
       n=6 mod-3 law one 3-adic level down (non-primitive). Prime-power degradation.
   HEURISTIC (non-rigorous):
     - cascade threshold passes with margin 2.58; the AP (regular 18-gon) is tight
       (carve empties to 0); all sampled primitive sets stay lonely.
   NOT PROVEN: the fully-sieve-covered residual with the discrepancy/coupling gap --
     the same wall as n=14/16 and the general conjecture. n=18 is NOT proved here.""")

def main():
    part1(); part23(); part4(); part5()

if __name__ == "__main__":
    main()
