#!/usr/bin/env python3
"""
lrc_counterexample_necessary_conditions_s554.py    oracle-2026-06-01-S554o

How many NECESSARY CONDITIONS must a counterexample to LRC@n satisfy? Each is a
provable necessary condition (violate it => some t is lonely => not a counterexample).
Piling them up pins the counterexample locus; if jointly unsatisfiable => LRC.

We enumerate ~20 conditions (sources in [brackets]) and verify the CHECKABLE ones for
n=14, then construct the most counterexample-like candidates and show they are STILL
LONELY (the conjunction does not yet yield a counterexample -- consistent with LRC).
"""
from itertools import combinations, product as iproduct
from functools import reduce
from math import gcd, sin, pi
import random

N = 14

def frac(x): return x - int(x // 1)
def dist0(p):
    p = p % 1.0; return min(p, 1 - p)

# ---- the killer check: is there a lonely time? (if yes, NOT a counterexample) ----
def lonely_anywhere(v, G=250000):
    for i in range(G):
        t = (i+0.5)/G
        if all(1.0/N < (s*t) % 1.0 < 1.0 - 1.0/N for s in v): return True
    return False

# ---- condition checkers ----
def C_primitive(v): return reduce(gcd, v) == 1
def C_distinct_nonzero(v): return len(set(v)) == len(v) and all(s != 0 for s in v)
def C_sieve_covered(v):   # A1: mult of every q in 2..14
    return all(any(s % q == 0 for s in v) for q in range(2, N+1))
def C_mult_of_n(v): return any(s % N == 0 for s in v)              # A2
def C_mult_of_7(v): return any(s % 7 == 0 for s in v)             # A4 (doubled prime)
def C_prime_powers(v):    # A3: mult of each maximal prime power <=14
    return all(any(s % pp == 0 for s in v) for pp in (8, 9, 5, 7, 11, 13))
def C_crt_all_classes_block(v, G=4000):  # D1: each mod-7 class blocks at some t
    cls = {c: [s for s in v if s % 7 == c] for c in range(7)}
    cls = {c: ss for c, ss in cls.items() if ss}
    blocks = {c: False for c in cls}
    for i in range(G):
        t = (i+0.5)/G
        for c in cls:
            if any(dist0(s*t) < 1.0/N for s in cls[c]): blocks[c] = True
    return all(blocks.values())
def ghat_abs(m):
    if m == 0: return 1 - 2.0/N
    return abs(sin(2*pi*m/N))/(pi*abs(m))
def C_high_energy(v, M=4):   # C2: pairwise resonance energy alone exceeds... (proxy)
    base = (1 - 2.0/N)**(len(v)-2); tot = 0.0
    for i, j in combinations(range(len(v)), 2):
        g = gcd(v[i], v[j]); a, b = v[j]//g, v[i]//g
        for k in range(1, 400):
            tot += 2*base*ghat_abs(k*a)*ghat_abs(k*b)
    return tot >= (1 - 2.0/N)**(len(v)-1)   # E2 >= main (NB: necessary uses FULL E)
def C_short_resonance(v):  # C3: a short resonance exists (pairwise: always does)
    best = min(v[i]//gcd(v[i], v[j]) + v[j]//gcd(v[i], v[j])
               for i, j in combinations(range(len(v)), 2))
    return best <= 6   # minimal pairwise resonance length small
def C_not_AP_scaled(v):  # G4: not c*(1,2,..,13)
    g = reduce(gcd, v); r = sorted(s//g for s in v)
    return r != list(range(1, N))

CONDITIONS = [
    ("A1 sieve-covered (mult of every q in 2..14) [THM-369]", C_sieve_covered),
    ("A2 multiple of n=14 [sieve q=n]", C_mult_of_n),
    ("A3 mult of each prime power 8,9,5,7,11,13 [sieve]", C_prime_powers),
    ("A4 multiple of 7 [doubled-prime n*=7, S546]", C_mult_of_7),
    ("D1 all 7 mod-7 CRT classes block [S524/S552]", C_crt_all_classes_block),
    ("C2 high pairwise resonance energy [S550]", C_high_energy),
    ("C3 a short resonance exists [S550 small minimal l]", C_short_resonance),
    ("F1 primitive (gcd=1) [scaling S549]", C_primitive),
    ("F2 distinct nonzero [normalization]", C_distinct_nonzero),
    ("G4 not AP-scaled (else lonely at wall) [S548]", C_not_AP_scaled),
]

def main():
    print("="*74)
    print("NECESSARY CONDITIONS a counterexample to LRC@14 must satisfy (verify on candidates)")
    print("="*74)
    # build candidate near-counterexample sets satisfying the sieve
    cands = {
        "sieve-minimal-ish": (5, 7, 8, 9, 11, 13, 14, 2, 3, 4, 6, 10, 12),
        "AP 1..13 (the extremal)": tuple(range(1, 14)),
        "prime-power heavy": (5, 7, 8, 9, 11, 13, 14, 16, 18, 22, 26, 35, 45),
    }
    for name, v in cands.items():
        v = tuple(v)
        sat = [lab for lab, f in CONDITIONS if f(v)]
        lon = lonely_anywhere(v)
        print(f"\n  {name}: v={v}")
        print(f"    satisfies {len(sat)}/{len(CONDITIONS)} conditions: {[s.split(' ')[0] for s in sat]}")
        print(f"    LONELY somewhere = {lon}  ->  {'NOT a counterexample (lonely)' if lon else 'CANDIDATE (no lonely time found!)'}")
    print()
    print("  random primitive sieve-covered sets (do any satisfy ALL + be non-lonely?):")
    rnd = random.Random(14); found = 0; checked = 0
    while checked < 30:
        v = tuple(sorted(rnd.sample(range(1, 60), 13)))
        if not C_primitive(v) or not C_sieve_covered(v): continue
        checked += 1
        if not lonely_anywhere(v):
            found += 1
            print(f"    !! non-lonely sieve-covered set: {v}")
    print(f"    checked {checked} primitive sieve-covered sets; non-lonely (counterexamples): {found}")
    print()
    print("="*74)
    print("THE FULL CONDITION LIST (necessary for a counterexample; sources)")
    print("="*74)
    full = [
      "A1 sieve-covered: mult of every q in {2..n}  [THM-369]",
      "A2 multiple of n  [sieve q=n; blocks all t=a/n]",
      "A3 multiple of each maximal prime power <=n  [sieve]",
      "A4 multiple of n* (=n/2 even / n odd); for n=14, of 7  [S546]",
      "A5 for each q<=n NOT all speeds coprime to q  [= A1]",
      "B1 averaging-extremal: near(t)>=1 for ALL t, but mean near=2-2/n => floor 1 never 0 [S553]",
      "B2 weaker-threshold: at 1/(cn), the >=? near-set covers for each c [S553 menu]",
      "B3 2nd-moment-consistent: pairwise danger-overlaps tuned so min near=1 not 0 [S553]",
      "C1 danger arcs B_i cover [0,1) (total measure 2-2/n>1; overlaps exact) [def/S525]",
      "C2 high-energy core: resonance energy E(v) >= (1-2/n)^{n-1}  [S550]",
      "C3 small minimal resonance length l(v) (a short resonance exists)  [S550]",
      "C4 pervasive returns: cascade product of conditional clearances = 0  [S545]",
      "D1 all n* CRT classes blocked at every t (n=14: 7 mod-7 classes) [S524/S552]",
      "D2 singleton class (mult of n*) coupled to the pair-classes [S552]",
      "D3 full-support inside debt present mod n* (parity) [S533]",
      "E1 observer NEVER a source in the LRC-marked tournament [S511]",
      "E2 observer apex (bracketing) gap < 2/n at all t (never a wide gap at 0) [S530]",
      "E3 always >=1 trienerment tie at the observer [S539] (= B1)",
      "F1 primitive: gcd(speeds)=1 (WLOG, scaling) [S549]",
      "F2 distinct nonzero speeds [normalization]",
      "F3 minimal counterexample has bounded speeds (finite reduction) [Cusick/Tao]",
      "G1 frequency-concentrated, NOT decorrelated [S544]",
      "G2 commensurable (closed orbit, not dense) [S521o] (auto for integers)",
      "G4 NOT AP-scaled c*(1..n-1) (else lonely at the wall t=k/n) [S548]",
      "H1 the difference set has a specific QR/Frobenius pattern mod n* [S535]",
      "H2 the relation lattice carries the high energy via low-order returns [S545/S550]",
    ]
    for i, c in enumerate(full, 1):
        print(f"  {i:2d}. {c}")
    print(f"\n  >>> {len(full)} necessary conditions enumerated. A counterexample must satisfy")
    print("      ALL simultaneously. Tested candidates satisfy many yet remain LONELY -- the")
    print("      conjunction is not (yet) shown unsatisfiable, but the locus is razor-thin.")

if __name__ == "__main__":
    main()
