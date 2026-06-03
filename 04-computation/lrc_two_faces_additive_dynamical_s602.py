#!/usr/bin/env python3
"""
claudebox-2026-06-03-S602 : everything from block-diagonalization — the TWO FACES of the worry-set.

The user's insight: the worry-set's STATIC MODULAR RIGIDITY (the 64 self-converse classes,
realized by antipodal transversals mod 2n-1 = mod 27 = 3^3 at n=14) is the ADDITIVE-FACE shadow and
SURVIVES; the DYNAMICAL DOUBLING face (x->2x mod n) FRAGMENTS at 2q — the same 2·7 seam from both
faces.

Block-diagonalization gives two decompositions of the same worry-set:
  ADDITIVE face   = residue structure mod 2n-1 (the antipodal transversal / 64-flip-lattice).
                    2n-1 is ALWAYS ODD => doubling is a permutation, NO 2-block, often max-mixing.
                    Static modular rigidity lives here and SURVIVES the seam.
  DYNAMICAL face  = the doubling map x->2x mod n (the phase dynamics, HYP-2117/S585/S596).
                    Even n => 2-adic COLLAPSE => the apex = the rank-1 2-block (S599) = where the
                    single-corrector obstruction lives. FRAGMENTS at the seam.

KEY (verified): at n=14 the dynamical face mod 14 collapses, but the additive face mod 27=3^3 is
MAXIMALLY MIXING (2 is a primitive root mod 3^3). The apex obstruction = the 2-block of the mod-n
gcd-divisor decomposition (S599); mod 2n-1 (odd) there is NO 2-block => the additive face is
APEX-FREE. The pair-sum multi-sieve (HYP-2075) works because it lives in this odd, apex-free face.
"""
from math import gcd

def units(m): return [u for u in range(1, m) if gcd(u, m) == 1]
def phi(m): return len(units(m))
def ord_mult(g, m):
    if gcd(g, m) != 1: return None
    x = g % m; k = 1
    while x != 1: x = (x * g) % m; k += 1
    return k
def doubling_cosets(m):
    if m % 2 == 0: return None
    seen = set(); sizes = []
    for u in units(m):
        if u in seen: continue
        orb = set(); x = u
        while x not in orb: orb.add(x); x = (2 * x) % m
        seen |= orb; sizes.append(len(orb))
    return sorted(sizes, reverse=True)
def divisors(m): return [d for d in range(1, m + 1) if m % d == 0]

def main():
    print("=" * 96)
    print("S602  the two faces of the worry-set: additive (mod 2n-1, survives) vs dynamical (mod n)")
    print("=" * 96)

    print("\n[1] THE DOUBLING DUALITY: dynamical face mod n vs additive face mod 2n-1")
    print("  n  | mod n (dynamical)              | mod 2n-1 (additive)            | seam")
    for n in range(6, 23):
        m2 = 2 * n - 1
        dyn = "COLLAPSE (even, 2-adic)" if n % 2 == 0 else f"perm, cosets {doubling_cosets(n)}"
        o2 = ord_mult(2, m2); p = phi(m2)
        add = "MAX-MIXING (2 prim root)" if o2 == p else f"ord(2)={o2}/{p}"
        seam = "<- 2-adic seam (apex)" if n % 2 == 0 else ""
        print(f"  {n:2d}  | {dyn:30s} | mod {m2:2d}: {add:24s} | {seam}")
    print("  => n=14: dynamical mod 14 DEAD, additive mod 27=3^3 MAXIMALLY ALIVE. The static")
    print("     rigidity survives in the additive face precisely where the dynamical one collapses.")

    print("\n[2] THE ADDITIVE FACE IS APEX-FREE: the apex = the 2-block of mod-n; mod 2n-1 has none")
    print("  n  | mod-n divisor blocks (2-block?) | mod-(2n-1) divisor blocks (no even ⇒ no 2-block)")
    for n in [10, 14, 18, 22]:
        m2 = 2 * n - 1
        dn = [d for d in divisors(n) if d > 1]
        d2 = [d for d in divisors(m2) if d > 1]
        has2_n = any(d % 2 == 0 for d in dn)
        has2_m = any(d % 2 == 0 for d in d2)
        print(f"  {n:2d}  | {str(dn):26s} 2-block={has2_n}  | mod {m2}: {str(d2):14s} 2-block={has2_m}")
    print("  => mod n (even) carries the rank-1 2-block = the apex obstruction (S599/HYP-2063);")
    print("     mod 2n-1 (odd) has NO even divisor ⇒ NO 2-block ⇒ APEX-FREE. The multi-sieve lives here.")

    print("\n[3] STATIC MODULAR RIGIDITY: the 64 = 2^((n-2)/2) self-converse flip-lattice (additive)")
    print("  the antipodal pairs {k, n-k} mod n (k=1..(n-2)/2) + apex n/2; orienting them = the worry classes")
    for n in [6, 8, 10, 14]:
        half = (n - 2) // 2
        pairs = [(k, n - k) for k in range(1, half + 1)]
        cls = 2 ** half
        print(f"  n={n:2d}: {half} antipodal pairs {pairs} + apex {n//2}  ⇒  2^{half} = {cls} self-converse classes "
              f"(AP at the tight bottom, all flips loosen)")
    print("  => the flip-lattice is the additive-face combinatorial skeleton; rigid and survives.")

    print("\n[4] THE SEAM FROM BOTH FACES at n=14")
    print("  DYNAMICAL (mod 14): doubling 2-adic collapse; apex 7 = rank-1 2-block; single-corrector dies.")
    print("  ADDITIVE  (mod 27=3^3): 2 is a primitive root ⇒ ONE max-mixing 18-cycle; no 2-block;")
    print("                          8191 transversals all lonely; 64-flip-lattice rigid. SURVIVES.")
    o = ord_mult(2, 27)
    print(f"  verify: ord_27(2) = {o} = φ(27) = {phi(27)}  ⇒  2 IS a primitive root mod 27 (max mixing).")
    print("  SAME 2·7 seam: the 2 kills the dynamical face; the additive face (odd 2n-1) dissolves the apex.")

if __name__ == "__main__":
    main()
