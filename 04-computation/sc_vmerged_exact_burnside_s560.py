"""
sc_vmerged_exact_burnside_s560.py — monad-researcher-2026-06-02-S560

EXACT, float-free Burnside computation of the THREE mirror sequences:
  A000568(n) = number of tournaments on n vertices (up to iso)
  SC(n)      = number of self-complementary (self-converse) tournaments on n
  V_merged(n)= (A000568(n) + SC(n)) / 2  = vertices of the merged metagraph G_n/Z_2

MOTIVATION:
  The exact values existed only in a STALE OUTPUT FILE
  (05-knowledge/results/sc_a000568_q_deformed.out, to n=30); the script that
  produced it was lost from 04-computation/. This script RECOVERS that
  computation from the canonical Burnside formulas, extends it to n=60, proves
  every value is an exact integer, and verifies two canon identities with fresh
  exact arithmetic:

    (I)  THM-283 mirror:  A000568 = Σ_{all parts odd} 2^c/z,
                          SC      = Σ_{parts ≡2 mod 4, +one 1 if n odd} 2^c/z.
    (II) MISTAKE-049 (corrected framing):  SC(2m) = A(m, 4),
         where A(m,q) = Σ_{odd λ ⊢ m} q^c(λ) / z(λ) is the q-deformed tournament
         count (A(m,2) = A000568(m)). This is just base-4 Burnside over odd
         partitions — the corrected identity that replaced the fabricated
         SC(n)=A000568(n-1) claim.

DAVIS / BURNSIDE EXPONENT (shared by all three):
  For a cycle type with cycles of lengths l_1,...,l_k (multiplicities m_i),
  the number of orbits of UNORDERED PAIRS {i,j} under the permutation is
    c = Σ_i  ⌊l_i / 2⌋                 # orbits of pairs INSIDE each cycle
      + Σ_{i<j} gcd(l_i, l_j)          # orbits of pairs BETWEEN two cycles
  NOTE the inside term is ⌊l/2⌋, which equals (l-1)/2 ONLY for odd l.
  A000568 uses odd parts only, so the sibling script
  a000568_exact_burnside.py legitimately writes (l-1)//2; for SC the parts are
  even (≡2 mod 4) and the ⌊l/2⌋ form is REQUIRED. This script uses l//2 = ⌊l/2⌋
  universally, which is correct for both parities.

  z(λ) = Π_i  l_i^{m_i} * m_i!   (the order of the centralizer of the cycle type)

Everything is accumulated with Fraction and asserted to have denominator 1.
"""
from fractions import Fraction
from math import gcd, factorial

# ---- known cross-check values (OEIS / repo canon) ----------------------------
KNOWN_A000568 = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536,
    10: 9733056, 11: 903753248, 12: 154108311168, 13: 48542114686912,
    14: 28401423719122304, 15: 31021002160355166848,
    16: 63530415842308265100288, 17: 244912778438520759443245824,
    18: 1783398846284777975419600287232,
    19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}
# THM-283 / anti_aut_edge_merge_s189.out, indexed by n
KNOWN_SC = {
    2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176, 9: 2752, 10: 8784,
    11: 279968, 12: 1492288, 13: 95458560, 14: 872687552, 15: 111698291584,
    16: 1787154671104, 17: 457509297625088, 18: 13013584213369088,
    19: 6662951988432581120,
}


# ---- exponent / weight machinery --------------------------------------------
def pair_orbit_exponent(lengths):
    """#orbits of unordered pairs under a permutation with given cycle lengths.
    Uses floor(L/2) inside each cycle (correct for odd AND even L)."""
    t = 0
    for L in lengths:
        t += L // 2                      # = floor(L/2) = (L-1)/2 for odd L
    for i in range(len(lengths)):
        for j in range(i + 1, len(lengths)):
            t += gcd(lengths[i], lengths[j])
    return t


def z_factor(lengths):
    """Returns Fraction(1, z(lambda)) where z = Π l^{m_l} * m_l!."""
    from collections import Counter
    denom = 1
    for L, m in Counter(lengths).items():
        denom *= (L ** m) * factorial(m)
    return Fraction(1, denom)


def burnside_sum(parts_lists, base):
    """Σ over the given cycle types of  base^c(λ) / z(λ)."""
    total = Fraction(0)
    for lengths in parts_lists:
        total += z_factor(lengths) * (base ** pair_orbit_exponent(lengths))
    return total


# ---- partition generators ----------------------------------------------------
def partitions_from(n, allowed, extra_fixed=0):
    """All multisets (as sorted tuples) of parts drawn from `allowed` summing to
    n - extra_fixed, with `extra_fixed` parts of size 1 appended.
    `allowed` is a sorted list of permissible part sizes."""
    target = n - extra_fixed
    results = []

    def rec(remaining, idx, current):
        if remaining == 0:
            results.append(tuple(current) + (1,) * extra_fixed)
            return
        for k in range(idx, len(allowed)):
            p = allowed[k]
            if p > remaining:
                break
            current.append(p)
            rec(remaining - p, k, current)
            current.pop()

    rec(target, 0, [])
    return results


def odd_partitions(n):
    """Partitions of n into odd parts {1,3,5,...}."""
    allowed = [p for p in range(1, n + 1) if p % 2 == 1]
    return partitions_from(n, allowed, extra_fixed=0)


def sc_partitions(n):
    """Cycle types contributing to SC(n) (THM-283):
       parts ≡ 2 mod 4, plus exactly one fixed point (part=1) iff n is odd."""
    allowed = [p for p in range(2, n + 1) if p % 4 == 2]   # {2,6,10,14,...}
    extra = 1 if n % 2 == 1 else 0
    return partitions_from(n, allowed, extra_fixed=extra)


# ---- the three sequences -----------------------------------------------------
def A000568(n):
    v = burnside_sum(odd_partitions(n), base=2)
    assert v.denominator == 1, f"A000568({n}) not integer: {v}"
    return v.numerator


def A_q(m, q):
    """q-deformed tournament count A(m,q) = Σ_{odd λ ⊢ m} q^c / z."""
    v = burnside_sum(odd_partitions(m), base=q)
    assert v.denominator == 1, f"A({m},{q}) not integer: {v}"
    return v.numerator


def SC(n):
    if n < 2:
        return 1  # SC(1)=1 by convention
    v = burnside_sum(sc_partitions(n), base=2)
    assert v.denominator == 1, f"SC({n}) not integer: {v}"
    return v.numerator


def V_merged(n):
    s = A000568(n) + SC(n)
    assert s % 2 == 0, f"A+SC odd at n={n}"
    return s // 2


# ---- main --------------------------------------------------------------------
if __name__ == "__main__":
    import time
    NMAX = 60
    print("=" * 72)
    print("EXACT BURNSIDE — A000568, SC, V_merged   (float-free, big-int only)")
    print("monad-researcher-2026-06-02-S560")
    print("=" * 72)

    print("\n[1] Cross-verification against known OEIS / canon values:")
    ok = True
    for n in range(1, 21):
        a, sc = A000568(n), SC(n)
        a_ok = (n not in KNOWN_A000568) or (a == KNOWN_A000568[n])
        sc_ok = (n not in KNOWN_SC) or (sc == KNOWN_SC[n])
        ok &= a_ok and sc_ok
        flag = "" if (a_ok and sc_ok) else "   <-- MISMATCH!"
        print(f"  n={n:2d}  A={a}   SC={sc}{flag}")
    print(f"\n  All known A000568 & SC values reproduced: {ok}")

    print("\n[2] Identity (I) V_merged = (A000568 + SC)/2 is an exact integer:")
    vm_ok = True
    vals = []
    for n in range(3, 21):
        vm = V_merged(n)
        vals.append(vm)
    print(f"  V_merged(3..20) = {vals}")
    print(f"  (all integer divisions exact — asserted in V_merged)")

    print("\n[3] Identity (II)  SC(2m) = A(m, 4)  [MISTAKE-049 corrected framing]:")
    ii_ok = True
    for m in range(1, 26):
        lhs = SC(2 * m)
        rhs = A_q(m, 4)
        match = (lhs == rhs)
        ii_ok &= match
        if m <= 8 or not match:
            print(f"  m={m:2d}:  SC({2*m}) = {lhs}   A({m},4) = {rhs}   "
                  f"[{'OK' if match else 'MISMATCH!'}]")
    print(f"  SC(2m) == A(m,4) for m=1..25: {ii_ok}")
    # sanity: A(m,2) must equal A000568(m)
    aq_ok = all(A_q(m, 2) == A000568(m) for m in range(1, 21))
    print(f"  Sanity A(m,2) == A000568(m) for m=1..20: {aq_ok}")

    print("\n[4] EXTENSION — exact values n=21..%d:" % NMAX)
    print(f"  {'n':>3} | {'A000568(n)':>12} ... | SC(n) ... | V_merged(n) ...")
    for n in range(21, NMAX + 1):
        t0 = time.time()
        a, sc, vm = A000568(n), SC(n), V_merged(n)
        dt = time.time() - t0
        print(f"  n={n}:")
        print(f"    A000568  = {a}")
        print(f"    SC       = {sc}")
        print(f"    V_merged = {vm}    ({dt:.3f}s)")

    print("\nDONE.  Every value asserted exact (Fraction denominator == 1).")
