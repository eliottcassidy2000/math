# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-10-S127: the difference-primitive tight-family search.
#
# Context: asked to prove TightRigidity (mu=0 => dilated) for the difference-primitive case.
# The provable piece is the COLLAPSE: for a primitive family, "dilated" forces c=1, so the
# family is literally {1,...,13} (dilated_primitive_eq_range, kernel-pure in LRCTightRigidity.lean).
# The restricted rigidity PrimitiveTightRigidity (mu=0 => {1,...,13} for primitive v) is STILL the
# conjecture. This script tests it: is there any PRIMITIVE (gcd=1) mu=0 family that is NOT {1,...,13}?
# If found, PrimitiveTightRigidity would be REFUTED. None found in range => evidence FOR it.
#
# mu(S)=0 (measure-zero safe set) is decided exactly: the safe set is a union of closed intervals in
# the residue partition; it has measure 0 iff no open partition cell is fully safe. We test the
# midpoint of each maximal gap between breakpoints {k/(14 v) : ...} for safety.
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

def mu_is_zero(S):
    """True iff the safe set {t in [0,1] : all |v t - round| >= 1/14} has measure 0."""
    bps = {F(0), F(1)}
    for v in S:
        for k in range(v):
            bps.add(F(14 * k + 1, 14 * v))
            bps.add(F(14 * k + 13, 14 * v))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        if all(F(1, 14) <= (v * mid) % 1 <= F(13, 14) for v in S):
            return False   # an open cell is fully safe => mu > 0
    return True

def main():
    print("Difference-primitive tight-family search (kps-S127).")
    print("Testing PrimitiveTightRigidity: is any primitive mu=0 family != {1,...,13}?\n")
    ap = tuple(range(1, 14))
    assert mu_is_zero(list(ap)), "sanity: {1,...,13} must be tight"
    print(f"sanity: {{1,...,13}} is tight (mu=0): OK")
    found = []
    for N in range(14, 21):          # Vmax = N; N=13 is the AP itself
        for S in combinations(range(1, N + 1), 13):
            if max(S) != N:
                continue             # only NEW families (Vmax exactly N)
            if reduce(gcd, S) != 1:
                continue             # PRIMITIVE only (dilates c>=2 are non-primitive by construction)
            if mu_is_zero(list(S)):
                found.append(S)
                print(f"  PRIMITIVE mu=0, Vmax={N}, != AP: {S}")
    print(f"\nprimitive mu=0 families beyond {{1,...,13}} in Vmax<=20: {len(found)}")
    if not found:
        print("=> ZERO. Only the AP {1,...,13} is a primitive tight family in range.")
        print("=> evidence FOR PrimitiveTightRigidity (NOT a proof; it is >= LRC-hard).")

if __name__ == "__main__":
    main()
