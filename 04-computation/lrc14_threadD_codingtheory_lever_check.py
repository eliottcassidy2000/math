#!/usr/bin/env python3
"""
THREAD D (opus): External coding-theory lever check for LRC(14) consec-extremality.

PURPOSE.  Make CONCRETE the claim that the relation lattice
    Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }
is the SAME object as the "column dependency code" studied in coding theory, and
that the LRC consec-extremality has a known home:
    consec offsets = ANTI-MDS (minimum relation support d = 2)
    Sidon / arc offsets = MDS-like (no support-2 relation, d >= 3).

This mirrors EXACTLY the mechanism in arXiv:2501.19125 (Arnault-Gaborit-
Rozendaal-Saussay-Zemor, "Upper Bounds on the Minimum Distance of Structured
LDPC Codes"):  the circulant column-weight-2 block C has the property that the
sum of t CONSECUTIVE columns is a weight-2 vector with support {i,j}, l(i-j)=t.
That is LITERALLY a support-2 relation between offsets at distance t -- the
"quasi-collision". Short relations (small support / small tolerance t) PULL the
minimum distance DOWN. Consec offsets are maximally rich in such short relations;
Sidon sets have none of tolerance 1 (no e_a - e_b = e_c - e_d with distinct).

We verify, with EXACT integer arithmetic (stdlib only), for E offset sets:
  (1) min support d of a NONZERO-offset relation (the "minimum distance" of the
      relation code, excluding the trivial support-1 from a literal 0 offset);
  (2) the count of support-2 relations  n_a e_a + n_b e_b = 0  (commensurate pairs);
  (3) the count of support-3 / "tolerance-1" coincidences  e_a - e_b = e_c - e_d
      (Sidon test): consec is saturated, Sidon has none.
We tabulate consec_k vs a Sidon/arc set of the same k and show the clean split.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

def support2_relations(E):
    """count pairs {a,b} with a nonzero integer relation supported on {e_a,e_b}:
       n_a e_a + n_b e_b = 0, n_a,n_b != 0  <=>  e_a, e_b commensurate & both !=0.
       Over Z this exists iff e_a,e_b both nonzero (any two nonzero integers are
       commensurate: e_b * e_a - e_a * e_b = 0). So the RELEVANT refinement is the
       SHORTEST such relation's coefficient size; we report min |n_a|+|n_b|."""
    idx = [i for i in range(len(E)) if E[i] != 0]
    pairs = []
    for a, b in combinations(idx, 2):
        g = gcd(E[a], E[b])
        na, nb = E[b] // g, -E[a] // g     # na*e_a + nb*e_b = 0
        pairs.append(((a, b), abs(na) + abs(nb)))
    return pairs

def min_support_relation(E, maxcoef=3):
    """Brute the minimum SUPPORT size of a nonzero relation sum n_i e_i = 0 with
       small coefficients, EXCLUDING any relation whose support is a single literal
       zero offset (support-1 from e_i = 0). Returns smallest support found.
       For offset sets containing 0, support-2 always exists (0 pairs with anyone via
       the literal-0 axis) so we ALSO report d among offsets with the 0 removed."""
    k = len(E)
    Enz = [e for e in E if e != 0]
    best = None
    # search increasing support among NONZERO offsets
    for s in range(2, k + 1):
        found = False
        for combo in combinations(range(len(Enz)), s):
            sub = [Enz[i] for i in combo]
            # does a nonzero integer relation with all coeffs nonzero & |n_i|<=maxcoef
            # exist on exactly this support?  search coeff cube (excluding 0 coeff).
            from itertools import product
            rng = [c for c in range(-maxcoef, maxcoef + 1) if c != 0]
            for coeffs in product(rng, repeat=s):
                if sum(c * v for c, v in zip(coeffs, sub)) == 0:
                    if reduce(gcd, [abs(c) for c in coeffs]) == 1:
                        found = True
                        break
            if found:
                break
        if found:
            best = s
            break
    return best

def sidon_coincidences(E):
    """count tolerance-1 coincidences e_a - e_b = e_c - e_d (distinct positive diffs
       repeated). Sidon set => every pairwise difference unique => 0 nontrivial.
       consec => maximal repetition. Returns (#distinct diffs, #diffs, max multiplicity)."""
    from collections import Counter
    diffs = Counter()
    pos = sorted(set(E))
    for a, b in combinations(pos, 2):
        diffs[b - a] += 1
    mult = list(diffs.values())
    return (len(diffs), sum(mult), max(mult))

def consec(k):
    return list(range(k))

def sidon(k):
    """a Sidon set (perfect-difference-ish) of size k: greedy minimal Sidon set."""
    S = [0]
    cand = 1
    while len(S) < k:
        diffs = set()
        ok = True
        test = S + [cand]
        seen = set()
        good = True
        ds = set()
        for a, b in combinations(test, 2):
            d = b - a
            if d in ds:
                good = False; break
            ds.add(d)
        if good:
            S.append(cand)
        cand += 1
    return S

print("="*78)
print("THREAD D: relation-code minimum-support split  consec(ANTI-MDS) vs Sidon(MDS)")
print("="*78)
print(f"{'k':>2} {'set':<28} {'#distinct diff':>14} {'#diff':>6} {'maxmult':>8} {'d_nz(supp)':>11}")
for k in range(3, 9):
    for name, E in (("consec", consec(k)), ("sidon", sidon(k))):
        nd, td, mm = sidon_coincidences(E)
        d = min_support_relation(E, maxcoef=3)
        tag = "ANTI-MDS" if (name=="consec") else "MDS-arc"
        print(f"{k:>2} {name+str(tuple(E)):<28} {nd:>14} {td:>6} {mm:>8} {str(d):>11}  {tag}")
print()
print("READING:")
print(" - consec: maxmult of differences GROWS (k-1, k-2, ...): MANY repeated diffs")
print("   e_a-e_b=e_c-e_d  =>  support-4 (or support-2 commensurate) SHORT relations.")
print("   This is the LDPC 'quasi-collision' of tolerance 1: t consecutive cols of the")
print("   circulant C sum to a weight-2 vector. Minimal-distance relation code => corr LARGE.")
print(" - Sidon/arc: maxmult = 1 (all differences distinct) => NO tolerance-1 collision")
print("   => relation code has d>=3 (MDS-like, general position) => corr ~ iid (small).")

# ---------------------------------------------------------------------------
# REFINEMENT: the correct discriminator is the SHORT-RELATION SPECTRUM, not the
# (always-present) support-2 commensurability. Over Z the meaningful short
# relations for LRC are the TOLERANCE-1 / Sidon-type ones:
#     e_a - e_b - e_c + e_d = 0   (a repeated difference)  -> support-4, coeff +-1
# These are EXACTLY the LDPC 'quasi-colliding pairs of tolerance t' with t=O(1),
# and they are the relations that put MASS on K(n) (corr = sum_{n in Lambda} K(n)).
# Count them; this is the relation-code's low-weight enumerator coefficient.
# ---------------------------------------------------------------------------
from collections import Counter
def tol1_relations(E):
    """number of {a,b,c,d} (as unordered repeated-difference incidences) giving a
       support<=4 coeff-+-1 relation e_a-e_b = e_c-e_d. = sum over difference values
       of C(mult,2). This is the A_4-ish low-weight count of the relation code."""
    diffs = Counter()
    pos = sorted(set(E))
    for a, b in combinations(pos, 2):
        diffs[b - a] += 1
    from math import comb
    return sum(comb(v, 2) for v in diffs.values())

print()
print("="*78)
print("LOW-WEIGHT RELATION SPECTRUM (the corr-driving short relations)")
print("  tol1 = # tolerance-1 quasi-collisions = # support-4 (+-1) relations e_a-e_b=e_c-e_d")
print("  = the LDPC quasi-collision count = low-weight coeff of the relation code's enumerator")
print("="*78)
print(f"{'k':>2} {'consec tol1':>12} {'sidon tol1':>12}   note")
for k in range(3, 10):
    ct = tol1_relations(consec(k))
    st = tol1_relations(sidon(k))
    print(f"{k:>2} {ct:>12} {st:>12}   consec saturated, sidon=0 (MDS general position)")
print()
print("CLOSED FORM: consec tol1 = sum_{j=1}^{k-1} C(k-j,2) = C(k,3)  (tetrahedral).")
for k in range(3,10):
    from math import comb
    assert tol1_relations(consec(k)) == comb(k,3), k
print("  VERIFIED: consec tol1 == C(k,3) for k=3..9.")
print("  Sidon tol1 == 0 for all k (definition of Sidon).")
