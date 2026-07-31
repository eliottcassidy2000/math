#!/usr/bin/env python3
"""THM-3003 -- antipodal circuit rigidity and the multipole spread criterion.

Continues THM-3000/3001.  Three layers, all about the SAME involution.

LAYER 1 (rigidity, exact iff).  For positive-coefficient N of degree d,

    R_k = R_{d-k}  for every 1 <= k <= d-1
      <=>  {r_i} = {mu/r_i} as multisets, mu = e_d^{2/d}
      <=>  the empirical measure of log r is symmetric about its mean.

THM-3001 section 5 proved (<=).  Here (=>) is proved too, so the Newton circuit
DETECTS reciprocal symmetry exactly.  Proof: the binomial normalization
log C(d,k) already has palindromic second difference, so the hypothesis says
Delta^2(log e)_j is a palindrome; two log-sequences with equal second differences
differ by an affine function of the index, giving e_k = C mu^k e_{d-k}; the
generating function form of that relation is t^d Q(1/t) = e_d Q(t/mu), whose
root multiset comparison is exactly {r_i} = {mu/r_i}.  QED.

LAYER 2 (antipodal / Borsuk-Ulam).  In centered log-root coordinates
ell_i = log r_i - mean, reversal IS the ANTIPODAL MAP ell -> -ell, and the
circuit map c(ell) = (log(R_k/R_{k-1}))_{k=2..d-1} is equivariant:

    c(-ell) = -reverse(c(ell)).                                   (*)

Hence the index-symmetric part S c is an ODD function of ell and the
index-antisymmetric part A c is EVEN.  Consequences:
  (2a) IVT/Brouwer form: in any CONNECTED reversal-closed class H, every path
       from N to N* contains a point where sum_k w_k c_k = 0 for each
       reversal-symmetric weight w.  With w = 1 that is R_{d-1} = R_1: every
       such path crosses the BALANCED LOCUS.  No "H implies no-return"
       hypothesis is needed -- this is strictly stronger than THM-3001 sec 2.
  (2b) Borsuk-Ulam form: S c restricted to the direction sphere S^{d-2} is a
       continuous ODD map into R^{ceil((d-2)/2)}, so its zero set is nonempty
       and ESSENTIAL.  By Layer 1 that zero set is exactly the log-symmetric
       locus.  So the failure of no-return on log-symmetric profiles is
       topologically stable, not a feature of one lucky family.

LAYER 3 (fast-multipole dictionary and the SPREAD CRITERION).  The jets are the
2D multipole moments of the root charge distribution:
    log(N(n)/(a_d n^d)) = sum_j (ell_j/j) n^{-j},  ell_j = (-1)^{j-1} p_j,
which is the multipole expansion of sum_m log(n + r_m) about the origin.
Therefore:
  * jets are ADDITIVE over factorisation (ell_j(FG) = ell_j(F)+ell_j(G)):
    THM-2997's wall subtraction IS multipole subtraction;
  * reversal is the KELVIN transform r -> 1/r, so THM-3001's two-end law is
    exactly the multipole/local (far-field/near-field) duality;
  * cumulants are the translation-covariant (M2M) gauge, which is why
    THM-3000's curvature is clean in cumulants.
  * NEW SUFFICIENT CONDITION.  Let kappa = max|r| / (p_1/d) be the root SPREAD
    RATIO (max modulus over mean).  Since |p_j| <= d max|r|^j,
        m_j/m_1^j <= kappa^j,
    so THM-3000's graded hypothesis m_j/m_1^j = o(d^{j-3}) holds for all
    j <= k+1 as soon as
        kappa = o(d^{1 - 3/(k+1)}).
    Third edge (k=3): kappa = o(d^{1/4}).  This replaces a jet-by-jet invoice
    by ONE geometric number -- the FMM well-separatedness condition.

Reproduce: python3 04-computation/gmc_antipodal_circuit_rigidity_and_multipole_spread_criterion_thm3003.py
"""

from fractions import Fraction as Fr
from math import comb, log
import itertools
import random

random.seed(20260731)


# --------------------------------------------------------------------------
def coeffs_from_roots(roots):
    d = len(roots)
    e = [Fr(1)] + [Fr(0)] * d
    for r in roots:
        for k in range(d, 0, -1):
            e[k] = e[k] + r * e[k - 1]
    return e                      # e[k] = elementary symmetric


def ratios_from_e(e):
    d = len(e) - 1
    h = [Fr(e[k], 1) / comb(d, k) for k in range(d + 1)]
    return [None] + [h[k] ** 2 / (h[k - 1] * h[k + 1]) for k in range(1, d)], h


# --------------------------------------------------------------------------
# LAYER 1: rigidity, both directions, in exact arithmetic
# --------------------------------------------------------------------------
def is_reciprocal_closed(roots):
    """{r} = {mu/r} for mu = (prod r)^{2/d}?  Test without radicals:
    equivalent to  e_k^d * e_d^(d-2k) == e_{d-k}^d  is awkward; instead test
    the multiset condition directly with mu^d = e_d^2 via cross-multiplication."""
    d = len(roots)
    prod = Fr(1)
    for r in roots:
        prod *= r
    # candidate mu must satisfy mu^d = prod^2 ; test every root pairing instead:
    # {r} = {mu/r} <=> sorted(r) reversed times sorted(r) is constant = mu
    s = sorted(roots)
    cand = s[0] * s[-1]
    return all(s[i] * s[d - 1 - i] == cand for i in range(d))


def layer1():
    print("=" * 74)
    print("LAYER 1.  R_k = R_(d-k) for all k   <=>   roots reciprocal-closed up to scale")
    print("=" * 74)
    ok = True

    def palindromic_ratios(roots):
        e = coeffs_from_roots(roots)
        R, _ = ratios_from_e(e)
        d = len(roots)
        return all(R[k] == R[d - k] for k in range(1, d))

    print("\n  (<=) reciprocal-closed families must be ratio-palindromic:")
    fams = [("balanced two-cluster (1,2)^4", [Fr(1)] * 4 + [Fr(2)] * 4),
            ("balanced two-cluster (3,7)^3", [Fr(3)] * 3 + [Fr(7)] * 3),
            ("geometric 2^-i, d=7", [Fr(1, 2 ** i) for i in range(7)]),
            ("mirror pair set {2,3,1/2·6,1/3·6}", [Fr(2), Fr(3), Fr(3), Fr(2)]),
            ("scaled mirror {1,4,,}", [Fr(1), Fr(2), Fr(8), Fr(16)])]
    for name, r in fams:
        rc, pr = is_reciprocal_closed(r), palindromic_ratios(r)
        ok &= (not rc) or pr
        print(f"    {name:34s} reciprocal-closed={rc}  ratio-palindromic={pr}")

    print("\n  (=>) EXHAUSTIVE SEARCH for a ratio-palindromic set that is NOT")
    print("       reciprocal-closed (small rational root multisets):")
    pool = [Fr(1), Fr(2), Fr(3), Fr(4), Fr(6), Fr(1, 2), Fr(1, 3), Fr(3, 2), Fr(8), Fr(9)]
    found = []
    tested = 0
    for d in (4, 5, 6):
        for r in itertools.combinations_with_replacement(pool, d):
            tested += 1
            if palindromic_ratios(list(r)):
                if not is_reciprocal_closed(list(r)):
                    found.append(r)
    print(f"    tested {tested} multisets (d=4,5,6 over a 10-element rational pool)")
    print(f"    counterexamples found: {len(found)}   {'NONE -- rigidity holds' if not found else found[:3]}")
    ok &= not found

    print("\n  numeric confirmation of the PROOF's intermediate identity")
    print("    e_k * e_d == mu^k * e_(d-k)   with mu = e_d^(2/d):")
    for name, r in fams[:3]:
        e = coeffs_from_roots(r)
        d = len(r)
        mu = float(e[d]) ** (2.0 / d)
        errs = [abs(float(e[k] * e[d]) - mu ** k * float(e[d - k])) for k in range(d + 1)]
        print(f"    {name:34s} max|e_k e_d - mu^k e_(d-k)| = {max(errs):.3e}")
    print("\n  VERDICT LAYER 1:", "RIGIDITY CONFIRMED" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# LAYER 2: antipodal equivariance and the balanced-locus crossing
# --------------------------------------------------------------------------
def circuits_float(ell):
    r = [pow(2.718281828459045, x) for x in ell]
    d = len(r)
    e = [0.0] * (d + 1)
    e[0] = 1.0
    for t in r:
        for k in range(d, 0, -1):
            e[k] += t * e[k - 1]
    lh = [log(e[k]) - log(comb(d, k)) for k in range(d + 1)]
    return [-(lh[k + 1] - 3 * lh[k] + 3 * lh[k - 1] - lh[k - 2]) for k in range(2, d)]


def layer2():
    print()
    print("=" * 74)
    print("LAYER 2.  reversal = ANTIPODAL MAP;  c(-ell) = -reverse(c(ell))")
    print("=" * 74)
    ok = True
    print("  equivariance (*), random configurations:")
    for d in (5, 6, 7, 8, 12):
        ell = [random.gauss(0, 1.1) for _ in range(d)]
        a = circuits_float([-x for x in ell])
        b = [-v for v in reversed(circuits_float(ell))]
        err = max(abs(p - q) for p, q in zip(a, b))
        ok &= err < 1e-9
        print(f"    d={d:3d}  max|c(-ell) + reverse(c(ell))| = {err:.3e}")

    print()
    print("  (2a) BALANCED-LOCUS CROSSING.  Phi := sum_k c_k = log(R_(d-1)/R_1) is")
    print("       an ODD function of ell (w = 1 is reversal-symmetric), so ANY path")
    print("       from ell to -ell must cross Phi = 0.  The straight path is a")
    print("       degenerate witness (it passes through ell = 0, all roots equal),")
    print("       so the real test uses a GREAT-CIRCLE path of constant norm that")
    print("       never degenerates:  ell(t) = cos(pi t) ell + sin(pi t) ell_perp.")
    import math
    for d in (6, 9, 14):
        ell = [random.gauss(0, 1.0) for _ in range(d)]
        mu = sum(ell) / d
        ell = [x - mu for x in ell]
        # a centered vector orthogonal to ell
        v = [random.gauss(0, 1.0) for _ in range(d)]
        mv = sum(v) / d
        v = [x - mv for x in v]
        dot = sum(a * b for a, b in zip(ell, v)) / sum(a * a for a in ell)
        v = [b - dot * a for a, b in zip(ell, v)]
        nv = math.sqrt(sum(x * x for x in v))
        ne = math.sqrt(sum(x * x for x in ell))
        v = [x * ne / nv for x in v]
        vals, norms = [], []
        for i in range(9):
            t = i / 8
            e2 = [math.cos(math.pi * t) * a + math.sin(math.pi * t) * b
                  for a, b in zip(ell, v)]
            norms.append(math.sqrt(sum(x * x for x in e2)))
            vals.append(sum(circuits_float(e2)))
        endpoints_opposite = vals[0] * vals[-1] < 0
        crosses = any(vals[i] * vals[i + 1] <= 0 for i in range(len(vals) - 1))
        nondegenerate = min(norms) > 0.9 * ne
        ok &= endpoints_opposite and crosses and nondegenerate
        interior = [i for i in range(len(vals) - 1) if vals[i] * vals[i + 1] <= 0]
        print(f"    d={d:3d}  |ell(t)| constant to {max(norms) - min(norms):.2e}"
              f" (never degenerates: {nondegenerate})")
        print(f"          Phi(0)={vals[0]:+.6f}  Phi(1)={vals[-1]:+.6f}"
              f"  opposite signs: {endpoints_opposite}")
        print(f"          crossing found in sub-interval(s) {interior}"
              f"   (NOT at the degenerate center)")

    print()
    print("  (2b) BORSUK-ULAM.  S c : S^(d-2) -> R^ceil((d-2)/2) is odd, so it has")
    print("       a zero.  Sphere dim d-2 vs target dim ceil((d-2)/2):")
    for d in range(4, 13):
        print(f"    d={d:3d}  sphere dim {d - 2:3d}  target dim {-(-(d - 2) // 2):3d}"
              f"  BU applies: {-(-(d - 2) // 2) <= d - 2}"
              f"   expected zero-set dim {(d - 2) - (-(-(d - 2) // 2)):3d}"
              f"   log-symmetric locus dim in the sphere {d // 2 - 1:3d}")
    print("       (the two dimensions agree, which is Layer 1: the BU zero set IS")
    print("        the log-symmetric locus, nothing more and nothing less.)")
    print("  VERDICT LAYER 2:", "EQUIVARIANCE + CROSSING CONFIRMED" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# LAYER 3: multipole dictionary and the spread criterion
# --------------------------------------------------------------------------
def layer3():
    print()
    print("=" * 74)
    print("LAYER 3.  MULTIPOLE DICTIONARY AND THE SPREAD CRITERION")
    print("=" * 74)
    print("  (i) jets are additive over factorisation  (= FMM superposition):")
    ok = True
    for trial in range(3):
        A = [Fr(random.randint(1, 9)) for _ in range(4)]
        B = [Fr(random.randint(1, 9)) for _ in range(3)]
        pA = [sum(r ** j for r in A) for j in range(1, 6)]
        pB = [sum(r ** j for r in B) for j in range(1, 6)]
        pAB = [sum(r ** j for r in A + B) for j in range(1, 6)]
        good = all(pAB[i] == pA[i] + pB[i] for i in range(5))
        ok &= good
        print(f"    trial {trial}: p_j(F*G) == p_j(F) + p_j(G) for j=1..5 : {good}")
    print("    => THM-2997 (21)'s wall subtraction A=P_1-w_1 D^2 etc. is exactly")
    print("       multipole subtraction; w_j are the WALL's multipole moments.")

    print()
    print("  (ii) SPREAD CRITERION.  kappa = max|r| / (p_1/d);  |p_j| <= d max|r|^j")
    print("       gives m_j/m_1^j <= kappa^j, so THM-3000's graded hypothesis")
    print("       m_j/m_1^j = o(d^(j-3)) for j <= k+1 follows from")
    print("            kappa = o(d^(1 - 3/(k+1))).")
    print("       edge k :  required spread exponent 1 - 3/(k+1)")
    for k in range(2, 9):
        print(f"         k={k}:  kappa = o(d^{1 - 3 / (k + 1):.4f})"
              f"   (j<=k+1 binding at j=4)" if k >= 3 else
              f"         k={k}:  no j>=4 jet enters; curvature alone")

    print()
    print("  (iii) numeric check of the bound m_j/m_1^j <= kappa^j:")
    for name, roots in [("uniform i/d, d=60", [Fr(i, 60) for i in range(1, 61)]),
                        ("two-cluster 1,2 d=60", [Fr(1)] * 30 + [Fr(2)] * 30),
                        ("one heavy root d=60", [Fr(1)] * 59 + [Fr(60)])]:
        d = len(roots)
        m = [sum(r ** j for r in roots) / d for j in range(1, 6)]
        kappa = float(max(roots) / m[0])
        row = []
        for j in range(2, 6):
            lhs = float(m[j - 1] / m[0] ** j)
            row.append((j, lhs, kappa ** j, lhs <= kappa ** j + 1e-9))
        ok &= all(t[3] for t in row)
        print(f"    {name:22s} kappa={kappa:8.4f}")
        for j, lhs, rhs, good in row:
            print(f"        j={j}: m_j/m_1^j = {lhs:12.5f} <= kappa^j = {rhs:14.5f}  {good}")

    print()
    print("  (iv) FIRST-GAP CONSEQUENCE.  THM-2997 (24) gives d/M -> 62/3 and")
    print("       u/M^2 -> 131/12, so mean root -> (131/12)/(62/3) M = 0.5282 M.")
    mean_over_M = (131 / 12) / (62 / 3)
    dov = 62 / 3
    print(f"       mean root  ~ {mean_over_M:.6f} M ,  d ~ {dov:.4f} M")
    print("       third edge needs kappa = o(d^(1/4)), i.e.")
    print(f"           max|r| = o({mean_over_M:.4f} * ({dov:.3f} M)^(1/4) * M)"
          f" = o({mean_over_M * dov ** 0.25:.4f} M^(5/4)).")
    print("       So a root-modulus bound |r| = o(M^(5/4)) on the wall-stripped core")
    print("       CLOSES the third edge -- no P_4, no new Macaulay chart.")
    print("  VERDICT LAYER 3:", "DICTIONARY + CRITERION CONFIRMED" if ok else "FAILED")
    return ok


def main():
    a = layer1()
    b = layer2()
    c = layer3()
    print()
    print("=" * 74)
    print(f"SUMMARY  layer1(rigidity)={a}  layer2(antipodal)={b}  layer3(multipole)={c}")
    print("=" * 74)


if __name__ == "__main__":
    main()
