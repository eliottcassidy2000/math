#!/usr/bin/env python3
"""THM-3001 -- Newton-circuit reversal involution and the two-end curvature law.

Universe.  N(n)=sum_{i=0}^d a_i n^i with every a_i > 0.  Normalized Newton
coefficients and ratios as in THM-2997 (1):

    h_k = a_{d-k}/(a_d C(d,k)),   h_0 = 1,
    R_k = h_k^2/(h_{k-1} h_{k+1}),   1 <= k <= d-1.

Coefficient reversal:  N*(n) = n^d N(1/n),  a*_i = a_{d-i}.

CLAIMS PROVED / VERIFIED HERE
  I.   (involution)  h*_k = (a_d/a_0) h_{d-k}, hence  R*_k = R_{d-k}  EXACTLY,
       for every 1 <= k <= d-1 and every positive-coefficient N.  Equivalently
       Delta^3 (log h*)_j = -Delta^3 (log h)_{d-3-j}: reversal negates every
       Newton circuit.  Unconditional, no asymptotics.
  II.  (class no-go)  If a hypothesis class H is closed under reversal and some
       N in H has a non-constant ratio sequence, then H cannot imply the global
       no-return property "R_k >= R_{k-1} for every 2 <= k <= d-1".  One line.
       Positive coefficients, positive real-rootedness, PF-infinity, Hurwitz
       stability and strict ULC are ALL reversal-closed -- checked here -- so
       none of them, in any combination, can prove no-return.
  III. (two-end curvature law)  Combining I with THM-3000: for fixed k,
         log(R_k/R_{k-1})           =  C(N)  d^{-2} + O(d^{-3}),
         log(R_{d-k+1}/R_{d-k})     = -C(N*) d^{-2} + O(d^{-3}),
       where C(N) = (3 kappa_2^2 - 2 kappa_1 kappa_3)/kappa_1^4 is THM-3000's
       cumulant curvature of the root measure and C(N*) is the same functional
       of the RECIPROCAL root measure.  So the two ends of the Newton circuit
       are governed by two explicit cumulant numbers, and C(N)>0, C(N*)>0
       forces an interior maximum.
  IV.  (Newton-ratio-palindromic class)  R_k is invariant under root scaling,
       so if the root multiset is closed under reciprocal UP TO A SCALAR then
       R_k = R_{d-k} exactly.  Every balanced two-cluster (n+a)^m (n+b)^m is in
       this class: its ratio sequence is exactly symmetric about k = m, so it
       rises then falls with an exact turn at the midpoint.

Relation to THM-2991.  THM-2991 is STRICTLY STRONGER in one direction: it adds
an arbitrarily long improving LEADING PREFIX (not a reversal-closed hypothesis)
and gets a return strictly BELOW R_1.  Claim II does not reach that; it explains
instead why the PF-infinity / Hurwitz / ULC decorations in THM-2991 are
automatically inert, and why the construction had to be one-sided.

Reproduce:  python3 04-computation/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.py
"""

from fractions import Fraction as Fr
from math import comb, log
import itertools


# --------------------------------------------------------------------------
def coeffs_from_roots(roots):
    """N = prod (n + r); returns a[i] = coefficient of n^i, i = 0..d."""
    d = len(roots)
    e = [Fr(1)] + [Fr(0)] * d
    for r in roots:
        for k in range(d, 0, -1):
            e[k] = e[k] + r * e[k - 1]
    return [e[d - i] for i in range(d + 1)]


def hseq(a):
    d = len(a) - 1
    return [Fr(a[d - k]) / (Fr(a[d]) * comb(d, k)) for k in range(d + 1)]


def ratios(a):
    h = hseq(a)
    d = len(a) - 1
    return [None] + [h[k] ** 2 / (h[k - 1] * h[k + 1]) for k in range(1, d)]


def reverse(a):
    return list(reversed(a))


# --------------------------------------------------------------------------
# I. the involution
# --------------------------------------------------------------------------
def claim_I():
    print("=" * 74)
    print("I. REVERSAL INVOLUTION:  h*_k = (a_d/a_0) h_{d-k}  and  R*_k = R_{d-k}")
    print("=" * 74)
    cases = [
        ("(n+1)^4(n+2)^4", [Fr(1)] * 4 + [Fr(2)] * 4),
        ("distinct 1,2,5,7,11", [Fr(1), Fr(2), Fr(5), Fr(7), Fr(11)]),
        ("rational mix", [Fr(1, 3), Fr(2), Fr(5, 7), Fr(9), Fr(1), Fr(4)]),
        ("THM-2991 control (1,1,3,20,20)", [Fr(1), Fr(1), Fr(3), Fr(20), Fr(20)]),
    ]
    ok = True
    # also a non-real-rooted positive-coefficient control
    cases.append(("non-real-rooted n^4+n^3+3n^2+n+1", None))
    for name, roots in cases:
        if roots is None:
            a = [Fr(1), Fr(1), Fr(3), Fr(1), Fr(1)]  # a_0..a_4
        else:
            a = coeffs_from_roots(roots)
        d = len(a) - 1
        h, hs = hseq(a), hseq(reverse(a))
        c = Fr(a[d], a[0])
        hok = all(hs[k] == c * h[d - k] for k in range(d + 1))
        R, Rs = ratios(a), ratios(reverse(a))
        Rok = all(Rs[k] == R[d - k] for k in range(1, d))
        ok &= hok and Rok
        print(f"  {name:34s} d={d}  h*_k=(a_d/a_0)h_(d-k): {hok}   R*_k=R_(d-k): {Rok}")
    print("  VERDICT I:", "EXACT INVOLUTION" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# II. reversal-closure of the standard hypothesis classes
# --------------------------------------------------------------------------
def is_pf_infinity_from_roots(roots):
    """prod (n+r), r>0, is PF-infinity (Aissen-Schoenberg-Whitney-Edrei)."""
    return all(r > 0 for r in roots)


def toeplitz_minors_nonneg(a, size=3):
    """Direct total-nonnegativity spot check of the Toeplitz matrix of a."""
    n = len(a)

    def A(i, j):
        k = i - j
        return Fr(a[k]) if 0 <= k < n else Fr(0)

    for s in range(1, size + 1):
        for rows in itertools.combinations(range(n + s), s):
            for cols in itertools.combinations(range(n + s), s):
                M = [[A(i, j) for j in cols] for i in rows]
                # exact Laplace for s<=3
                if s == 1:
                    det = M[0][0]
                elif s == 2:
                    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
                else:
                    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
                if det < 0:
                    return False
    return True


def strict_ulc(a):
    h = hseq(a)
    d = len(a) - 1
    return all(h[k] ** 2 > h[k - 1] * h[k + 1] for k in range(1, d))


def claim_II():
    print()
    print("=" * 74)
    print("II. REVERSAL-CLOSURE OF EVERY STANDARD HYPOTHESIS CLASS")
    print("=" * 74)
    print("  reversal maps roots r -> 1/r, so:")
    print("   * positive coefficients      : reversal permutes the coefficients      -> closed")
    print("   * all roots real and negative: 1/r stays real and negative             -> closed")
    print("   * PF-infinity                : ASWE = product of (1+r t), r>0; 1/r>0   -> closed")
    print("   * Hurwitz stable             : Re(root)<0 <=> Re(1/root)<0             -> closed")
    print("   * strict ULC                 : h*_k = c h_(d-k), log-concavity is      -> closed")
    print("                                  reversal-symmetric")
    print()
    tests = [
        ("(n+1)^4(n+2)^4", [Fr(1)] * 4 + [Fr(2)] * 4),
        ("(n+1)^2(n+3)(n+20)^2", [Fr(1), Fr(1), Fr(3), Fr(20), Fr(20)]),
        ("(n+1)(n+2)(n+4)(n+8)", [Fr(1), Fr(2), Fr(4), Fr(8)]),
    ]
    ok = True
    for name, roots in tests:
        a = coeffs_from_roots(roots)
        ar = reverse(a)
        rr = sorted(Fr(1) / r for r in roots)
        ar2 = coeffs_from_roots(rr)
        scal = Fr(ar[len(ar) - 1], ar2[len(ar2) - 1])
        same = all(ar[i] == scal * ar2[i] for i in range(len(ar)))
        pf = toeplitz_minors_nonneg(a) and toeplitz_minors_nonneg(ar)
        ulc = strict_ulc(a) and strict_ulc(ar)
        ok &= same and pf and ulc
        print(f"  {name:24s}: reversal == reciprocal-root poly (up to scalar) {same};"
              f"  PF3 both {pf};  strict ULC both {ulc}")
    print()
    print("  ONE-LINE NO-GO.  Let H be reversal-closed and suppose every N in H had")
    print("  R_k >= R_(k-1) for all k.  For N in H also N* in H, so R*_k >= R*_(k-1),")
    print("  i.e. R_(d-k) >= R_(d-k+1), i.e. R_j >= R_(j+1).  Both give R_j = R_(j+1)")
    print("  for every j: the ratio sequence is CONSTANT.  Any H containing one member")
    print("  with a non-constant ratio sequence therefore cannot imply no-return.")
    print("  VERDICT II:", "ALL LISTED CLASSES ARE REVERSAL-CLOSED" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# III. two-end curvature law
# --------------------------------------------------------------------------
def curvature(roots):
    d = len(roots)
    m1 = sum(roots) / d
    m2 = sum(r * r for r in roots) / d
    m3 = sum(r ** 3 for r in roots) / d
    return 3 * m2 ** 2 / m1 ** 4 - 2 * m3 / m1 ** 3 - 1


def claim_III():
    print()
    print("=" * 74)
    print("III. TWO-END CURVATURE LAW  (bottom = +C(N), top = -C(N*))")
    print("=" * 74)
    fams = {
        "two-cluster 1,2 balanced": lambda d: [Fr(1)] * (d // 2) + [Fr(2)] * (d - d // 2),
        "two-cluster 1,2 (1:3 split)": lambda d: [Fr(1)] * (d // 4) + [Fr(2)] * (d - d // 4),
        "three-cluster 1,3,20": lambda d: ([Fr(1)] * (d // 3) + [Fr(3)] * (d // 3)
                                           + [Fr(20)] * (d - 2 * (d // 3))),
        "geometric 2^{-(i mod 6)}": lambda d: [Fr(1, 2 ** (i % 6)) for i in range(d)],
    }
    ok = True
    for name, gen in fams.items():
        print(f"\n  -- {name}")
        for d in (120, 240, 480):
            r = gen(d)
            rr = [Fr(1) / t for t in r]
            CN, CNs = curvature(r), curvature(rr)
            a = coeffs_from_roots(r)
            h = hseq(a)

            def circ(k):   # log(R_k/R_(k-1)) = -Delta^3(log h)_(k-2)
                v = h[k + 1] * h[k - 1] ** 3 / (h[k] ** 3 * h[k - 2])
                return -log(float(v))
            bot = d * d * circ(2)
            top = d * d * circ(d - 1)
            good = (abs(bot - float(CN)) < 0.25 * max(1.0, abs(float(CN)))
                    and abs(top + float(CNs)) < 0.25 * max(1.0, abs(float(CNs))))
            ok &= good
            print(f"     d={d:4d}  C(N)={float(CN):12.6f} vs d^2*bottom={bot:12.6f} | "
                  f"-C(N*)={-float(CNs):12.6f} vs d^2*top={top:12.6f}  {'OK' if good else 'X'}")
    print("\n  VERDICT III:", "TWO-END LAW CONFIRMED" if ok else "MISMATCH")
    return ok


# --------------------------------------------------------------------------
# IV. Newton-ratio-palindromic class
# --------------------------------------------------------------------------
def claim_IV():
    print()
    print("=" * 74)
    print("IV. RECIPROCAL-CLOSED-UP-TO-SCALAR ROOTS  =>  R_k = R_(d-k) EXACTLY")
    print("=" * 74)
    print("  R_k is invariant under r -> lambda r (h_k -> lambda^k h_k cancels), so if")
    print("  the root multiset equals lambda * {1/r} for some lambda>0 then R*_k=R_k,")
    print("  and with I this gives R_k = R_(d-k).")
    ok = True
    for (a_, b_, m) in [(1, 2, 4), (1, 3, 6), (2, 5, 8), (1, 10, 10), (3, 7, 5)]:
        roots = [Fr(a_)] * m + [Fr(b_)] * m
        d = 2 * m
        R = ratios(coeffs_from_roots(roots))
        sym = all(R[k] == R[d - k] for k in range(1, d))
        ups = [k for k in range(2, d) if R[k] > R[k - 1]]
        downs = [k for k in range(2, d) if R[k] < R[k - 1]]
        turn = min(downs) if downs else None
        ok &= sym and (turn == m + 1)
        print(f"  (n+{a_})^{m}(n+{b_})^{m}: R_k=R_(d-k) exactly {sym};  rises k<= {max(ups)};"
              f"  first fall k={turn} (= m+1 = {m + 1});  monotone? {not downs}")
    print("  => the simplest non-degenerate real-rooted positive family already has")
    print("     an EXACT interior maximum: full no-return is impossible on it.")
    print("  VERDICT IV:", "EXACT PALINDROMY" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------
# V. the two-number sign classifier: proved necessary condition + census
# --------------------------------------------------------------------------
def shape_of(roots):
    d = len(roots)
    R = ratios(coeffs_from_roots(roots))
    up = all(R[k] > R[k - 1] for k in range(2, d))
    dn = all(R[k] < R[k - 1] for k in range(2, d))
    if up:
        return "INCREASING"
    if dn:
        return "DECREASING"
    interior_max = R[2] > R[1] and R[d - 1] < R[d - 2]
    interior_min = R[2] < R[1] and R[d - 1] > R[d - 2]
    return "INTERIOR-MAX" if interior_max else ("INTERIOR-MIN" if interior_min else "MIXED")


def predicted(roots):
    c1 = curvature(roots)
    c2 = curvature([Fr(1) / r for r in roots])
    if c1 > 0 and c2 < 0:
        return "INCREASING"
    if c1 < 0 and c2 > 0:
        return "DECREASING"
    if c1 > 0 and c2 > 0:
        return "INTERIOR-MAX"
    if c1 < 0 and c2 < 0:
        return "INTERIOR-MIN"
    return "DEGENERATE"


def claim_V():
    print()
    print("=" * 74)
    print("V. TWO-NUMBER SIGN CLASSIFIER  (proved necessary condition + census)")
    print("=" * 74)
    print("  PROVED (from I + III): asymptotic global no-return REQUIRES")
    print("        C(N) >= 0 >= C(N*).")
    print("  OBSERVED (census below, NOT proved): on these families the pair of")
    print("  signs (C(N), -C(N*)) also determines the global shape.")
    print()
    d = 60
    tests = []
    for num, den in [(1, 6), (1, 4), (1, 3), (2, 5), (1, 2), (3, 5), (2, 3), (3, 4), (5, 6)]:
        for b_ in (2, 3, 5, 10):
            m = d * num // den
            tests.append((f"1^{num}/{den} x {b_}^rest", [Fr(1)] * m + [Fr(b_)] * (d - m)))
    for q in (Fr(1, 2), Fr(2, 3), Fr(9, 10)):
        tests.append((f"geometric q={q}", [q ** i for i in range(d)]))
    for a_, b_, c_ in [(1, 3, 20), (1, 2, 4), (1, 5, 25)]:
        tests.append((f"three-cluster {a_},{b_},{c_}",
                      [Fr(a_)] * (d // 3) + [Fr(b_)] * (d // 3) + [Fr(c_)] * (d - 2 * (d // 3))))
    agree = 0
    nec_ok = True
    rows = []
    for name, roots in tests:
        s, p = shape_of(roots), predicted(roots)
        agree += (s == p)
        if s == "INCREASING":
            nec_ok &= (curvature(roots) >= 0 >= curvature([Fr(1) / r for r in roots]))
        rows.append((name, s, p, float(curvature(roots)),
                     float(curvature([Fr(1) / r for r in roots]))))
    for name, s, p, c1, c2 in rows[:14]:
        print(f"   {name:24s} C(N)={c1:+9.5f} C(N*)={c2:+9.5f}  actual={s:12s} "
              f"predicted={p:12s} {'OK' if s == p else 'MISMATCH'}")
    print(f"   ... {len(rows)} families total")
    print()
    print(f"  classifier agreement: {agree}/{len(rows)}")
    print(f"  proved necessary condition C(N)>=0>=C(N*) held on every INCREASING case: {nec_ok}")
    print("  VERDICT V:", "CLASSIFIER MATCHES ON EVERY CENSUS FAMILY"
          if agree == len(rows) and nec_ok else
          f"PARTIAL ({agree}/{len(rows)}); necessary condition {nec_ok}")
    return nec_ok, agree == len(rows)


def main():
    a = claim_I()
    b = claim_II()
    c = claim_III()
    e = claim_IV()
    nec, cls = claim_V()
    print()
    print("=" * 74)
    print("SUMMARY  I=%s  II=%s  III=%s  IV=%s  V(necessary)=%s  V(classifier census)=%s"
          % (a, b, c, e, nec, cls))
    print("=" * 74)


if __name__ == "__main__":
    main()
