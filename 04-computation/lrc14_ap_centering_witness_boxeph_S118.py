#!/usr/bin/env python3
"""
boxeph-2026-07-18-S118 — The inverse-scaling centering witness for AP loneliness.

CLAIM (new, closes the S117 spread case for n=12):
For a primitive AP  C = {a + d*k : k=0..11}  with gcd(a,d)=1 and d ODD, put
    q = 2a + 11d,     p = d^{-1} mod q.
Then t = p/q is a witness giving
    min_k || (a+d k) * t ||  =  (q-11)/(2q)  =  1/2 - 11/(2q).
Hence  M(C) >= 1/2 - 11/(2q),  which is > 1/13 for every q > 13, and = 1/13
only at q = 13  <=>  a=d=1  <=>  C = {1,...,12}.

Mechanism: (a+d k)*p = a p + k (mod q) since d p = 1 (mod q).  With q=2a+11d,
2 a p = -(11) (mod q), so s := a p mod q = (q-11)/2  (q odd since d odd, n-1=11 odd),
and the 12 residues {s,...,s+11} = {(q-11)/2, ..., (q+11)/2} are 12 consecutive
integers symmetric about q/2.  Min distance to {0,q} is (q-11)/2 at the endpoints.

For d EVEN with gcd(a,d)=1 (=> a odd): t=1/2 gives min_k ||.|| = 1/2 (all terms odd).

This script:
  (A) verifies the witness formula exactly (Fraction) over many (a,d),
  (B) computes the TRUE M(C) by exact maximization over rational t and checks
      witness <= true M, and that M > 1/13 for every AP except {1,...,12},
  (C) confirms {1,...,12} is the unique primitive tight AP (M=1/13),
  (D) checks the general-n extension q = 2a+(n-1)d for even n-1.
"""
from fractions import Fraction as F
from math import gcd


def frac_dist(x: F) -> F:
    """distance from rational x to nearest integer, exact."""
    r = x - (x.numerator // x.denominator)   # fractional part in [0,1)
    return min(r, 1 - r)


def min_over_speeds(speeds, t: F) -> F:
    return min(frac_dist(F(v) * t) for v in speeds)


def true_M(speeds):
    """Exact M = max_t min_i ||v_i t||.  For a finite speed set the maximizer is a
    rational t=j/Q; candidates: all j/Q with Q up to 2*max(speed)+1 (safe superset
    of the standard denominator bound for this min-max)."""
    Qmax = 2 * max(speeds) + 2
    best = F(0)
    bestt = F(0)
    for Q in range(1, Qmax + 1):
        for j in range(0, Q):
            t = F(j, Q)
            m = min_over_speeds(speeds, t)
            if m > best:
                best = m
                bestt = t
    return best, bestt


def witness_center(a, d, n=12):
    """Return (q, p, s, witness_min) for the centering witness, d odd."""
    q = 2 * a + (n - 1) * d
    p = pow(d, -1, q)                       # d^{-1} mod q  (needs gcd(d,q)=1)
    s = (a * p) % q
    residues = [(s + k) % q for k in range(n)]
    wit = min(F(min(r, q - r), q) for r in residues)
    return q, p, s, residues, wit


def report(msg):
    print(msg, flush=True)


report("=" * 74)
report("(A) EXACT witness verification, n=12, d odd, gcd(a,d)=1")
report("=" * 74)
report(f"{'a':>4} {'d':>4} {'q=2a+11d':>9} {'s=(q-11)/2?':>12} "
       f"{'witness':>10} {'(q-11)/2q':>12} {'>1/13?':>7}")
thirteen = F(1, 13)
allok = True
for d in range(1, 60, 2):            # d odd
    for a in range(1, 60):
        if gcd(a, d) != 1:
            continue
        q, p, s, res, wit = witness_center(a, d, 12)
        formula = F(q - 11, 2 * q)
        s_ok = (s == (q - 11) // 2) and (q % 2 == 1)
        match = (wit == formula)
        gt = wit > thirteen
        if not (s_ok and match):
            allok = False
            report(f"  MISMATCH a={a} d={d}: s={s} expected={(q-11)//2} "
                   f"q_odd={q%2==1} wit={wit} formula={formula}")
        # print a sample
        if a <= 3 and d <= 9:
            report(f"{a:>4} {d:>4} {q:>9} {str(s_ok):>12} "
                   f"{str(wit):>10} {str(formula):>12} {str(gt):>7}")
report(f"\n  witness == (q-11)/(2q) AND s==(q-11)/2 for ALL tested: {allok}")

report("")
report("=" * 74)
report("(B)+(C) TRUE M vs witness; uniqueness of the tight AP  (small ranges)")
report("=" * 74)
report(f"{'a':>3} {'d':>3} {'AP':>26} {'trueM':>9} {'witness':>9} "
       f"{'wit<=M':>7} {'M>1/13':>7} {'tight?':>7}")
tights = []
for d in list(range(1, 12)) + [17, 19, 23, 41]:
    for a in range(1, 12):
        if gcd(a, d) != 1:
            continue
        speeds = [a + d * k for k in range(12)]
        if max(speeds) > 240:          # keep true_M search bounded
            continue
        M, tstar = true_M(speeds)
        # witness (choose parity-appropriate)
        if d % 2 == 1:
            q, p, s, res, wit = witness_center(a, d, 12)
        else:
            wit = min_over_speeds(speeds, F(1, 2))   # t=1/2, all-odd terms
        wit_le = wit <= M
        gt13 = M > thirteen
        tight = (M == thirteen)
        if tight:
            tights.append((a, d, speeds))
        if (a <= 4 and d <= 7) or tight or d in (17, 19, 23, 41) and a <= 3:
            apstr = "{" + ",".join(map(str, speeds[:3])) + ",..,%d}" % speeds[-1]
            report(f"{a:>3} {d:>3} {apstr:>26} {str(M):>9} {str(wit):>9} "
                   f"{str(wit_le):>7} {str(gt13):>7} {str(tight):>7}")

report(f"\n  TIGHT primitive APs found (M==1/13): "
       f"{[ (a,d) for a,d,_ in tights ]}")
report(f"  (expected: exactly [(1, 1)] = the set {{1,...,12}})")

report("")
report("=" * 74)
report("(D) general-n extension  q = 2a+(n-1)d  (even n-1 => q odd, clean)")
report("=" * 74)
report(f"{'n':>3} {'a':>3} {'d':>3} {'trueM':>9} {'1/2-(n-1)/2q':>13} "
       f"{'wit<=M':>7} {'=1/(n+1)@{1..n}?':>16}")
for n in [4, 6, 8, 10, 12]:
    for (a, d) in [(1, 1), (2, 1), (1, 3), (2, 3), (3, 5)]:
        if gcd(a, d) != 1:
            continue
        if (n - 1) % 2 != 1:     # need n-1 odd for q odd (n even). n even here.
            pass
        speeds = [a + d * k for k in range(n)]
        if max(speeds) > 160:
            continue
        M, _ = true_M(speeds)
        q = 2 * a + (n - 1) * d
        formula = F(1, 2) - F(n - 1, 2 * q)
        wit_le = formula <= M
        note = ""
        if a == 1 and d == 1:
            note = f"M={M} vs 1/(n+1)={F(1,n+1)} -> {M==F(1,n+1)}"
        report(f"{n:>3} {a:>3} {d:>3} {str(M):>9} {str(formula):>13} "
               f"{str(wit_le):>7}  {note}")

report("\nDONE.")
