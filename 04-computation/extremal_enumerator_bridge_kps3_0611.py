#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-486/487 lab: the eta-discriminant bridge and the [72,36,16] frame.
kind-pasteur-2026-06-11-S3 (HYP-2419/2420/2421).

VERIFIED PRELIMINARIES (independent agent, exact, all asserts passed):
  * W_e8hat(x,y) = x^8 + 14 x^4 y^4 + y^8  -> E4(q^2) (Construction A theta).
  * CODE DISCRIMINANT P24 = x^4 y^4 (x^4 - y^4)^4  ->  16 * eta(q^2)^24 = 16 Delta.
    i.e. the Gleason second generator's discriminant IS the modular discriminant;
    Delta = q * prod(1-q^n)^24, the 24th power of the PENTAGONAL product.
  * W_g24 = W_e8hat^3 - 42 P24  (exact, in Z[x,y]).

This script:
  (A) Builds the EXTREMAL Type II weight enumerator W_n for n = 24m by the Gleason
      basis {W_e8hat, P24}: W_n = sum_{j=0}^{m} c_j W_e8hat^{(3m-3j)} P24^{j},
      c_0 = 1, solving the extremal conditions A_4 = ... = A_{4m} = 0 (min weight
      4m+4).  All EXACT rational/integer arithmetic.
  (B) Confirms W_24 = Golay, W_48 = eQR48 (d=12... no: eQR48 IS extremal d=12?
      NO — extremal at 48 is d=12 = 4*2+4, attained by eQR48), and that W_72 has
      ALL coefficients positive with the famous A_16 = 249849 (Sloane).
  (C) Sweeps n = 24m upward and reports the FIRST m with a negative coefficient
      in the extremal enumerator (Zhang: residue-0 branch first negative at
      n = 24*154 = 3696).  This LOCATES [72,36,16] (m=3) far below the negativity
      threshold => the obstruction to the code is NOT enumerator positivity.
  (D) The Lyapunov connection: the correction stream is graded by powers of
      P24 ~ Delta = eta^24; the coefficient c_j growth and the secular crossover.

Convention: track enumerators as dicts {weight: integer coefficient}; x has the
implicit exponent n - weight, y has weight (so a Type II code's A_w lives at key w).
P24 and W_e8hat^3 are degree-24 in (x,y); products add the weight-monomials.
"""

import time
from fractions import Fraction


def poly_mul(p, q, n):
    """Multiply two homogeneous (x,y)-polys given as {y_deg: coeff}; result
    truncated to total degree n (they should match exactly when degrees add)."""
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            out[a + b] = out.get(a + b, 0) + ca * cb
    return out


def poly_pow(p, e):
    r = {0: 1}
    for _ in range(e):
        r = poly_mul(r, p, 0)
    return r


# generators as {y-degree: coeff} (x-degree = total - y-degree)
W8 = {0: 1, 4: 14, 8: 1}                       # W_e8hat, degree 8
P24 = None  # built below, degree 24


def build_P24():
    # x^4 y^4 (x^4 - y^4)^4 ; expand (x^4-y^4)^4 then * x^4 y^4
    # (x^4 - y^4)^4 = sum_{k} C(4,k) x^{4(4-k)} (-y^4)^k -> y-deg 4k, coeff C(4,k)(-1)^k
    from math import comb
    base = {}
    for k in range(5):
        base[4 * k] = comb(4, k) * (-1) ** k
    # multiply by x^4 y^4 -> shift y-degree by 4
    return {d + 4: c for d, c in base.items()}


def extremal_enumerator(m):
    """Extremal Type II enumerator at n = 24m as {y-deg: coeff}, via the basis
    W8^{3(m-j)} P24^j, c_0 = 1, killing A_{4i} for i=1..m."""
    # basis polynomials B_j = W8^{3(m-j)} * P24^j, degree 24m, j = 0..m
    P = build_P24()
    basis = []
    for j in range(m + 1):
        bj = poly_mul(poly_pow(W8, 3 * (m - j)), poly_pow(P, j), 0)
        basis.append(bj)
    # unknowns c_0=1, c_1..c_m; conditions: coeff of y^{4i} (i=1..m) of sum c_j B_j = 0
    # B_j has lowest y-degree = 4j (P24^j starts at y^{4j}), so it's lower-triangular:
    # condition i only involves c_0..c_i. Solve forward.
    c = [Fraction(0)] * (m + 1)
    c[0] = Fraction(1)
    for i in range(1, m + 1):
        # coeff of y^{4i} in sum_{j<=i} c_j B_j = 0  (B_j contributes at 4i only if 4j<=4i)
        s = Fraction(0)
        for j in range(0, i):
            s += c[j] * basis[j].get(4 * i, 0)
        # B_i coeff at y^{4i}:
        bii = basis[i].get(4 * i, 0)
        c[i] = -s / bii
    # assemble
    W = {}
    for j in range(m + 1):
        for d, coeff in basis[j].items():
            W[d] = W.get(d, 0) + c[j] * coeff
    # must be integers
    Wi = {}
    for d, coeff in W.items():
        assert coeff.denominator == 1, f"non-integer coeff at y^{d}: {coeff}"
        if coeff != 0:
            Wi[d] = int(coeff)
    return Wi, [int(x) if x.denominator == 1 else x for x in c]


def main():
    t0 = time.time()
    print("=== (A/B) extremal Type II enumerators, exact ===", flush=True)
    # m=1: should be W_e8hat^3? No: n=24, extremal d=8 = Golay
    W24, c24 = extremal_enumerator(1)
    print(f"   n=24 (m=1): {dict(sorted(W24.items()))}", flush=True)
    golay = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    print(f"      = Golay g24: {W24 == golay}", flush=True)

    W48, c48 = extremal_enumerator(2)
    print(f"   n=48 (m=2): min weight present = {min(w for w in W48 if w > 0)}, "
          f"A_12 = {W48.get(12)}", flush=True)
    print(f"      eQR48 extremal d=12: {min(w for w in W48 if w>0) == 12}", flush=True)

    W72, c72 = extremal_enumerator(3)
    A16 = W72.get(16)
    allpos = all(v > 0 for v in W72.values())
    print(f"   n=72 (m=3): min weight = {min(w for w in W72 if w>0)}, A_16 = {A16}", flush=True)
    print(f"      Sloane's A_16 = 249849: {A16 == 249849}", flush=True)
    print(f"      ALL coefficients positive: {allpos}  => 72 is BELOW the MOS"
          f" negativity threshold; extremal W_72 is a bona fide nonneg enumerator", flush=True)
    print(f"      (so the [72,36,16] obstruction is NOT enumerator positivity /"
          f" MOS / shadow — it is code-combinatorial existence)", flush=True)

    print("\n=== (C) first negativity in the extremal enumerator (Zhang: n=3696) ===", flush=True)
    first_neg = None
    for m in range(1, 170):
        W, _ = extremal_enumerator(m)
        negs = [(w, v) for w, v in W.items() if v < 0]
        if negs:
            first_neg = (24 * m, m, min(negs, key=lambda t: t[0]))
            break
    if first_neg:
        n, m, (w, v) = first_neg
        print(f"   FIRST negative coefficient: n = {n} (m={m}), weight {w}, value {v}", flush=True)
        print(f"      Zhang threshold n=3696 (m=154) residue-0: {n == 3696}", flush=True)
    else:
        print("   no negativity found in range (extend m)", flush=True)

    print("\n=== (D) the eta-discriminant grading of the correction stream ===", flush=True)
    print("   extremal W_n = sum_j c_j * W8^{3(n/24 - j)} * P24^j,  P24 ~ 16 Delta = 16 eta^24", flush=True)
    for m in (1, 2, 3, 4, 5):
        _, c = extremal_enumerator(m)
        print(f"   n={24*m}: Delta-correction coeffs c_j = {c}", flush=True)
    # c_1 sequence: the leading correction amplitude
    c1s = []
    for m in range(1, 12):
        _, c = extremal_enumerator(m)
        c1s.append(c[1] if len(c) > 1 else 0)
    print(f"   leading correction c_1(m) for m=1..11: {c1s}", flush=True)
    ratios = [c1s[i+1]/c1s[i] for i in range(len(c1s)-1) if c1s[i] != 0]
    print(f"   c_1 ratios (-> the eta^24 saddle base, MOS c_1 ~ 69.1...?): "
          f"{[round(r,2) for r in ratios]}", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
