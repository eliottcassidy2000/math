#!/usr/bin/env python3
"""
THM-438 ADDENDUM-8 exploration (monad-explorer-2026-06-07, 10th session).

Builds on ADDENDUM-7 (column denominator (1-x)^{2m-1} PROVED).

THE BINOMIAL REFRAMING.
-----------------------
The cycle-rank column GF is  T_m(x) = Sum_e R_s(m,e) (x/(1-x))^e  with e=m..2m-1.
Since [x^k] (x/(1-x))^e = C(k-1, e-1), this is EXACTLY a binomial expansion of the count:

        t(k,m) = Sum_{e=m}^{2m-1} R_s(m,e) * C(k-1, e-1).            (B)

Immediate consequences (this script verifies all of them):
  * t(k,m) = 0 automatically for k = 1, ..., m-1  (binomial support: C(k-1,e-1)=0
    when 0 <= k-1 <= e-2, i.e. 1 <= k <= e-1, and every e >= m).
  * t(0,m) [polynomial extension] = Sum_e R_s(m,e) C(-1,e-1)
                                  = Sum_e R_s(m,e) (-1)^{e-1} = -Q_m(-1).
    So  deg P_m = m-2  <=>  Q_m(-1)=0  <=>  t(0,m)=0  <=>  T_m(x) -> 0 as x->oo.
    (handoff #1, restated as a SINGLE transparent statement.)
  * t(-1,m) [poly] = Sum_e R_s(m,e) C(-2,e-1) = Sum_e R_s(m,e)(-1)^{e-1} e = Q_m'(-1).
    Given Q_m(-1)=0 this equals lead(P_m).  Conjecture (handoff #2): = 2^m - 1.
  * Reading P_m TOP-DOWN is the Taylor expansion of V(t,y)=Sum_m Q_m(t) y^m at t=-1:
        d^j/dt^j V(t,y)|_{t=-1} / j!  <->  coeff of x^{m-2-j} in P_m, summed over m.

We verify (B) against the known triangle, recompute deg/lead, tabulate the dual
sequence t(-k,m), the top-down P_m coefficients, and emit sequences for OEIS lookup.
"""
from fractions import Fraction
from math import comb

# ---- Known data (THM-438 ADD-6/7) ------------------------------------------
# R_s(m,e), e=m..2m-1  (the s-expansion coeffs; ADD-7 honesty: NOT reduced counts)
Rs = {
    1: {1: 1},
    2: {2: 3, 3: 3},
    3: {3: 13, 4: 33, 5: 20},
    4: {4: 69, 5: 304, 6: 416, 7: 181},
}
# Known cycle-rank triangle t(k,m), k=1..6  (ADD-3 rows + ADD-6 k=6 row)
TRI = {
    1: {1: 1},
    2: {1: 1, 2: 3},
    3: {1: 1, 2: 9, 3: 13},
    4: {1: 1, 2: 18, 3: 72, 4: 69},
    5: {1: 1, 2: 30, 3: 230, 4: 580, 5: 421},
    6: {1: 1, 2: 45, 3: 560, 4: 2626, 5: 4845, 6: 2867},
}
# Known numerator polynomials P_m (coeff list low->high degree)
P = {1: [1], 2: [3], 3: [13, 7], 4: [69, 97, 15]}
A088368 = {1: 1, 2: 3, 3: 13, 4: 69, 5: 421, 6: 2867}


def gen_binom(k, e):
    """C(k-1, e-1) as a POLYNOMIAL in k, valid for all integer k (incl. negative)."""
    # C(n, r) = n(n-1)...(n-r+1)/r!  (polynomial extension)
    n = k - 1
    r = e - 1
    if r < 0:
        return Fraction(0)
    num = Fraction(1)
    for i in range(r):
        num *= (n - i)
    return num / Fraction(1) * Fraction(1, 1) / _fact(r)


def _fact(r):
    f = 1
    for i in range(2, r + 1):
        f *= i
    return f


def t_poly(k, m):
    """t(k,m) via the binomial reframing (B), polynomial extension to any integer k."""
    return sum(Rs[m][e] * gen_binom(k, e) for e in Rs[m])


print("=" * 78)
print("PART 1 — verify the binomial reframing  t(k,m) = Sum_e R_s(m,e) C(k-1,e-1)")
print("=" * 78)
ok = True
for k in range(1, 7):
    for m in range(1, k + 1):
        if m not in Rs:
            continue
        pred = t_poly(k, m)
        true = TRI[k][m]
        flag = "" if pred == true else "  <-- MISMATCH"
        if pred != true:
            ok = False
        print(f"  t({k},{m}): binom={pred}  triangle={true}{flag}")
print(f"\n  ALL MATCH: {ok}")

print()
print("=" * 78)
print("PART 2 — forced zeros at k=1..m-1, and the single nontrivial value t(0,m)")
print("=" * 78)
for m in Rs:
    zeros = [t_poly(k, m) for k in range(1, m)]
    t0 = t_poly(0, m)
    Qm_at_m1 = sum(Rs[m][e] * (-1) ** e for e in Rs[m])
    print(f"  m={m}: t(1..{m-1},m)={[str(z) for z in zeros]}  t(0,m)={t0}  "
          f"Q_m(-1)={Qm_at_m1}  (deg P_m=m-2  <=> both 0)")

print()
print("=" * 78)
print("PART 3 — deg P_m and lead P_m from the binomial form")
print("=" * 78)
for m in Rs:
    # lead P_m = t(-1, m) (given Q_m(-1)=0)
    lead = t_poly(-1, m)
    twopow = 2 ** m - 1
    # also recompute deg via the explicit P
    degP = len(P[m]) - 1 if m in P else None
    print(f"  m={m}: t(-1,m)=lead P_m = {lead}   2^m-1 = {twopow}   "
          f"match={lead == twopow}   (deg P_m={degP}, m-2={m-2})")

print()
print("=" * 78)
print("PART 4 — the DUAL sequence t(-j, m) for j>=1  (reciprocity candidates)")
print("=" * 78)
for m in Rs:
    dual = [int(t_poly(-j, m)) for j in range(0, 9)]
    print(f"  m={m}: t(-j,m), j=0..8 = {dual}")
print("  (t(0,m)=0 expected for m>=2; t(-1,m)=2^m-1)")

print()
print("=" * 78)
print("PART 5 — P_m read TOP-DOWN = Taylor of V(t,y) at t=-1 (per-column coeffs)")
print("=" * 78)
print("  P_m coefficients high->low degree (top = lead = 2^m-1, bottom = A088368(m)):")
for m in P:
    hi2lo = list(reversed(P[m]))
    print(f"   m={m}: {hi2lo}   [lead={hi2lo[0]}=2^m-1?, const={P[m][0]}=A088368({m})={A088368[m]}]")

print()
print("  Top-down coefficient columns across m (the Taylor-at-(-1) sequences):")
maxlen = max(len(P[m]) for m in P)
for j in range(maxlen):       # j = distance from top
    seq = []
    for m in sorted(P):
        if len(P[m]) > j:
            seq.append((m, P[m][len(P[m]) - 1 - j]))
    label = {0: "lead (2^m-1)", 1: "2nd-from-top", 2: "3rd-from-top"}.get(j, f"{j+1}-from-top")
    print(f"   distance {j} [{label}]: {seq}")

print()
print("=" * 78)
print("PART 6 — sequences for OEIS lookup")
print("=" * 78)
print(f"  top pole residue R_s(m,2m-1)=P_m(1):  {[sum(P[m]) for m in sorted(P)]}")
print(f"  diagonal A088368(m):                  {[A088368[m] for m in sorted(A088368)]}")
print(f"  lead P_m = 2^m-1:                     {[2**m-1 for m in range(1,7)]}")
print(f"  GF check Sum(2^m-1)y^m = y/((1-y)(1-2y))")
# verify GF y/((1-y)(1-2y)) coefficients
from itertools import islice
def gf_coeffs(N):
    # 1/((1-y)(1-2y)) = Sum (2^{n+1}-1) y^n ; times y shifts
    return [ (2**(n)-1) for n in range(1,N+1)]
print(f"     y/((1-y)(1-2y)) coeffs m=1..6: {gf_coeffs(6)}")
print(f"  dual t(-j,3): {[int(t_poly(-j,3)) for j in range(1,9)]}")
print(f"  dual t(-j,4): {[int(t_poly(-j,4)) for j in range(1,9)]}")
