#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_apex_d0_proofaudit_macmini.py  (mac-mini 2026-06-21, THREAD D audit)

ADVERSARIAL audit of the apex-prime D=0 law (HYP-2733) step (D).

The claimed mechanism: row0 flat <=> band-start phases {7pk/q mod 7} tile Z/7 evenly
<=> 7|p. I want to verify that the TRUE reason row0 is flat is a clean statement and
that the law is robust well beyond the [1,40) box, AND give a CLEAN INDEPENDENT proof
that mu_{p,q} factorizes iff 7|pq via the Weyl/character criterion, which is rigorous
for ALL coprime (p,q), not just a finite box.

CHARACTER PROOF (the rigorous one I am checking):
  mu factorizes to uniform iff all nonzero Fourier coefficients of (sector(qv),sector(pv))
  joint indicator vanish. The 2D sector occupancy mu(i,j) has Fourier transform over
  Z/7 x Z/7:
    hatmu(a,b) = int_0^1 e(-a sector(qv)/7 - b sector(pv)/7) dv   [discrete sector chars]
  Cleaner: use the *continuous* phases. mu is uniform iff for all (a,b) != (0,0) in (Z/7)^2,
    sum over the joint sector law of chi_a(i) chi_b(j) = 0.
  Equivalent characterization (what I test): mu(i,j) = 1/49 for all i,j.

  Instead of Fourier I give the DIRECT measure-theoretic reason and test it exhaustively
  on a MUCH larger coprime box + random large (p,q) to stress the "all coprime" claim.
"""
from fractions import Fraction as Fr
from math import gcd
import random

P = 7
def sector(yf): return int(P * yf)

def mu_full(p, q):
    bp = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(0, P * f):
            bp.add(Fr(t, P * f))
    vs = sorted(bp); cell = {}
    for a, b in zip(vs, vs[1:]):
        mid = (a + b) / 2
        key = (sector((q*mid) % 1), sector((p*mid) % 1))
        cell[key] = cell.get(key, Fr(0)) + (b - a)
    return cell

def D_pq(p, q):
    cell = mu_full(p, q)
    tot = Fr(0)
    for i in range(P):
        for j in range(P):
            tot += abs(cell.get((i, j), Fr(0)) - Fr(1, 49))
    return tot

def factorizes(p, q):
    """Is mu(i,j) = mu_i(row marg) * mu_j(col marg)?  marginals are always 1/7."""
    cell = mu_full(p, q)
    return all(cell.get((i, j), Fr(0)) == Fr(1, 49) for i in range(P) for j in range(P))

def main():
    print("#" * 80)
    print("# APEX D=0 LAW AUDIT (HYP-2733)")
    print("#" * 80)

    # 1) Extend the verification box well beyond [1,40): test all coprime (p,q), p,q<=60
    print("\n=== 1) law D=0 <=> 7|pq, EXTENDED box p,q in [1,60] ===")
    viol = 0; checked = 0
    worst_examples = []
    for q in range(1, 61):
        for p in range(1, 61):
            if gcd(p, q) != 1: continue
            checked += 1
            d = D_pq(p, q)
            law = (p % 7 == 0) or (q % 7 == 0)
            iszero = (d == 0)
            if iszero != law:
                viol += 1
                if len(worst_examples) < 10:
                    worst_examples.append((p, q, d, law))
    print(f"  checked {checked} coprime pairs; violations of (D=0 <=> 7|pq): {viol}")
    if worst_examples:
        print("  EXAMPLES:", worst_examples)

    # 2) Random LARGE coprime pairs (stress "all coprime", since proof claims universality)
    print("\n=== 2) random LARGE coprime pairs (p,q up to ~500) ===")
    random.seed(108)
    viol2 = 0
    for _ in range(120):
        p = random.randint(1, 500); q = random.randint(1, 500)
        if gcd(p, q) != 1: continue
        d = D_pq(p, q)
        law = (p % 7 == 0) or (q % 7 == 0)
        if (d == 0) != law:
            viol2 += 1
            print(f"  VIOLATION p={p} q={q} D={d} law={law}")
    print(f"  random large pairs: violations = {viol2}")

    # 3) The mechanism: is D=0 EXACTLY equivalent to mu factorizing? (claim in HYP-2733)
    print("\n=== 3) D=0 <=> mu factorizes to uniform 1/49 (equivalence claim) ===")
    bad = 0
    for q in range(1, 41):
        for p in range(1, 41):
            if gcd(p, q) != 1: continue
            if (D_pq(p, q) == 0) != factorizes(p, q):
                bad += 1
    print(f"  D=0 <=> factorizes mismatches: {bad}")

    # 4) ADVERSARIAL: does 7|p alone (not 7|q) suffice, and is q irrelevant to D=0?
    #    Test the asymmetric claim: 7|p, any q coprime -> D=0 ; 7nmid p,q -> D>0.
    print("\n=== 4) asymmetric form: 7|p (7nmid q) => D=0 across q ===")
    cnt = 0
    for p in [7, 14, 21, 28]:
        for q in range(1, 30):
            if gcd(p, q) != 1: continue
            if q % 7 == 0: continue
            if D_pq(p, q) != 0:
                print(f"  FAIL 7|p but D>0: p={p} q={q} D={D_pq(p,q)}")
            else:
                cnt += 1
    print(f"  confirmed 7|p => D=0 cases: {cnt}")

    # 5) sup of D*p and D*q on the non-apex locus (the discrepancy bound check 14/p, 20/7)
    print("\n=== 5) sup D*p and D*q over 7nmid pq, p,q<=50 (ledger: sup D*p=20/7) ===")
    sup_Dp = Fr(0); arg_Dp = None
    sup_Dq = Fr(0); arg_Dq = None
    for q in range(1, 51):
        if q % 7 == 0: continue
        for p in range(1, 51):
            if p % 7 == 0: continue
            if gcd(p, q) != 1: continue
            d = D_pq(p, q)
            if d * p > sup_Dp: sup_Dp = d * p; arg_Dp = (p, q)
            if d * q > sup_Dq: sup_Dq = d * q; arg_Dq = (p, q)
    print(f"  sup D*p = {sup_Dp} = {float(sup_Dp):.5f} at {arg_Dp}  (ledger 20/7={float(Fr(20,7)):.5f})")
    print(f"  sup D*q = {sup_Dq} = {float(sup_Dq):.5f} at {arg_Dq}  (ledger 12/7={float(Fr(12,7)):.5f})")
    # check D <= 14/p (the PROVED Koksma bound) on the whole box
    viol_koksma = 0
    for q in range(1, 51):
        for p in range(1, 51):
            if gcd(p, q) != 1: continue
            if D_pq(p, q) > Fr(14, p):
                viol_koksma += 1
    print(f"  D <= 14/p violations over p,q<=50: {viol_koksma}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
