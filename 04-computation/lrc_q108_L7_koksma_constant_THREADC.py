#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_koksma_constant_THREADC.py  (THREAD C audit, 2026-06-21)

The single most load-bearing inequality in the elementary D<=14/p proof is:

   | (1/q) sum_{m=0}^{q-1} h_j(a_m)  -  int_0^1 h_j(a) da |  <=  Var(h_j) * (1/q),

applied to the q EQUALLY-SPACED-with-shift points a_m = {c + m/q} on the circle,
where h_j is the trapezoid overlap function (a function on the CIRCLE [0,1)).

This is the Koksma / Koksma-Hlawka inequality on the torus:
   |error| <= V(h_j) * D*_N(points),
and for N equally-spaced points {c+m/N} the star-discrepancy is D* = 1/N
(the standard fact: equally spaced points have discrepancy 1/N; on the circle
the *extreme* discrepancy is 1/N and the relevant 'one-sided' constant works out
to exactly 1/N here because h_j is periodic and the sum is a closed Riemann sum).

ADVERSARIAL CHECK: compute the EXACT LHS error and the claimed RHS Var*1/q for
EVERY (p,q,i,j), and report the worst ratio LHS/RHS. If any ratio > 1 the
proof's Koksma step is WRONG. We also test the sharper claim that for a periodic
function the equally-spaced Riemann sum error is bounded by V * (1/q) with the
shift c making no difference to the bound.

We also directly stress the CIRCULAR vs INTERVAL discrepancy distinction:
h_j is periodic (a function of a on the circle), and the standard Koksma-Hlawka
for periodic f uses the *extreme* discrepancy D_N which for {c+m/N} equals 1/N.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def frac(x):
    return x - (x.numerator // x.denominator)

def overlap_arc_bin(a, L, j):
    a = frac(a)
    end = a + L
    segs = []
    if end <= 1:
        segs.append((a, end))
    else:
        segs.append((a, Fr(1)))
        segs.append((Fr(0), end - 1))
    lo = Fr(j, P); hi = Fr(j + 1, P)
    tot = Fr(0)
    for (s, e) in segs:
        tot += max(Fr(0), min(e, hi) - max(s, lo))
    return tot

def var_hj(L, j):
    crit = {Fr(0), Fr(1)}
    for c in (Fr(j, P), Fr(j + 1, P)):
        crit.add(frac(c)); crit.add(frac(c - L))
    allpts = sorted(crit)
    vals = [overlap_arc_bin(a, L, j) for a in allpts]
    tv = Fr(0)
    for k in range(len(vals) - 1):
        tv += abs(vals[k + 1] - vals[k])
    return tv

def int_hj(L, j):
    return L / 7  # verified exactly in micro-audit

def main():
    print("=" * 78)
    print("Koksma constant audit: |Riemann-sum-error| <= Var(h_j)*(1/q) ?")
    print("=" * 78)
    worst = (0, 0, 0, 0, 0.0)
    viol = 0
    checked = 0
    for q in range(1, 40):
        pmax = int(Fr(43, 20) * q)
        for p in range(q + 1, pmax + 1):
            if gcd(p, q) != 1:
                continue
            r = Fr(p, q)
            if not (Fr(1) < r <= Fr(43, 20)):
                continue
            L = Fr(p, P * q)
            checked += 1
            for i in range(P):
                for j in range(P):
                    # the q equally-spaced (shifted) points a_m = {p(i+7m)/(7q)}
                    a = [frac(Fr(p * (i + P * m), P * q)) for m in range(q)]
                    avg = sum(overlap_arc_bin(am, L, j) for am in a) / q
                    err = abs(avg - int_hj(L, j))
                    rhs = var_hj(L, j) * Fr(1, q)
                    if err > rhs:
                        viol += 1
                        if viol <= 10:
                            print(f"  KOKSMA VIOLATION p/q={p}/{q} i={i} j={j}: |err|={err}({float(err):.6f}) > Var/q={rhs}({float(rhs):.6f})")
                    ratio = float(err / rhs) if rhs != 0 else 0.0
                    if ratio > worst[4]:
                        worst = (p, q, i, j, ratio)
    print(f"\n  checked {checked} ratios; 49 cells each")
    print(f"  Koksma |err| <= Var*(1/q) VIOLATIONS: {viol}  (0 == the Koksma step is valid as stated)")
    print(f"  worst LHS/RHS ratio = {worst[4]:.4f} at p/q={worst[0]}/{worst[1]} cell=({worst[2]},{worst[3]})")
    print("  (ratio<=1 => the proof's Koksma bound holds with constant 1/q; ratio<1 => slack)")

    # Also: confirm the *general* fact used -- for ANY shift c and ANY function of
    # bounded variation V on the circle, the N-point equally spaced average error
    # is <= V/N. Test with the actual h_j and 200 random shifts.
    print("\n  Independent stress of 'equally-spaced Riemann sum error <= V/N' with random shifts:")
    import random
    random.seed(1)
    badshift = 0
    worst_shift = 0.0
    for _ in range(3000):
        q = random.randint(1, 30)
        p = random.randint(q + 1, max(q + 1, int(2.15 * q)))
        if gcd(p, q) != 1:
            continue
        j = random.randint(0, 6)
        L = Fr(p, P * q)
        c = Fr(random.randint(0, 9999), 10000)  # arbitrary shift, not the structured one
        a = [frac(c + Fr(m, q)) for m in range(q)]
        avg = sum(overlap_arc_bin(am, L, j) for am in a) / q
        err = abs(avg - int_hj(L, j))
        rhs = var_hj(L, j) * Fr(1, q)
        if err > rhs:
            badshift += 1
        if rhs != 0:
            worst_shift = max(worst_shift, float(err / rhs))
    print(f"    random-shift equally-spaced violations: {badshift} (0 == general fact holds)")
    print(f"    worst random-shift ratio = {worst_shift:.4f}")

if __name__ == "__main__":
    main()
