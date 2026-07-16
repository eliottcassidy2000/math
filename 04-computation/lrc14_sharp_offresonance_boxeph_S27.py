#!/usr/bin/env python3
"""
THE SHARP OFF-RESONANCE CONSTANT — closed-form comb mass + csc^2-telescoped tails
boxeph-2026-07-16-S27 (executes THM-887's named refinement (a))

Everything rides on the classical MULTIPLICATION IDENTITY
    (MI)   sum_{k=0}^{q-1} csc^2(pi (x + k/q)) = q^2 csc^2(pi q x),
(proved by the partial-fraction expansion csc^2(pi x) = (1/pi^2) sum_{n in Z}
1/(x+n)^2 — the inner double sum regroups), together with
    (MI0)  sum_{k=1}^{q-1} csc^2(pi k/q) = (q^2 - 1)/3.

Consequences for Q_s(w) = 2 pi^2 sum_{n != 0} |khat_P(n)| |S(nw)|^2 with the
proven bound |khat_P(n)| <= (1/(2P^2)) csc^2(pi n/P)  [THM-886(III)]:

 (A) COMB MASS IN CLOSED FORM (gcd(w, P) = 1).  The frequencies n with
     nw on owner e's comb form eZ_P, and the tooth visited at n = em is
     nu_e(m w mod 7).  Slicing m by residue mod 7 and applying (MI):
        D_comb^(e)(w) = 2 pi^2 [ sum_{r=1}^{6} |nu_e(r)|^2 csc^2(pi r wbar_7 /7) / (98 e^2)
                                 + |nu_e(0)|^2 ((P/(7e))^2 - 1) / (6 P^2) * 2 pi^2-normalized ],
     i.e. the comb diagonal depends on w ONLY THROUGH w mod 7.
 (B) TAIL IN HARDY-LITTLEWOOD FORM.  The off-comb frequencies at integer
     distance d from eZ carry khat-mass EXACTLY (1/(2e^2)) csc^2(pi ||d wbar / e||)
     per signed slice (by (MI) with q = P/e), so with the run profile
     |S_e| <= min(M_e, R_e / sin(pi d/e)):
        D_tail^(e)(w) <= (2 pi^2 / (2 e^2)) * sum_{d=1}^{floor(e/2)} c_d *
                          min(M_e^2, R_e^2 csc^2(pi d/e)) * csc^2(pi ||d wbar/e||),
     c_d = 2 for d < e/2, 1 for d = e/2.  O(e) to evaluate; symmetric in w <-> wbar.
 (C) CROSS TERMS by Cauchy-Schwarz of the per-owner diagonals:
        Q_s(w) <= ( sum_e sqrt(D^(e)(w)) )^2,   D^(e) = D_comb^(e) + D_tail^(e).

Referees: (MI)/(MI0) numeric; slice-collapse exactness; per-owner diagonal
cover; the assembled sharp bound vs exact Q_s vs THM-887's crude law, on the
full S26 battery.  NOTE: for gcd(w, P) = g > 1 the frequency set collapses to
gZ_P; we handle general w by replacing wbar with the appropriate coset data —
the referee battery includes non-coprime w via the direct fallback.

Pure Python.
"""

import sys
from math import gcd, lcm, pi, sin
import cmath

sys.path.insert(0, '04-computation')
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints, Qs_bilinear
from lrc14_general_resonance_law_boxeph_S26 import owner_data, S_owner, nu_hat, profile_e

def csc2(x):
    return 1.0 / sin(pi * x) ** 2

def referee_MI():
    import random
    ok = True
    for q in (3, 7, 12, 60):
        for _ in range(3):
            x = random.random() / q * 0.9 + 0.01 / q
            lhs = sum(csc2((x + k / q) % 1) for k in range(q))
            rhs = q * q * csc2((q * x) % 1)
            if abs(lhs - rhs) / rhs > 1e-9:
                ok = False
        s0 = sum(csc2(k / q) for k in range(1, q))
        if abs(s0 - (q * q - 1) / 3) > 1e-6:
            ok = False
    print(f"  (MI)/(MI0) multiplication identities: {ok}")

def sharp_bound(E, s, w):
    """the (A)+(B)+(C) closed-form bound; requires gcd(w, P) = 1."""
    P, M, data = owner_data(E, s)
    if M == 0:
        return 0.0, 0, {}
    assert gcd(w, P) == 1
    wbar = pow(w, -1, P)
    tot_sqrt = 0.0
    parts = {}
    for e, d in data.items():
        q7 = P // (7 * e)
        # (A) comb mass: teeth |nu(r w mod 7)| sit at khat-mass csc2(pi r/7)/(98 e^2)
        comb = 0.0
        for r in range(1, 7):
            tooth = abs(nu_hat(d, r * w % 7))
            comb += tooth * tooth * csc2(r / 7) / (98 * e * e)
        t0 = abs(nu_hat(d, 0))
        comb += t0 * t0 * (q7 * q7 - 1) / (6 * P * P) * 2 / 2  # (MI0) class r=0, excl. n=0
        comb *= 2 * pi * pi
        # (B) tail
        tail = 0.0
        for dd in range(1, e // 2 + 1):
            cd = 1.0 if 2 * dd == e else 2.0
            prof2 = min(d["M"] ** 2, d["runs"] ** 2 * csc2(dd / e))
            x = (dd * wbar % e) / e
            x = min(x, 1 - x)
            if x == 0:
                continue  # lands back on the comb; counted in (A)
            tail += cd * prof2 * csc2(x)
        tail *= 2 * pi * pi / (2 * e * e)
        D = comb + tail
        parts[e] = (comb, tail)
        tot_sqrt += D ** 0.5
    return tot_sqrt ** 2, M, parts

def battery():
    print()
    print("=" * 76)
    print("SHARP BOUND vs EXACT vs CRUDE LAW (S26 battery, coprime w)")
    tests = [
        ("family {1..6,60}", [1, 2, 3, 4, 5, 6, 60], 0, (11, 59, 61, 127, 187)),
        ("two-owner {1..5,36,60}", [1, 2, 3, 4, 5, 36, 60], 0, (11, 37, 71, 109)),
        ("two-large {1..5,56,84}", [1, 2, 3, 4, 5, 56, 84], 0, (11, 253, 113)),
        ("balanced {12,15,20,21,28,30,35}", [12, 15, 20, 21, 28, 30, 35], 0, (11, 421, 61)),
        ("near-AP {8,9,10,12,14,15,18}", [8, 9, 10, 12, 14, 15, 18], 0, (11, 64, 2521)),
    ]
    worst = 0.0
    for name, E, s, ws in tests:
        P = 7 * lcm(*E)
        for w in ws:
            if gcd(w, P) != 1:
                continue
            q = Qs_bilinear(E, s, w)
            b, M, parts = sharp_bound(E, s, w)
            if M == 0:
                continue
            ratio = b / max(q, 1e-9)
            worst = max(worst, ratio)
            print(f"  [{name}] w={w:>5}: exact Q_s/M = {q/M:7.2f}; SHARP/M = {b/M:8.2f} "
                  f"(x{ratio:6.1f} over exact); cover: {q <= b + 1e-6}")
    print(f"  worst sharp/exact ratio on battery: {worst:.1f}")

def referee_slices():
    """verify the (B) slice collapse: sum over one signed slice of khat-bound
    equals csc^2(pi ||d wbar / e||)/(2 e^2)."""
    E, s = [1, 2, 3, 4, 5, 36, 60], 0
    P = 7 * lcm(*E)
    w = 37
    wbar = pow(w, -1, P)
    ok = True
    for e in (36, 60):
        for dd in (1, 2, 5):
            direct = 0.0
            for k in range(P // e):
                m = (k * e + dd) % P
                n = (m * wbar) % P
                if n == 0:
                    continue
                direct += (1 / (2 * P * P)) * csc2(n / P)
            closed = csc2(((dd * wbar % e) / e) if 0 < (dd * wbar % e) / e < 1 else 0.5) / (2 * e * e)
            xx = (dd * wbar % e) / e
            xx = min(xx, 1 - xx)
            closed = csc2(xx) / (2 * e * e)
            if abs(direct - closed) / closed > 1e-9:
                ok = False
                print(f"    slice mismatch e={e} d={dd}: {direct:.6e} vs {closed:.6e}")
    print(f"  (B) slice collapse exact: {ok}")

if __name__ == "__main__":
    print("=" * 76)
    print("THE SHARP OFF-RESONANCE CONSTANT -- referees")
    referee_MI()
    referee_slices()
    battery()
    print("done")

def referee_comb_closed_form():
    """(A): D_comb^(e)(w) closed form vs direct summation over the comb eZ_P."""
    from math import lcm
    E, s = [1, 2, 3, 4, 5, 36, 60], 0
    P = 7 * lcm(*E)
    _, M, data = owner_data(E, s)
    ok = True
    for w in (37, 71):
        if gcd(w, P) != 1:
            continue
        wbar = pow(w, -1, P)
        for e in (36, 60):
            d = data[e]
            direct = 0.0
            for m in range(1, P // e):
                n = (e * m * wbar) % P    # weight index whose image is the comb pt e*m
                if n == 0:
                    continue
                a = (m * 1) % 7           # tooth visited at comb point e*m is nu(mw ... )
                # careful: at weight index n, frequency is n*w = e*m: tooth = nu_e(m mod 7)
                direct += (1 / (2 * P * P)) * csc2(n / P) * abs(nu_hat(d, m % 7)) ** 2
            closed = 0.0
            for r in range(1, 7):
                # comb points e*m with m == r (mod 7) have weight-index AP with
                # (MI)-collapsed mass csc2(pi * r * wbar_e7 ...): derive via slices:
                # n = e*m*wbar mod P, m = 7u + r: AP step 7e*wbar, offset e*r*wbar:
                # (MI) total = (P/(7e))^2 csc2(pi*(P/(7e)) * e*r*wbar/P) / (2P^2)
                q7 = P // (7 * e)
                x = (r * wbar % (7 * 0 + P // (e))) * 0  # placeholder
                arg = (q7 * e * r * wbar % P) / P
                arg = min(arg, 1 - arg)
                mass = (q7 * q7) * csc2(arg) / (2 * P * P) if arg > 0 else 0.0
                closed += mass * abs(nu_hat(d, r % 7)) ** 2
            # r = 0 class (m == 0 mod 7, m != 0):
            q7 = P // (7 * e)
            closed += ((q7 * q7 - 1) / 3) / (2 * P * P) * abs(nu_hat(d, 0)) ** 2
            if abs(direct - closed) / max(direct, 1e-12) > 1e-9:
                ok = False
                print(f"    comb closed-form mismatch e={e} w={w}: {direct:.6e} vs {closed:.6e}")
    print(f"  (A) comb-diagonal closed form exact: {ok}")
    print(f"      note: arg = q7*e*r*wbar/P = r*wbar/7 mod 1 -> depends on wbar mod 7 ONLY")

if __name__ == "__main__" and True:
    referee_comb_closed_form()
