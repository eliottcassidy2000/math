#!/usr/bin/env python3
"""
THE INVARIANT MEAN <Q_s> — closed form (boxeph-2026-07-16-S29)

THEOREM (the frame-averaged second moment). For any cluster E, section s,
P = 7 lcm E, endpoints (p_k, eps_k), and the signed coincidence counts
    N(h) := sum_{k,k'} eps_k eps_k' [ p_k == p_k' (mod h) ]   (h | P;
    N(1) = (sum eps)^2 = 0, N(P) = M for distinct endpoints):

 (K)  khat_P(n) = -csc^2(pi n/P)/(2P^2) EXACTLY for n != 0 (discrete second
      difference: Delta^2 K = -2/P^2 + (2/P) delta_0), so
 (Q)  Q_s(w) = (pi^2/P^2) sum_{n != 0} csc^2(pi n/P) |S(nw)|^2 EXACTLY.
 (C*) sum_{u in (Z/q)^*} csc^2(pi u/q) = J_2(q)/3 (Moebius + MI0),
      J_2 = Jordan totient.
 (P)  sum_{m in dZ_P} |S(m)|^2 = (P/d) N(P/d) (subgroup Parseval).
 (<>) averaging over the frame group (Z/P)^*:
      <Q_s> = (pi^2/P^2) sum_{g | P, g < P} (J_2(P/g)/3) (1/phi(P/g))
              sum_{g | d | P} mu(d/g) (P/d) N(P/d)
      -- an EXACT closed form: pi^2 times a rational functional of the
      endpoint coincidence spectrum {N(h)}_{h | P} on the divisor lattice.

Every step machine-refereed; final referee = direct average over all phi(P)
units vs the closed form.
"""

import sys
from math import gcd, lcm, pi, sin
import cmath

sys.path.insert(0, '04-computation')
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints, Qs_bilinear

def csc2(x):
    return 1.0 / sin(pi * x) ** 2

def divisors(n):
    ds = [d for d in range(1, n + 1) if n % d == 0]
    return ds

def mu(n):
    m, p, r = n, 2, 1
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    if m > 1:
        r = -r
    return r

def J2(n):
    r = n * n
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            r = r // (p * p) * (p * p - 1)
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        r = r // (m * m) * (m * m - 1)
    return r

def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            r -= r // p
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        r -= r // m
    return r

def cluster_data(E, s):
    P = 7 * lcm(*E)
    pts = endpoints(E, s)
    pos = [int(p * P) for p, sg, o in pts]
    sgn = [sg for p, sg, o in pts]
    return P, pos, sgn

def S_spec(pos, sgn, P):
    """full spectrum |S(n)|^2, n = 0..P-1 (O(PM))."""
    return [abs(sum(sg * cmath.exp(2j * pi * n * q / P)
                    for sg, q in zip(sgn, pos))) ** 2 for n in range(P)]

def Ncoin(pos, sgn, h):
    tot = 0
    from collections import Counter
    c = Counter()
    for q, sg in zip(pos, sgn):
        c[q % h] += sg
    return sum(v * v for v in c.values())

def run(E, s, name, full_avg=True):
    P, pos, sgn = cluster_data(E, s)
    M = len(pos)
    if M == 0:
        print(f"  [{name}] empty R_s; skip")
        return
    spec = S_spec(pos, sgn, P)

    # (K) exact kernel
    okK = True
    for n in (1, 2, 7, P // 2):
        kh = sum((j / P) * (1 - j / P) * cmath.exp(-2j * pi * n * j / P)
                 for j in range(P)) / P
        if abs(kh - (-csc2(n / P) / (2 * P * P))) > 1e-12:
            okK = False
    # (Q) exact spectral form vs bilinear
    okQ = True
    for w in (1, 11):
        qs_spec = (pi * pi / (P * P)) * sum(
            csc2(n / P) * spec[(n * w) % P] for n in range(1, P))
        qs_bil = Qs_bilinear(E, s, w)
        if abs(qs_spec - qs_bil) > 1e-6 * max(1, abs(qs_bil)):
            okQ = False
    # (C*) unit csc^2 sums
    okC = all(abs(sum(csc2(u / q) for u in range(1, q) if gcd(u, q) == 1)
                  - J2(q) / 3) < 1e-6 for q in range(2, 40))
    # (P) subgroup Parseval
    okP = True
    for d in divisors(P):
        if d == P or P // d > 2000:
            continue
        direct = sum(spec[(d * u) % P] for u in range(1, P // d))
        # sum_{m in dZ, m != 0} |S|^2 = (P/d) N(P/d)  (S(0) = 0)
        closed = (P // d) * Ncoin(pos, sgn, P // d)
        if abs(direct - closed) > 1e-6 * max(1, closed):
            okP = False
    print(f"  [{name}] P={P} M={M}: (K) {okK}; (Q) {okQ}; (C*) {okC}; (P) {okP}")

    # (<>) closed form
    total = 0.0
    for g in divisors(P):
        if g == P:
            continue
        q = P // g
        inner = 0
        for d in divisors(P):
            if d % g == 0:
                inner += mu(d // g) * (P // d) * Ncoin(pos, sgn, P // d)
        total += (J2(q) / 3) * inner / phi(q)
    closed_mean = (pi * pi / (P * P)) * total

    if full_avg:
        units = [w for w in range(1, P) if gcd(w, P) == 1]
        direct_mean = 0.0
        vals = []
        for w in units:
            qsw = (pi * pi / (P * P)) * sum(
                csc2(n / P) * spec[(n * w) % P] for n in range(1, P))
            vals.append(qsw)
        direct_mean = sum(vals) / len(vals)
        vmin, vmax = min(vals), max(vals)
        err = abs(direct_mean - closed_mean) / max(direct_mean, 1e-12)
        print(f"      <Q_s> closed = {closed_mean:.6f}; direct = {direct_mean:.6f} "
              f"(rel err {err:.2e}); fluctuation range [{vmin:.3f}, {vmax:.3f}]; "
              f"<Q_s>/M = {closed_mean/M:.4f}")
    else:
        print(f"      <Q_s> closed = {closed_mean:.6f}; <Q_s>/M = {closed_mean/M:.4f}")
    # the coincidence spectrum (the arithmetic content)
    Ns = {h: Ncoin(pos, sgn, h) for h in divisors(P) if h <= 84 or h == P}
    print(f"      coincidence spectrum N(h) (h | P, h <= 84): "
          f"{[(h, Ns[h]) for h in sorted(Ns) if h <= 84]}")

if __name__ == "__main__":
    print("=" * 74)
    print("THE INVARIANT MEAN <Q_s> -- closed form referees")
    run([1, 2, 3, 4, 5, 6, 60], 0, "family {1..6,60}")
    run([1, 2, 3, 4, 5, 36, 60], 0, "two-owner {1..5,36,60}")
    run([12, 15, 20, 21, 28, 30, 35], 0, "balanced {12,15,20,21,28,30,35}")
    run([1, 2, 3, 4, 5, 6, 120], 0, "family {1..6,120}")
    print("done")
