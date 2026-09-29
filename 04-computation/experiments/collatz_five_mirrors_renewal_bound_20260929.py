#!/usr/bin/env python3
"""Lemma R' (renewal bound) checked term by term, the weighted negative-power norm, and the bound at h <= 300 (S22).

Lemma R'.  For every h >= 1 and every exponent s,
    |mu_hat_h(2^s)| <= sum_(J=0)^h mass_h(J) * Ntilde_(J-1),
where mass_h(J) = P(T_(J+1) <= s < T_J) is the probability that exactly the bottom J levels have negative exponent
(T_j = a_j + ... + a_h, T_(h+1) = 0), Ntilde_n = sum_(m>=1) 2^-m |mu_hat_n(2^-m mod 3^n)| (Ntilde_(-1) = Ntilde_0 = 1).
Chernoff: mass_h(J) <= P(T_(J+1) <= s) <= exp(-(h-J) I) exp(-theta (s - (h-J) log2 3)) with theta = ln(2 log_3 2) < 0,
I = I(log2 3) = -ln 3^(h*-1); at s = h log2 3 - delta this is exp(-h I) e^(theta delta) (log2 3 - 1)^-J.
Checks: (a) Ntilde_n to n = 80 (m <= 60 exact, tail <= 2^-60) and its rate; (b) the per-J inequality |c_J| <= mass_h(J)
Ntilde_(J-1) at h = 40, 80 (c_J = contribution of the words with exactly J negative-exponent levels, from the exact DP);
(c) mass_h(J) against the Chernoff bound; (d) the bound sum_J mass_h(J) Ntilde_(J-1) at h = 40..300 (Ntilde beyond 80
extrapolated as 3.6 * 3^(-n/2), the largest observed constant) against the measured maxima M(h) of the level-300 run.
Run: python 04-computation/experiments/collatz_five_mirrors_renewal_bound_20260929.py
"""
from __future__ import annotations

import math
import re
import sys

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402
from collatz_five_mirrors_renewal_20260929 import closed_family_all, decompose  # noqa: E402

Q = 1 - 1 / LOG23                      # tilt parameter: geometric weights (e^theta / 2)^a with e^theta = 2 Q
THETA = math.log(2 * Q)                # < 0
IRATE = THETA * LOG23 - math.log(Q / (1 - Q))   # I(log2 3) = -ln 3^(h* - 1)


def mass_law(h: int, s: int, JMAX: int):
    """Exact mass_h(J) = P(T_(J+1) <= s < T_J) by negative-binomial convolution (cost distribution of N geometrics)."""
    smax = s + 1
    # p[N][t] = P(S_N = t) for t <= smax
    p = np.zeros(smax + 1)
    p[0] = 1.0
    dist = {0: p.copy()}
    for N in range(1, h + 1):
        new = np.zeros(smax + 1)
        for a in range(1, smax + 1):
            new[a:] += 2.0 ** (-a) * p[: smax + 1 - a]
        p = new
        dist[N] = p.copy()
    out = []
    for J in range(0, JMAX + 1):
        N = h - J
        if N < 0:
            out.append(0.0)
            continue
        pN = dist[N]
        if J == 0:
            out.append(float(pN[: s + 1].sum()))
        else:
            ts = np.arange(0, s + 1)
            out.append(float((pN[: s + 1] * 2.0 ** (-(s - ts))).sum()))
    return out


if __name__ == "__main__":
    print(f"constants: q = {Q:.6f}, e^theta = {2 * Q:.6f}, theta = {THETA:.6f}, I = {IRATE:.6f}, e^-I = {math.exp(-IRATE):.6f}, 1/(log2 3 - 1) = {1 / (LOG23 - 1):.6f}, product with 3^-1/2: {3 ** -0.5 / (LOG23 - 1):.6f}")
    print("== (a) weighted negative-power norm Ntilde_n = sum_(m>=1) 2^-m |mu_hat_n(2^-m)| (m <= 60) ==")
    HM = 80
    fam = closed_family_all(HM, -60, 0)
    Nt = {0: 1.0}
    for n in range(1, HM + 1):
        Nt[n] = sum(2.0 ** (-m) * abs(fam[n][-m]) for m in range(1, 61))
    for n in list(range(1, 13)) + list(range(16, HM + 1, 4)):
        print(f"   n={n:2d}: Ntilde_n = {Nt[n]:.3e}, Ntilde_n 3^(n/2) = {Nt[n] * 3 ** (n / 2):.3f}, (1/(log2 3 - 1))^n Ntilde_n = {Nt[n] / (LOG23 - 1) ** n:.3e}")
    ns = np.arange(20, HM + 1)
    ys = np.array([math.log(Nt[n]) for n in ns])
    slope, icpt = np.polyfit(ns, ys, 1)
    print(f"   least-squares rate of Ntilde_n over n = 20..80: {math.exp(slope):.4f} per level (3^-1/2 = {3 ** -0.5:.4f}, log2 3 - 1 = {LOG23 - 1:.4f}); prefactor {math.exp(icpt):.3f}")
    ys2 = np.array([math.log(Nt[n]) for n in range(40, HM + 1)])
    slope2, _ = np.polyfit(np.arange(40, HM + 1), ys2, 1)
    print(f"   rate over n = 40..80: {math.exp(slope2):.4f}; max_n Ntilde_n 3^(n/2) over 20..80 = {max(Nt[n] * 3 ** (n / 2) for n in range(20, HM + 1)):.3f}")
    print("== (b) per-J inequality |c_J| <= mass_h(J) Ntilde_(J-1) and (c) Chernoff mass law, h = 40, 80 ==")
    for h in (40, 80):
        s0 = int(round(h * LOG23)) - 6
        famh = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(famh, key=lambda j: abs(famh[j]))
        full = abs(famh[s_star])
        delta = h * LOG23 - s_star
        mass, contrib = decompose(h, s_star, "J")
        exact_mass = mass_law(h, s_star, 30)
        bound = 0.0
        worst = 0.0
        print(f"   h={h}: s* = {s_star} (delta = {delta:.2f}), |full| = {full:.4e}")
        for J in range(0, 31):
            Ntm = Nt[J - 1] if J >= 2 else 1.0
            cJ = abs(contrib[J])
            rhs = exact_mass[J] * Ntm
            bound += rhs
            ratio = cJ / rhs if rhs > 0 else 0.0
            worst = max(worst, ratio)
            chern = math.exp(-(h - J) * IRATE) * math.exp(-THETA * (s_star - (h - J) * LOG23))
            if J <= 12 or J % 3 == 0:
                print(f"      J={J:2d}: DP mass {mass[J]:.3e}, exact mass {exact_mass[J]:.3e}, Chernoff {min(1.0, chern):.3e}, |c_J| {cJ:.3e}, mass*Ntilde {rhs:.3e}, ratio {ratio:.3f}")
        print(f"   h={h}: worst |c_J|/(mass Ntilde) = {worst:.3f} (must be <= 1); bound sum = {bound:.4e}; |full|/bound = {full / bound:.4f}")
    print("== (d) the bound sum_J mass_h(J) Ntilde_(J-1) at h = 40..300 against the measured maxima M(h) (level-300 run) ==")
    M = {}
    with open("05-knowledge/results/collatz_five_mirrors_powers_of_two300_20260929.out", encoding="utf-8") as f:
        for line in f:
            m = re.match(r"\s*h=\s*(\d+):\s*([0-9.e+-]+)\s+argmax j=(-?\d+)", line)
            if m:
                M[int(m.group(1))] = (float(m.group(2)), int(m.group(3)))
    CN = max(Nt[n] * 3 ** (n / 2) for n in range(20, HM + 1))
    for h in (40, 60, 80, 100, 120, 150, 200, 250, 300):
        if h not in M:
            continue
        Mh, sh = M[h]
        JMAX = min(h, 120)
        em = mass_law(h, sh, JMAX)
        b = 0.0
        for J in range(0, JMAX + 1):
            if J <= 1:
                Ntm = 1.0
            elif J - 1 <= HM:
                Ntm = Nt[J - 1]
            else:
                Ntm = CN * 3 ** (-(J - 1) / 2)
            b += em[J] * Ntm
        b_exp = math.exp(-h * IRATE)
        print(f"   h={h:3d}: s = {sh} (delta = {h * LOG23 - sh:.2f}), M(h) = {Mh:.3e}, bound = {b:.3e}, M/bound = {Mh / b:.4f}, bound/e^(-hI) = {b / b_exp:.4f}, M/e^(-hI) = {Mh / b_exp:.4f}")
    print("DONE")
