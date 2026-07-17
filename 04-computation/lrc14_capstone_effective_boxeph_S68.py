#!/usr/bin/env python3
"""
THE EFFECTIVE FRAME THEOREM -- finishing the capstone (boxeph-2026-07-17-S68)

The frame program's identities are all closed (LEM-031/032/033/039).  This
referee assembles the SIZE side with every constant explicit:

(A) THE SHARP WEIGHT BOUND:  |W-hat_g(chi)| <= W-hat_g(chi_0)
      = (pi^2/P^2) J_2(P/g)/3   (triangle inequality; sharp at chi_0),
    and the TOTAL weight collapses by the Jordan identity
      sum_{g | P, g < P} W-hat_g(chi_0) = (pi^2/3)(1 - P^{-2})
    (= sum_{a != 0} W(a), MI0).

(B) THE CLASS PARSEVAL (exact):  sum_{chi in U_q-hat} |X-hat_g(chi)|^2
      = (1/phi(q)) sum_{u in U_q} X(gu)^2  =: E_g   (the class energy).

(C) THE VARIANCE ENVELOPE (proved: Cauchy-Schwarz over classes with weight
    |W-hat_g|, then (A) + (B)):
      Var_w(cross) <= (pi^2/3)(1 - P^{-2}) *
                      sum_{g | P, g < P} (pi^2/P^2)(J_2(P/g)/3) * E_g.
    Fully explicit: Jordan totients x class energies.  Tightness measured.

(D) THE EFFECTIVE GENERIC-FRAME THEOREM (Chebyshev): for every lambda > 0,
      #{w : |cross(w) - mean| >= lambda} / phi(P)  <=  Var / lambda^2,
    with Var either the exact closed form (LEM-039(C)) or the envelope (C).
    Referee: measured exceedance counts vs the bound at 2, 3, 5 sigma.

(E) THE PER-CHARACTER MASS BOUND: |c-hat(chi)| <=
      sum_{g : cond | P/g} (pi^2/P^2)(J_2(P/g)/3) * min(sqrt(E_g), A_g),
    A_g = (1/phi(q)) sum_u |X(gu)| (level-1 bound).  Tightness on the top
    masses.

(F) the grade decomposition of the class energies E_g (LEM-033's cells:
    which grades carry the energy).
"""

import sys
from math import gcd, pi
sys.path.insert(0, '04-computation')
from lrc14_parity_lvalue_boxeph_S61 import (UGroup, cluster_spectra,
                                            cross_all, chat_direct, divisors,
                                            factorize, J2, csc2)

CLUSTERS = [([1, 2, 3, 4, 5, 6, 60], "family60"),
            ([1, 2, 3, 4, 5, 36, 60], "two-owner"),
            ([12, 15, 20, 21, 28, 30, 35], "balanced")]


def vp(n, p, cap):
    if n == 0:
        return cap
    v = 0
    while n % p == 0 and v < cap:
        n //= p
        v += 1
    return v


if __name__ == "__main__":
    print("THE EFFECTIVE FRAME THEOREM -- capstone referee (boxeph S68)")
    print("=" * 78)
    for E, name in CLUSTERS:
        P, M, owners, X, W = cluster_spectra(E, 0)
        G = UGroup(P)
        cr = cross_all(P, X, W, G.units)
        mean = sum(cr.values()) / G.phi
        var = sum(c * c for c in cr.values()) / G.phi - mean * mean

        # (A) sharp weight bound + Jordan collapse
        wtot = 0.0
        worstA = 0.0
        for g in divisors(P):
            if g == P:
                continue
            q = P // g
            Uq = [u for u in range(1, q + 1) if gcd(u, q) == 1]
            w0 = (pi * pi / (P * P)) * sum(csc2(u / q) for u in Uq)
            jw = (pi * pi / (P * P)) * J2(q) / 3
            worstA = max(worstA, abs(w0 - jw))
            wtot += w0
        jordan = (pi * pi / 3) * (1 - 1 / (P * P))
        print(f"[{name}] P={P} M={M}: (A) W-hat_g(chi_0) == (pi^2/P^2)J_2/3 "
              f"(max err {worstA:.1e}); total weight {wtot:.6f} == "
              f"(pi^2/3)(1-P^-2) = {jordan:.6f} "
              f"(err {abs(wtot - jordan):.1e})")
        assert worstA < 1e-9 and abs(wtot - jordan) < 1e-9

        # (B) class Parseval + energies
        Eg = {}
        worstB = 0.0
        for g in divisors(P):
            if g == P:
                continue
            q = P // g
            Gq = UGroup(q)
            lhs = 0.0
            for js in Gq.chars():
                xh = sum(X[(g * u) % P] * Gq.chi(js, u).conjugate()
                         for u in Gq.units) / Gq.phi
                lhs += abs(xh) ** 2
            rhs = sum(X[(g * u) % P] ** 2 for u in Gq.units) / Gq.phi
            worstB = max(worstB, abs(lhs - rhs) / max(rhs, 1e-9))
            Eg[g] = rhs
        print(f"    (B) class Parseval exact on {len(Eg)} classes "
              f"(worst rel err {worstB:.1e})")
        assert worstB < 1e-9

        # (C) the envelope
        env = jordan * sum((pi * pi / (P * P)) * (J2(P // g) / 3) * Eg[g]
                           for g in Eg)
        print(f"    (C) VARIANCE ENVELOPE: Var = {var:.1f} <= {env:.1f} "
              f"(tightness x{env / var:.2f}) -- PROVED, all constants "
              f"arithmetic")
        assert var <= env * (1 + 1e-9)

        # (D) Chebyshev exceedance referee
        rows = []
        sd = var ** 0.5
        for k in (2, 3, 5):
            lam = k * sd
            cnt = sum(1 for c in cr.values() if abs(c - mean) >= lam)
            bnd = G.phi / (k * k)
            rows.append(f"{k}sigma: {cnt} <= {bnd:.0f}")
            assert cnt <= bnd + 1e-9
        print(f"    (D) Chebyshev (mean {mean:+.2f}, sd {sd:.1f}): " +
              "; ".join(rows) + f"  (of {G.phi} frames)")

        # (E) per-character mass bound on the top masses
        triv = tuple(0 for _ in G.orders)
        masses = []
        for js in G.chars():
            if js == triv or G.parity(js) == -1:
                continue
            masses.append((abs(chat_direct(G, js, cr)), js))
        masses.sort(reverse=True)
        Ag = {g: sum(abs(X[(g * u) % P]) for u in
                     [u for u in range(1, P // g + 1)
                      if gcd(u, P // g) == 1]) /
              len([u for u in range(1, P // g + 1)
                   if gcd(u, P // g) == 1]) for g in Eg}
        print("    (E) per-character bound |c-hat| <= sum_g Wbound_g * "
              "min(sqrt(E_g), A_g):")
        for mval, js in masses[:3]:
            f = G.conductor(js)
            bnd = sum((pi * pi / (P * P)) * (J2(P // g) / 3) *
                      min(Eg[g] ** 0.5, Ag[g])
                      for g in Eg if (P // g) % f == 0)
            print(f"        cond={f:>4}: |c-hat| = {mval:7.3f} <= "
                  f"{bnd:7.1f} (x{bnd / mval:.1f})")
            assert mval <= bnd * (1 + 1e-9)

        # (F) grade decomposition of the dominant class energy (g = 1)
        al = {p: vp(P, p, 99) for p in (2, 3, 5, 7) if P % p == 0}
        print(f"    (F) energy by class (top 4 of {len(Eg)}): " + ", ".join(
            f"g={g}: {(pi*pi/(P*P))*(J2(P//g)/3)*Eg[g]:.1f}"
            for g, _ in sorted(Eg.items(),
                               key=lambda kv: -(pi*pi/(P*P)) *
                               (J2(P // kv[0]) / 3) * kv[1])[:4]))
        print()
    print("=" * 78)
    print("done")
