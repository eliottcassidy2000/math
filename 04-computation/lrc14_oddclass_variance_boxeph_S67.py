#!/usr/bin/env python3
"""
THE ODD CLASS LAW AND THE ASSEMBLED VARIANCE (boxeph-2026-07-17-S67)
Finishing sweep, part 1: LEM-038's two named opens.

(A) THE ODD CLASS LAW (LEM-039(A)).  For any modulus D | P let
    T_r(D) = sum_{p_i == r (mod D)} eps_i  (the signed class sums; the
    coincidence spectrum is N(D) = sum_r T_r^2).  The s = 3 reflection
    pairing (p, eps) <-> (P-p, -eps) gives, in one line:
        T_{-r}(D) = -T_r(D)          (the class spectrum is ODD),
    hence T_0 = 0, T_{D/2} = 0 (D even) -- self-mirror classes carry zero
    signed mass -- and N(D) = 2 sum_{half} T_r^2 is EVEN.  This is the
    mechanism behind S66(G)'s "N = M mostly": generic classes are mirror
    pairs of singletons (T = +-1).  GENERAL-s FORM: T_r^{(s)}(D) =
    -T_{-r}^{(6-s)}(D)  (the cross-section law).
(B) SUB-GRID SINE COHERENCE: S((P/D) nu) = 2i sum_{half} T_r sin(2 pi nu
    r / D) -- the imaginary law restricted to every divisor grid.
(C) THE ASSEMBLED VARIANCE (closing LEM-031..033 quantitatively):
        Var_w(cross) = sum_{chi != chi_0} | sum_g (2/g^2) L_{P/g}(2,chi)
                                            X-hat_g(chi) |^2
    -- the frame-cross variance reproduced ENTIRELY from the closed form
    (L-values x twisted coincidence sums), no direct frame sweep on the
    right side.  Referee on family60 and two-owner at s = 0.
(D) THE VARIANCE-DROP QUANTIFICATION (LEM-038 named (b)): per-conductor
    mass tables at s = 0 vs s = 3 (balanced): WHICH character masses
    shrink at the sine section.

Pure Python.  Reuses S61 (cluster_spectra, cross_all, UGroup, chat_direct,
chat_factorized, divisors) and S26 (owner_data).
"""

import sys
from math import pi, sin
import cmath

sys.path.insert(0, '04-computation')
from lrc14_parity_lvalue_boxeph_S61 import (UGroup, cluster_spectra,
                                            cross_all, chat_direct,
                                            chat_factorized, divisors)
from lrc14_general_resonance_law_boxeph_S26 import owner_data


def class_sums(E, s, D):
    P, M, data = owner_data(E, s)
    T = [0] * D
    for e, d in data.items():
        for sg, q in zip(d["sgn"], d["pos"]):
            T[q % D] += sg
    return P, M, T


CLUSTERS3 = [([12, 15, 20, 21, 28, 30, 35], "balanced"),
             ([1, 2, 3, 4, 5, 36, 60], "two-owner"),
             ([1, 2, 3, 4, 5, 6, 70], "family70"),
             ([8, 9, 10, 12, 14, 15, 18], "near-AP"),
             ([1, 2, 3, 4, 5, 56, 84], "two-large")]

if __name__ == "__main__":
    print("THE ODD CLASS LAW AND THE ASSEMBLED VARIANCE (boxeph S67)")
    print("=" * 78)
    print("PART A -- the odd class law at s = 3; the cross-section law s/6-s")
    for E, name in CLUSTERS3:
        P0, M0, _ = owner_data(E, 3)
        if M0 == 0:
            continue
        nchk = 0
        for D in divisors(P0):
            if D < 2 or D > 3000:
                continue
            P, M, T = class_sums(E, 3, D)
            for r in range(D):
                assert T[(D - r) % D] == -T[r], (name, D, r)
                nchk += 1
            assert T[0] == 0
            if D % 2 == 0:
                assert T[D // 2] == 0
            N = sum(t * t for t in T)
            assert N % 2 == 0
        # cross-section law at s = 1 vs 5 (sampled D)
        for s in (0, 1, 2):
            Ps, Ms, _ = owner_data(E, s)
            if Ms == 0:
                continue
            for D in divisors(Ps)[:8]:
                if D < 2:
                    continue
                _, _, T1 = class_sums(E, s, D)
                _, _, T2 = class_sums(E, 6 - s, D)
                for r in range(D):
                    assert T1[r] == -T2[(D - r) % D], (name, s, D, r)
                    nchk += 1
        print(f"  [{name}] odd class law + cross-section law: {nchk} exact "
              f"class checks; T_0 = T_D/2 = 0; N(D) even throughout")

    print()
    print("PART B -- sub-grid sine coherence")
    E, name = CLUSTERS3[0]
    P, M, data = owner_data(E, 3)
    for D in (7, 14, 20, 84, 140):
        if P % D:
            continue
        _, _, T = class_sums(E, 3, D)
        worst = 0.0
        for nu in range(1, D):
            S1 = sum(sg * cmath.exp(2j * pi * (P // D) * nu * q / P)
                     for e, d in data.items()
                     for sg, q in zip(d["sgn"], d["pos"]))
            S2 = 2j * sum(T[r] * sin(2 * pi * nu * r / D)
                          for r in range(1, (D + 1) // 2))
            worst = max(worst, abs(S1 - S2))
        print(f"  [{name}] D={D}: S((P/D)nu) == 2i sum T_r sin: max err "
              f"{worst:.1e}")
        assert worst < 1e-9 * max(M, 1)

    print()
    print("PART C -- the assembled variance (pure closed form)")
    for E, name in [([1, 2, 3, 4, 5, 6, 60], "family60"),
                    ([1, 2, 3, 4, 5, 36, 60], "two-owner")]:
        Pc, Mc, ownc, X, W = cluster_spectra(E, 0)
        G = UGroup(Pc)
        cr = cross_all(Pc, X, W, G.units)
        mean = sum(cr.values()) / G.phi
        var_direct = sum(c * c for c in cr.values()) / G.phi - mean * mean
        triv = tuple(0 for _ in G.orders)
        var_asm = 0.0
        for js in G.chars():
            if js == triv or G.parity(js) == -1:
                continue
            var_asm += abs(chat_factorized(G, js, Pc, X)) ** 2
        err = abs(var_asm - var_direct) / var_direct
        print(f"  [{name}] Var direct = {var_direct:.4f}; assembled from "
              f"L-values x X-hat = {var_asm:.4f} (rel err {err:.1e})")
        assert err < 1e-6

    print()
    print("PART D -- variance drop by conductor (balanced, s = 0 vs 3)")
    E = [12, 15, 20, 21, 28, 30, 35]
    for s in (0, 3):
        Pc, Mc, ownc, X, W = cluster_spectra(E, s)
        G = UGroup(Pc)
        cr = cross_all(Pc, X, W, G.units)
        mean = sum(cr.values()) / G.phi
        var = sum(c * c for c in cr.values()) / G.phi - mean * mean
        triv = tuple(0 for _ in G.orders)
        bycond = {}
        for js in G.chars():
            if js == triv:
                continue
            m2 = abs(chat_direct(G, js, cr)) ** 2
            f = G.conductor(js)
            bycond[f] = bycond.get(f, 0.0) + m2
        top = sorted(bycond.items(), key=lambda kv: -kv[1])[:6]
        print(f"  s={s}: M={Mc} var={var:9.1f}; mass by conductor: " +
              ", ".join(f"{f}: {v / var * 100:.1f}%" for f, v in top))
    print("=" * 78)
    print("done")
