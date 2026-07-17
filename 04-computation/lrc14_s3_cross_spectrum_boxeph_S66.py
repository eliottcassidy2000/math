#!/usr/bin/env python3
"""
THE s = 3 CROSS-SPECTRUM (boxeph-2026-07-17-S66)
Owner directive: the s = 3 cross-spectrum (LEM-037's named open); find other
things to prove.

THE MASTER REFLECTION IDENTITY (LEM-038(A)).  The reflection x -> 1-x maps
the endpoint system of R_s to that of R_{6-s} (positions p -> P-p, signs
flip, owners preserved -- LEM-037), so for every owner e and every n:
    S_e^{(6-s)}(n) = - conj( S_e^{(s)}(n) ).
This subsumes LEM-037(D) (|S| equal) and specializes at the CENTRAL SECTION:

  (B) THE IMAGINARY LAW: Re S_e^{(3)}(n) = 0 for ALL owners and ALL n --
      the s = 3 system is a PURE SINE SYSTEM: S_e(n) = 2i sum_{pairs}
      eps sin(2 pi n p / P).
  (C) THE QUADRUPLE LAW: R_3 has no endpoint at 0 or 1/2 (at x = 1/2 an odd
      runner sits in the forbidden section 3; an all-even cluster collapses
      occupancy to {0}; at x = 0 occupancy is {0}), so endpoints pair with
      no fixed points and intervals pair I <-> 1-I with no self-symmetric
      interval: M(3) == 0 (mod 4).
  (D) THE PHASE-FREE CROSS: writing S_e^{(3)} = i a_e (a real), every owner
      cross-term S_e conj S_e' = a_e a_e' is REAL: the s = 3 cross-spectrum
      X(m) = sum_{e != e'} a_e(m) a_e'(m) is a Gram form -- interference
      without phases.
  (E) THE BASELINE COLLAPSE: sigma_e(3) = 0 kills the -(pi^2/3) sum sigma^2
      term of the LEM-030 baseline: the s = 3 cross mean is pure
      (K - 1/6)-torsion.  Measured next to s = 0..2 means.
  (F) THE CENSUS VERDICT: full character census at s = 3 (parity law holds
      verbatim -- its proof never used s), Parseval, variance; the
      "are s = 3 frames better" comparison: mean/var of cross, <Q_s>/M and
      frame-extremes of Q_s/M at s = 0 vs s = 3.
  (G) N(h) coincidence spectrum at s = 3 (measured; the pairing forces
      structure noted in the table).

Pure Python.  Reuses S61 (cluster_spectra, cross_all, UGroup, chat_direct)
and S26 (owner_data).
"""

import sys
from math import gcd, lcm, pi, sin
from fractions import Fraction as Fr
import cmath

sys.path.insert(0, '04-computation')
from lrc14_parity_lvalue_boxeph_S61 import (UGroup, cluster_spectra,
                                            cross_all, chat_direct, divisors)
from lrc14_general_resonance_law_boxeph_S26 import owner_data
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints

CLUSTERS = [([12, 15, 20, 21, 28, 30, 35], "balanced"),
            ([1, 2, 3, 4, 5, 36, 60], "two-owner"),
            ([1, 2, 3, 4, 5, 6, 70], "family70"),
            ([8, 9, 10, 12, 14, 15, 18], "near-AP")]


def S_all(E, s):
    """P, owners, per-owner spectra arrays S_e[n] (full Z_P), M."""
    P, M, data = owner_data(E, s)
    root = [cmath.exp(2j * pi * t / P) for t in range(P)]
    Se = {}
    for e, d in data.items():
        arr = [0j] * P
        for sg, q in zip(d["sgn"], d["pos"]):
            for n in range(P):
                arr[n] += sg * root[(n * q) % P]
        Se[e] = arr
    return P, M, data, Se


if __name__ == "__main__":
    print("THE s = 3 CROSS-SPECTRUM -- LEM-038 referee (boxeph S66)")
    print("=" * 78)
    print("PART A/B -- master reflection identity; the s = 3 imaginary law")
    for E, name in CLUSTERS:
        P0, M3, d3 = owner_data(E, 3)[0], *owner_data(E, 3)[1:3]
        if M3 == 0:
            print(f"  [{name}] R_3 empty; skipped")
            continue
        # master identity at s = 0..2 (sampled n), imaginary law at 3 (all n)
        worstA = 0.0
        for s in (0, 1, 2):
            Ps, Ms, ds = owner_data(E, s)
            Ps2, Ms2, ds2 = owner_data(E, 6 - s)
            if Ms == 0:
                continue
            for e in ds:
                for n in (1, 7, 11, Ps // 7 + 1, 3 * Ps // 5 + 1):
                    S1 = sum(sg * cmath.exp(2j * pi * n * q / Ps)
                             for sg, q in zip(ds[e]["sgn"], ds[e]["pos"]))
                    S2 = sum(sg * cmath.exp(2j * pi * n * q / Ps)
                             for sg, q in zip(ds2[e]["sgn"], ds2[e]["pos"]))
                    worstA = max(worstA, abs(S2 + S1.conjugate()))
        P, M, data, Se = S_all(E, 3)
        worstB = max(abs(Se[e][n].real) for e in Se for n in range(P))
        print(f"  [{name}] master identity (s=0,1,2 x owners x n): max err "
              f"{worstA:.1e}; s=3 imaginary law: max |Re S_e(n)| = "
              f"{worstB:.1e} over {len(Se)} owners x {P} freqs")
        assert worstA < 1e-9 and worstB < 1e-9 * max(M, 1)

    print()
    print("PART C -- the quadruple law: M(3) == 0 (mod 4); pairing fixed-point-free")
    for E, name in CLUSTERS + [([1, 2, 3, 4, 5, 36, 60], "two-owner-dup"),
                               ([1, 2, 3, 4, 5, 56, 84], "two-large"),
                               ([12, 15, 20, 21, 28, 30, 35], "balanced-dup")]:
        pts = endpoints(E, 3)
        if not pts:
            continue
        M = len(pts)
        pos = [p for p, sg, o in pts]
        assert Fr(0) not in pos and Fr(1, 2) not in pos, name
        # pairing: multiset {1-p} == {p} with signs/owners matched opposite
        look = {(p, sg, o) for p, sg, o in pts}
        for p, sg, o in pts:
            assert (1 - p if p != 0 else Fr(0), -sg, o) in look, (name, p)
        print(f"  [{name}] M(3) = {M} == 0 mod 4: {M % 4 == 0}; "
              f"no endpoint at 0 or 1/2; pairing exact")
        assert M % 4 == 0

    print()
    print("PART D/E -- phase-free cross; baseline collapse (means by s)")
    for E, name in CLUSTERS[:3]:
        P, M, data, Se = S_all(E, 3)
        owners = sorted(Se)
        # phase-free: X(m) two ways
        worst = 0.0
        for m in range(0, P, max(1, P // 500)):
            tot = sum(Se[e][m] for e in owners)
            X1 = abs(tot) ** 2 - sum(abs(Se[e][m]) ** 2 for e in owners)
            a = [Se[e][m].imag for e in owners]
            X2 = sum(a) ** 2 - sum(x * x for x in a)
            worst = max(worst, abs(X1 - X2))
            if len(owners) >= 2:
                worst = max(worst, max(abs((Se[e][m] *
                            Se[f][m].conjugate()).imag) for e in owners
                            for f in owners if f != e))
        print(f"  [{name}] s=3 phase-free (X = Gram of Im-parts; cross-terms "
              f"real): max err {worst:.1e}")
        assert worst < 1e-7 * M * M
        row = []
        for s in range(4):
            Ps, Ms, ds = owner_data(E, s)
            if Ms == 0:
                row.append(f"s={s}: empty")
                continue
            Pc, Mc, ownc, X, W = cluster_spectra(E, s)
            G = UGroup(Pc)
            cr = cross_all(Pc, X, W, G.units)
            mean = sum(cr.values()) / G.phi
            sig2 = (pi * pi / 3) * sum(sum(ds[e]["N"]) ** 2 for e in ds)
            row.append(f"s={s}: mean={mean:+.2f} (sigma^2-term={-sig2:+.2f})")
        print(f"      baselines: " + "  ".join(row))

    print()
    print("PART F -- the s=0 vs s=3 census: variance, masses, Q-frame stats")
    for E, name in CLUSTERS[:3]:
        stats = {}
        for s in (0, 3):
            Pc, Mc, ownc, X, W = cluster_spectra(E, s)
            G = UGroup(Pc)
            cr = cross_all(Pc, X, W, G.units)
            mean = sum(cr.values()) / G.phi
            var = sum(c * c for c in cr.values()) / G.phi - mean * mean
            # parity assert on 20 odd characters
            nodd = 0
            ptol = 1e-6 * max(max(abs(c) for c in cr.values()), 1e-9)
            for js in G.chars():
                if G.parity(js) == -1:
                    assert abs(chat_direct(G, js, cr)) < ptol, (name, s)
                    nodd += 1
                    if nodd >= 20:
                        break
            # Q_s over all frames from the total spectrum
            P2, M2, d2, Se = S_all(E, s)
            S2 = [abs(sum(Se[e][n] for e in Se)) ** 2 for n in range(P2)]
            csc = [0.0] + [(pi * pi / (P2 * P2)) / sin(pi * n / P2) ** 2
                           for n in range(1, P2)]
            Qs = []
            for w in G.units:
                q = 0.0
                nw = 0
                for n in range(1, P2):
                    nw += w
                    if nw >= P2:
                        nw -= P2 * (nw // P2)
                    q += csc[n] * S2[nw]
                Qs.append(q / M2)
            Qs.sort()
            stats[s] = (M2, mean, var, sum(Qs) / len(Qs), Qs[0], Qs[-1])
        r0, r3 = stats[0], stats[3]
        print(f"  [{name}] s=0: M={r0[0]:>3} mean={r0[1]:+8.2f} var={r0[2]:9.1f} "
              f"<Q>/M={r0[3]:6.2f} Q/M in [{r0[4]:.2f}, {r0[5]:.2f}]")
        print(f"  [{name.replace(name, ' ' * len(name))}] s=3: M={r3[0]:>3} "
              f"mean={r3[1]:+8.2f} var={r3[2]:9.1f} "
              f"<Q>/M={r3[3]:6.2f} Q/M in [{r3[4]:.2f}, {r3[5]:.2f}]")

    print()
    print("PART G -- N(h) at s = 3 (signed coincidence spectrum; mod-4 scan)")
    for E, name in CLUSTERS[:3]:
        P, M, data, Se = S_all(E, 3)
        pos = [q for e in data for q in data[e]["pos"]]
        sgn = [sg for e in data for sg in data[e]["sgn"]]
        vals = []
        for h in divisors(P):
            if h == 1 or h > 100:
                continue
            N = sum(si * sj for i, si in enumerate(sgn)
                    for j, sj in enumerate(sgn)
                    if (pos[i] - pos[j]) % (P // h) == 0)
            vals.append((h, N))
        print(f"  [{name}] M={M}; (h, N(P/h-classes)): " + ", ".join(
            f"({h},{N})" for h, N in vals[:10]) +
            f"; all N == M mod 2: {all((N - M) % 2 == 0 for _, N in vals)}; "
            f"all N == 0 mod 4: {all(N % 4 == 0 for _, N in vals)}")
    print("=" * 78)
    print("done")
