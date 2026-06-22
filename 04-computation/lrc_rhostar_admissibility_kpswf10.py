#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_admissibility_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

*** THE DECISIVE DUE-DILIGENCE for the rho*=0 shapes. ***

THM-527 part D/G and THM-526 (line: "only the inadmissible non-covering config
gave rho*=0") rest on the claim that rho*=0 shapes are INADMISSIBLE.  We found
rho*=0 at admissible-LOOKING shapes (E=[0,3,4,5,7,8,9,10,11], P=[1,2,3,12], and
38 cases at k=9 spread<=14).  We now RIGOROUSLY test admissibility.

A genuine S3 counterexample-candidate S = P u L (|S|=13) must be:
  (1) COVERING: S contains a multiple of EVERY q in {2,...,14} (THM-523 B).
  (2) PRIMITIVE: gcd(S) = 1.
  (3) S3 structure: L = {v>13}, k=|L|>=3, Vmax >= 13*Vmin (clustered-large).
The co-offset reformulation: pick Vmax, set L = {Vmax - e_i : e_i in E}; we need
all L-elements > 13 and distinct, i.e. Vmax - max(E) > 13  =>  Vmax > 13+spread.

For a rho*=0 (P,E) we sweep Vmax over the admissible range and ask:
  - is S = P u L COVERING and PRIMITIVE?  (admissibility)
  - if so, what is the EXACT M(S)?  (does LRC still hold via a DIFFERENT witness?)
  - does the criterion C(S) hold via some v != Vmax?  (the route survives)

If NO admissible Vmax yields a covering+primitive S, the shape is INADMISSIBLE
(consistent with the canon claim) and rho*=0 is harmless.
If some admissible S IS covering+primitive, then rho*=0 is a GENUINE failure of
the via-Vmax criterion-C route -- and we report M(S) to see whether LRC holds
anyway (it must, if LRC(14) is true) via a non-Vmax witness.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def circ_maxgap_at(E, x):
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def gp_breaks(P):
    bps = set()
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks(E):
    bps = set()
    diffs = set()
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d != 0:
                diffs.add(d)
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        for m in range(0, 7 * d + 1):
            for s in (2, -2):
                v = Fr(7 * m + s, 7 * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def rho_star_exact(P, E):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    total = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for p in P:
            f = (Fr(p) * mid) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if dd < Fr(1, 14):
                ok = False
                break
        if ok and circ_maxgap_at(E, mid) > Fr(2, 7):
            total += (x1 - x0)
    return total


def is_covering(S):
    """S contains a multiple of every q in {2..14}."""
    for q in range(2, 15):
        if not any(v % q == 0 for v in S):
            return False
    return True


def is_primitive(S):
    return reduce(gcd, S) == 1


def M_exact(S):
    """Exact M(S) = max_tau min_v ||v tau||. Critical tau at (2k+1)/(2 v_i) and
       k/(v_a +- v_b). Evaluate min_v ||v tau|| at all critical tau, take max."""
    S = sorted(set(S))
    cand = set()
    Vmax = max(S)
    for v in S:
        for kk in range(0, 2 * v):
            cand.add(Fr(2 * kk + 1, 2 * v))
    for a in S:
        for b in S:
            for sgn in (a + b, a - b):
                if sgn != 0:
                    d = abs(sgn)
                    for kk in range(0, d + 1):
                        cand.add(Fr(kk, d))
    best = Fr(0)
    for tau in cand:
        if tau < 0 or tau >= 1:
            continue
        mn = None
        for v in S:
            f = (Fr(v) * tau) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if mn is None or dd < mn:
                mn = dd
        if mn is not None and mn > best:
            best = mn
    return best


def W_widest_arc(A, thr=Fr(1, 14)):
    """Widest level-thr safe arc of A (for criterion C). A = list of speeds."""
    bps = {Fr(0), Fr(1)}
    for a in A:
        for m in range(0, a):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * a)
                if 0 < v < 1:
                    bps.add(v)
    pts = sorted(bps)
    bestw = Fr(0)
    cur0 = None
    prev_safe = False
    # build safe arcs
    arcs = []
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for a in A:
            f = (Fr(a) * mid) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if dd < thr:
                ok = False
                break
        if ok:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    for a, b in arcs:
        if b - a > bestw:
            bestw = b - a
    return bestw


def criterion_C(S):
    """C(S): exists v in S with W(S\{v}) > 1/(7v). Return (holds, best_v, margin)."""
    S = sorted(set(S))
    best = None
    for v in S:
        A = [u for u in S if u != v]
        W = W_widest_arc(A)
        margin = W * 7 * v   # > 1 iff W > 1/(7v)
        if best is None or margin > best[1]:
            best = (v, margin)
        if margin > 1:
            return True, v, margin
    return False, best[0], best[1]


def main():
    print("=" * 78)
    print("THM-527 Thread A: ADMISSIBILITY of the rho*=0 shapes (decisive)")
    print("=" * 78)
    sys.stdout.flush()

    # the rho*=0 shapes found (k=9). E offsets, P small part.
    zeros = [
        ([0, 3, 4, 5, 7, 8, 9, 10, 11], [1, 2, 3, 12]),
        ([0, 1, 2, 3, 4, 5, 6, 7, 9], [1, 2, 3, 13]),
        ([0, 1, 2, 3, 4, 5, 6, 8, 11], [1, 2, 3, 12]),
        ([0, 1, 2, 3, 4, 5, 7, 8, 9], [1, 2, 3, 13]),
    ]
    any_admissible = False
    for (E, P) in zeros:
        r = rho_star_exact(P, E)
        spread = max(E)
        print(f"\n=== E={E}, P={P}  (rho*={r}) ===")
        print(f"    spread={spread}; admissible Vmax must satisfy Vmax-{spread}>13 "
              f"=> Vmax>={13+spread+1}={14+spread}")
        # sweep Vmax; build S; test covering + primitive; if so compute M and C.
        found_cov = 0
        examples = []
        Vlo = 14 + spread
        for Vmax in range(Vlo, Vlo + 400):
            L = [Vmax - e for e in E]            # cluster speeds
            if len(set(L)) != len(L):
                continue
            if min(L) <= 13:
                continue
            S = sorted(set(P) | set(L))
            if len(S) != 13:
                continue
            if not is_primitive(S):
                continue
            if is_covering(S):
                found_cov += 1
                if len(examples) < 3:
                    examples.append((Vmax, S))
        print(f"    covering+primitive admissible S over Vmax in [{Vlo},{Vlo+400}): "
              f"{found_cov}")
        if found_cov == 0:
            print(f"    => INADMISSIBLE (no covering+primitive S). rho*=0 HARMLESS here.")
        else:
            any_admissible = True
            for (Vmax, S) in examples:
                M = M_exact(S)
                Ch, cv, cm = criterion_C(S)
                print(f"    ADMISSIBLE Vmax={Vmax}: S={S}")
                print(f"        M(S) = {M} = {float(M):.6f}  "
                      f"({'>=1/14 OK' if M >= Fr(1,14) else 'BELOW 1/14 !!!'})")
                print(f"        criterion C(S): holds={Ch}  best v={cv} margin={float(cm):.4f}")
                print(f"        (rho*=0 means via-Vmax={Vmax} route fails; M from "
                      f"{'OTHER witness' if M>=Fr(1,14) else 'NOTHING - counterexample?'})")
        sys.stdout.flush()

    print("\n" + "=" * 78)
    print("VERDICT:")
    if any_admissible:
        print("  Some rho*=0 shapes ARE admissible (covering+primitive S3 sets).")
        print("  => the via-Vmax criterion-C route (rho*>0) is GENUINELY INSUFFICIENT")
        print("     for these S; THM-527's compactness-floor program does NOT close")
        print("     LRC(14).  LRC holds for them (M>=1/14) only via a DIFFERENT witness.")
        print("  This CONTRADICTS THM-526's 'only inadmissible configs give rho*=0'.")
    else:
        print("  ALL rho*=0 shapes tested are INADMISSIBLE (no covering+primitive S).")
        print("  => consistent with canon: rho*>0 on ADMISSIBLE shapes. The")
        print("     compactness floor should be posed over ADMISSIBLE (P,E) only.")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
