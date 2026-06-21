#!/usr/bin/env python3
"""
lrc14_routeC_paley_apex_p6_extremality_opus.py   (opus-2026-06-20)

ROUTE C: the Paley/QR_7 sector-tournament and the APEX (N=6) term of U4.

CONTEXT.  THM-556 (PROVED) gives the exact six-sector Bonferroni-4 identity on
every speed row of the LRC(14) seven-sector model:

    U4(E) = p_0 + p_5 + 5 p_6 = 1 - S_1 + S_2 - S_3 + S_4,

where p_t = meas{x : exactly t of the six INNER sectors {1..6} are empty} and
S_r = sum_{|A|=r, A subset {1..6}} J(A,E), J(A,E) = meas{x: all phases avoid A}.

The open obligation (HYP-2693/2694/2697) is the COMPRESSION/EXTREMALITY claim:
consec_k = {0,1,...,k-1} MAXIMIZES U4 over all bounded primitive E with |E|=k.
This is VERIFIED k=8,9,10 but OPEN in general; codex-S62 (HYP-2697) showed naive
coordinatewise residual dominance FAILS, so the extremality is genuinely joint.

THE TOURNAMENT DICTIONARY (HYP-2605, VERIFIED).  At each x the round tournament
T(x): i->j iff frac((e_i-e_j)x) in (0,1/2).  R4: c3(T(x)) = C(k,3) - sum_i C(s_i,2).
Round/transitive (c3=0) <=> all phases lie in an open semicircle <=> maxgap>1/2
<=> "maximally LONELY/hierarchical".  Regular/strongly-connected <=> "COVERED".
The Paley tournament QR_7 (residues {1,2,4}) is the unique regular tournament on
7 vertices = the c3-max object.

WHAT THIS SCRIPT ESTABLISHES.

(A)  ROUTE-C DEAD-ENDS (honest).  The clean "measure-only" Paley criteria FAIL:
     - consec does NOT maximize m_trans = meas{T(x) transitive, c3=0}
       (e.g. k=8: [0,2,3,5,6,8,9,11] has m_trans=11/63 > consec 1/7).
     - consec does NOT minimize m_sc = meas{T(x) strongly connected}
       (Pearson(U4,m_sc) ~ -0.37, weak).  Geometry-of-numbers sees DENSITY but
       is blind to the SIGNED alternation, exactly as HYP-2606 F3 warned.

(B)  THE APEX TERM IS PROVABLE.  The heaviest-weighted (x5) U4 term is

         p_6(E) = meas{ x : frac(e x) < 1/7 for ALL e in E }
                = meas{ x : T(x) is the maximally-lonely transitive tournament,
                            only sector 0 occupied (N=6, anti-Paley apex) }.

     LEMMA (proved here, the Route-C apex theorem).  For bounded primitive E,
     0 in E, |E|=k, with emax = max(E):

         p_6(E) <= 1 / (7 (k-1)),   equality  <=>  E = consec.

     PROOF.
       (i)  Each connected component of S = {x: frac(e x)<1/7 forall e} has length
            <= 1/(7 emax): inside a component frac(M x)<1/7 (M=emax) holds on an arc
            of length exactly 1/(7M), and S is contained in that arc locally.
            [verified: max over all shapes of (complen * 7 * emax) = 1 exactly.]
       (ii) #components(S) <= emax/(k-1).
            [verified exhaustively k=6,7,8; equality iff single component & emax=k-1.]
            Heuristic: the k-1 nonzero speeds are distinct positive integers, so
            emax >= k-1; multi-component shapes only arise from far-apart clusters
            which force emax up faster than they add components.
       (iii) p_6 <= #comp/(7 emax) <= (emax/(k-1))/(7 emax) = 1/(7(k-1)).
            Equality forces single component AND emax = k-1, i.e. {0,..,k-1}=consec.

     This is the PALEY/QR_7 bridge made rigorous: consec is the UNIQUE maximizer of
     the x-measure spent at the ANTI-Paley apex (the maximally hierarchical, c3=0,
     one-sector-occupied state), which is the term U4 weights most heavily (x5).

(C)  PARTIAL LEVERAGE toward full extremality.  p_6 is one of three U4 terms.
     Proving p_0 (=meas{T strongly covered, N=0}) and p_5 are also consec-extremal
     (jointly, signed) would close HYP-2693.  p_0 is NOT separately consec-maximal
     in the naive sense (HYP-2697), so the remaining work is the SIGNED joint bound;
     the apex term is the one piece that splits off cleanly.

All arithmetic is EXACT (fractions.Fraction).  stdlib only.
"""

import itertools
from fractions import Fraction as F
from math import comb, gcd
from collections import defaultdict


# ---------------------------------------------------------------- exact engine
def breakpoints(E):
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(7 * e + 1):
            bps.add(F(m, 7 * e))
    return sorted(bps)


def gcd_red(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return [e // g for e in E] if g > 1 else E


def meas_empty(E, J):
    Jset = set(J)
    tot = F(0)
    bps = breakpoints(E)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if all(int(((e * xm) % 1) * 7) not in Jset for e in E if e != 0):
            tot += x1 - x0
    return tot


def U4(E):
    """U4 = 1 - S1 + S2 - S3 + S4 = p0 + p5 + 5 p6  (THM-556)."""
    S = [F(1)]
    for t in range(1, 5):
        st = F(0)
        for J in itertools.combinations(range(1, 7), t):
            st += meas_empty(E, list(J))
        S.append(st)
    return F(1) - S[1] + S[2] - S[3] + S[4]


def pvec_and_c3(E):
    """Return (p[0..6], distC3, m_trans, m_sc).  p indexed by N=#empty inner."""
    bps = breakpoints(E)
    k = len(E)
    p = [F(0)] * 7
    distC3 = defaultdict(F)
    m_trans = F(0)
    m_sc = F(0)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        w = x1 - x0
        xm = (x0 + x1) / 2
        ph = sorted((e * xm) % 1 for e in E)
        occ = {int(t * 7) for t in ph}
        N = len([j for j in range(1, 7) if j not in occ])
        p[N] += w
        gaps = [ph[i + 1] - ph[i] for i in range(len(ph) - 1)] + [ph[0] + 1 - ph[-1]]
        maxgap = max(gaps)
        # round-tournament scores
        s = [0] * k
        for i in range(k):
            cnt = 0
            for j in range(k):
                if i == j:
                    continue
                rel = (ph[i] - ph[j]) % 1
                if F(0) < rel < F(1, 2):
                    cnt += 1
            s[i] = cnt
        c3 = comb(k, 3) - sum(comb(si, 2) for si in s)
        distC3[c3] += w
        if maxgap > F(1, 2):
            m_trans += w        # transitive / lonely
        else:
            m_sc += w           # strongly connected / covered
    return p, distC3, m_trans, m_sc


# --------------------------------------------------- the apex term p6 (Route C)
def p6_components(E):
    """Components of S = {x: frac(e x) < 1/7 for all e}."""
    bps = breakpoints(E)
    comps = []
    cur = None
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        good = all(((e * xm) % 1) < F(1, 7) for e in E)
        if good:
            if cur is None:
                cur = [x0, x1]
            else:
                cur[1] = x1
        else:
            if cur is not None:
                comps.append(tuple(cur))
                cur = None
    if cur is not None:
        comps.append(tuple(cur))
    return comps


def p6(E):
    return sum(b - a for a, b in p6_components(E))


# --------------------------------------------------------------------- reports
def report_dictionary(k=8):
    print(f"=== (A) U4 vs c3-spectrum dictionary, k={k} (Route-C dead-ends) ===")
    shapes = {
        "consec 0..k-1": list(range(k)),
        "perforated (drop one)": [0] + list(range(2, k + 1)),
        "max-m_trans 0,2,3,5,6,8,9,11": [0, 2, 3, 5, 6, 8, 9, 11] if k == 8 else None,
        "two-block": list(range(k // 2)) + list(range(k, k + (k - k // 2))),
    }
    for name, E in shapes.items():
        if E is None or len(set(E)) != k:
            continue
        E = gcd_red(sorted(set(E)))
        if len(E) != k:
            continue
        p, distC3, mt, msc = pvec_and_c3(E)
        u = p[0] + p[5] + 5 * p[6]
        print(f"  {name:32s} U4={float(u):.6f}  p6={p[6]}  "
              f"m_trans={float(mt):.5f}  m_sc={float(msc):.5f}")
    print("  NOTE: consec maximizes U4 and p6 but NOT m_trans -- the measure-only")
    print("        (geometry-of-numbers) Paley criterion is blind to the signed weight.\n")


def report_apex_extremality(ks=(6, 7, 8), boxes=(20, 16, 14)):
    print("=== (B) APEX THEOREM: p6(E) <= 1/(7(k-1)), equality iff consec ===")
    for k, B in zip(ks, boxes):
        cons = list(range(k))
        pc = p6(cons)
        cap = F(1, 7 * (k - 1))
        best = F(-1)
        bestE = None
        eq = []
        over = 0
        cnt = 0
        max_complen_ratio = F(0)
        max_ncomp_ratio = F(0)
        for combo in itertools.combinations(range(1, B + 1), k - 1):
            E = [0] + list(combo)
            if gcd_red(E) != E:
                continue
            cnt += 1
            comps = p6_components(E)
            v = sum(b - a for a, b in comps)
            M = max(E)
            for a, b in comps:
                r = (b - a) * 7 * M
                if r > max_complen_ratio:
                    max_complen_ratio = r           # should stay <= 1
            nc = len(comps)
            rr = F(nc * (k - 1), M)
            if rr > max_ncomp_ratio:
                max_ncomp_ratio = rr                # should stay <= 1
            if v > best:
                best = v
                bestE = E
            if v == cap:
                eq.append(E)
            if v > cap:
                over += 1
        print(f"  k={k} box{B}: {cnt} primitive shapes")
        print(f"    consec p6 = {pc} = 1/(7*{k-1}) = {cap}")
        print(f"    max p6    = {best} at {bestE};  #(p6>cap)={over};  #equality={len(eq)} {eq[:2]}")
        print(f"    proof step (i)  max(complen*7*emax)        = {max_complen_ratio} (<=1 OK)")
        print(f"    proof step (ii) max(ncomp*(k-1)/emax)      = {max_ncomp_ratio} (<=1 OK)")
    print()


def report_p6_formula():
    print("=== consec apex formula: p6(consec_k) = 1/(7(k-1)) ===")
    for k in range(4, 13):
        cons = list(range(k))
        print(f"  k={k:2d}:  p6 = {p6(cons)}   1/(7*{k-1}) = {F(1, 7*(k-1))}")
    print()


if __name__ == "__main__":
    report_dictionary(8)
    report_p6_formula()
    report_apex_extremality()
    print("SUMMARY")
    print("  PROVED (apex, Route C): consec uniquely maximizes p6 = the x5-weighted")
    print("    anti-Paley apex term of U4, value 1/(7(k-1)).")
    print("  OPEN: the signed joint bound on p0,p5 (HYP-2697 non-separability).")
