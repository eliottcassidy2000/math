#!/usr/bin/env python3
"""
klein-2026-07-07-S158 -- (1) the TRIPLE-MASS UPPER BOUND (conditioning lemma, measured
constant); (2) MIXED-LIFT EROSION; (3) the HARD-CORE DOSSIER.

(1) TRIPLE BOUND PLAN (HYP-4821): order q1 <= q2 <= q3 (positive, the endpoint's same-sign
    differences after gcd-reduction per pair is NOT available jointly -- work with raw d's,
    the law reads through reduced pairs). Conditioning on the smallest frequency:
      m123 = sum over the q1 intervals I (length w = theta/q1) of meas(I ∩ A2 ∩ A3).
    CONDITIONAL ERROR LEMMA (candidate): meas(I ∩ A2 ∩ A3) <= w*m23 + C*theta/q3,
    giving m123 <= theta*m23 + C*theta*q1/q3  [with m23 = theta^2 + G23/(q2'q3') the law].
    THIS SCRIPT: measure the sharp per-interval constant C = max_I (meas(I∩S) - w*m23)*q3/theta
    over many (q1,q2,q3) and shapes; then assemble S3 <= 35*theta^3 + 5*theta*(S2-21*theta^2)
    + C*theta*T(D) with T(D) = sum_triples d_(a)/d_(c), and find which diameter/spread makes
    Bonf3 = S2 - S3 >= 0.197 (the R-route bar) rigorous-given-C.

(2) MIXED-LIFT EROSION: screen c1*{1..a1} u c2*{1..a2} u B (c1<c2 in {2,3,4,5}, |F|=13,
    small B) + non-AP bases (near-AP perforations lifted) vs T* = 56291/294294; exact top-5.

(3) HARD-CORE DOSSIER (monad-S3's two shapes): P8 = {9,10,11,12,13} (k=8), P9 = {10,11,12,13}
    (k=9). Adversarial min over E of rho*(P,E), with R, mu, and the per-vertex W stats at the
    minimizers -- the concrete target numbers for opus's avoidance-kernel and mac-mini's
    PZ-on-U machinery, and first tournament-bridge test cases (mean theta-degree at k=8/9).
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

TH = 1/7
THF = F(1, 7)
TSTAR = F(1,7) + F(6,7)*F(14249,252252)
M_P = float(F(14249,252252))

def grid(NG): return (np.arange(NG)+0.5)/NG

def law_pair(d1, d2):
    g = gcd(d1, d2); q1, q2 = d1//g, d2//g
    r1, r2 = q1 % 7, q2 % 7
    return THF*THF + F(min(r1,r2)*(7-max(r1,r2)), 49*q1*q2)

if __name__ == "__main__":
    rng = np.random.default_rng(158158)
    x = grid(60013)

    print("=== (1a) conditional error constant C: meas(I ∩ A2 ∩ A3) - w*m23 <= C*theta/q3 ===")
    Cmax = 0.0; worst = None
    for trial in range(300):
        q1 = int(rng.integers(1, 15)); q2 = int(rng.integers(q1, 40)); q3 = int(rng.integers(q2, 120))
        if q2 == 0 or q3 == 0: continue
        m23 = float(law_pair(q2, q3)) if gcd(q2,q3) else 0
        m23 = float(law_pair(q2, q3))
        w = TH/q1
        # exact-ish numeric: measure I ∩ A2 ∩ A3 for each of the q1 intervals via fine local grids
        for j in range(q1):
            lo = j/q1
            t = lo + (np.arange(4001)+0.5)/4001*w
            in2 = (((q2*t) % 1.0) > 0) & (((q2*t) % 1.0) <= TH)
            in3 = (((q3*t) % 1.0) > 0) & (((q3*t) % 1.0) <= TH)
            m_loc = (in2 & in3).mean()*w
            Cval = (m_loc - w*m23) * q3 / TH
            if Cval > Cmax: Cmax, worst = Cval, (q1, q2, q3, j)
    print(f"  sharp measured C over 300 random triples x all q1-intervals: C = {Cmax:.3f} at (q1,q2,q3,j)={worst}")
    print(f"  => CONDITIONAL LEMMA candidate: meas(I∩A2∩A3) <= w*m23 + {np.ceil(Cmax*10)/10}*theta/q3")

    print("\n=== (1b) assembled Bonferroni-3 with the C-chain vs exact, on spread 8-sets ===")
    def S2S3_exact_numeric(D, xg):
        hits = np.stack([(((d*xg) % 1.0) > 0) & (((d*xg) % 1.0) <= TH) for d in D])
        S2 = sum((hits[a] & hits[b]).mean() for a in range(7) for b in range(a+1,7))
        S3 = sum((hits[a] & hits[b] & hits[c]).mean()
                 for a in range(7) for b in range(a+1,7) for c in range(b+1,7))
        return S2, S3
    C_use = max(2.0, np.ceil(Cmax*10)/10)
    for nm, E in [("spread", [0,5,11,17,26,33,41,50]), ("geometric", [0,3,8,17,31,52,80,118]),
                  ("mild", [0,2,5,9,14,20,27,35])]:
        et = max(E); D = sorted(et - e for e in E if e != et)
        S2, S3 = S2S3_exact_numeric(D, x)
        # chain bound: S3 <= sum_triples [theta*m23(mid,max) + C*theta*dmin/dmax]
        S3_chain = 0.0
        for a in range(7):
            for b in range(a+1,7):
                for c in range(b+1,7):
                    d1,d2,d3 = D[a],D[b],D[c]
                    S3_chain += TH*float(law_pair(d2,d3)) + C_use*TH*d1/d3
        print(f"  {nm:>10}: S2={S2:.4f} S3={S3:.4f} Bonf3={S2-S3:.4f} | chain S3<= {S3_chain:.4f} -> floor >= {S2-S3_chain:.4f} (bar 0.197)")

    print("\n=== (2) MIXED-LIFT EROSION vs T* = %.6f ===" % float(TSTAR))
    def Emg_mu(Efam, xg):
        A = np.asarray(Efam, float)
        P = np.sort((A[None,:]*xg[:,None]) % 1.0, axis=1)
        G = np.diff(P, axis=1); wrap = P[:, :1] + 1.0 - P[:, -1:]
        mg = np.maximum(G.max(axis=1), wrap[:,0])
        return float(mg.mean()), float((mg > TH).mean())
    def primitive(E):
        g = reduce(gcd, E); return tuple(sorted(e//g for e in E))
    cands = set()
    # mixed lifts c1*{1..a1} u c2*{1..a2} u B
    for c1 in (1,2,3):
        for c2 in (2,3,4,5):
            if c2 <= c1: continue
            for a1 in range(3, 12):
                for a2 in range(3, 12):
                    base = sorted(set([c1*i for i in range(1,a1+1)] + [c2*i for i in range(1,a2+1)]))
                    if not (10 <= len(base) <= 13): continue
                    nb = 13 - len(base)
                    if nb == 0:
                        cands.add(primitive(base)); continue
                    pool = [v for v in range(1, max(base)+3) if v not in base]
                    for B in combinations(pool, nb):
                        cands.add(primitive(sorted(base + list(B))))
    # perforated-AP bases 2-lifted: 2*({1..12}\{j}) u {two defects}
    for j in range(2, 12):
        base = [2*i for i in range(1, 13) if i != j]
        pool = [v for v in range(1, 27) if v not in base]
        for B in combinations(pool, 2):
            cands.add(primitive(sorted(base + list(B))))
    cands = [E for E in cands if len(E) == 13]
    print(f"  candidates: {len(cands)}")
    xs = grid(24001)
    vals = sorted((Emg_mu(list(E), xs)[0], E) for E in cands)
    print("  numeric leaders:")
    for emg, E in vals[:6]:
        print(f"    E[mg]~{emg:.6f}  {E}")
    rec = float(F(12907,65520))
    below = [v for v in vals if v[0] < rec - 5e-5]
    print(f"  candidates numerically below the record {rec:.6f}: {len(below)}"
          f"  -> mixed lifts {'ERODE' if below else 'DO NOT erode'} the record on this sweep")

    print("\n=== (3) HARD-CORE DOSSIER: P8={9..13} (k=8), P9={10..13} (k=9) ===")
    def gp_mask(P, xg):
        m = np.ones_like(xg, bool)
        for p in P: m &= np.abs(((p*xg) % 1.0) - 0.5) <= 0.5 - 1/14
        return m
    for (Pfix, k, nm) in [([9,10,11,12,13], 8, "P8={9..13}"), ([10,11,12,13], 9, "P9={10..13}")]:
        g = gp_mask(Pfix, x); mG = g.mean()
        wr = (2.0, None)
        for trial in range(30):
            H = int(rng.choice([k+2, 2*k, 3*k, 5*k]))
            E = sorted(rng.choice(np.arange(0, H+1), size=k, replace=False).tolist())
            def rho(EE, xg, gm):
                A = np.asarray(EE, float)
                Pp = np.sort((A[None,:]*xg[:,None]) % 1.0, axis=1)
                Gg = np.diff(Pp, axis=1); wrp = Pp[:, :1] + 1.0 - Pp[:, -1:]
                mg = np.maximum(Gg.max(axis=1), wrp[:,0])
                return float((gm & (mg > TH)).mean()), float((mg > TH).mean())
            xs2 = grid(8009); gs = gp_mask(Pfix, xs2)
            cur = rho(E, xs2, gs)[0]
            for step in range(40):
                i = int(rng.integers(k)); new = int(rng.integers(0, int(rng.choice([2*k, 4*k, 7*k]))+1))
                if new in E: continue
                cand = sorted(set(E) - {E[i]} | {new})
                if len(cand) != k: continue
                c = rho(cand, xs2, gs)[0]
                if c < cur - 1e-4: E, cur = cand, c
            v, mu = rho(E, x, g)
            if v < wr[0]: wr = (v, tuple(E), mu)
        v, Emin, mu = wr
        R = v/(mG*mu) if mG*mu > 0 else float('inf')
        print(f"  {nm}: meas(G_P) = {mG:.4f}; adversarial min rho* = {v:.4f} at E={Emin}")
        print(f"     mu(E) = {mu:.4f}; R = {R:.3f}; rho*/m_P = {v/M_P:.2f}x; mean theta-degree = {(k-1)/7.:.3f} ({'critical' if k==8 else 'supercritical'})")
