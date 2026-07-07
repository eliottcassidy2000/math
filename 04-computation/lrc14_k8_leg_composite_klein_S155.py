#!/usr/bin/env python3
"""
klein-2026-07-07-S155 -- TARGET 1: the k=8 window-avoidance leg (HYP-2602 binding case),
attacked as a COMPOSITE mirroring the k=13 architecture one level down.

The k=8 leg (THM-530 union bound, bar = need_8 = thr_8 + m_P = 1 - 2243/5880 + 14249/252252):
  mu_{1/7}(E) >= need_8 ~ 0.67499  for every 8-element integer set E.

COMPOSITE:
 (a) diam <= D8*: PROVED by the subset lemma (kps-S59/monad-S2) + exact Farey roof:
     mu(E) >= mu(AP_{diam+1}); D8* = largest n-1 with mu(AP_n) >= need_8.
     THIS SCRIPT: exact crossings D_k* for ALL legs k=8..13 (new numbers for 8..12).
 (b) diam in [D8*+1, D0]: FINITE exact-screen enumeration over affine-normalized shapes
     (extends codex HYP-3530's span<=13 scan). Report per-diameter minima + margin.
 (c) '7 tame + 1 far': bisection identity (klein-S154) gives mu(E) = 1 - Bis_far
     (k=7 base PROVED); far-element crossing law Bis <= 1/7 + C*T0/gap. Measure the
     effective constant C_eff over tame-core x far-distance grids (the explicit-constant
     lemma's empirical side).
 (d) residual (all 7-subsets spread, no far split): adversarial mu minimization
     restricted OUT of (a)-(c) domains; plus direct rho* = meas(G_P n Good) and
     quasi-independence ratio R at k=8 over joint (P,E) adversaries -- measures how much
     the union bound wastes (THM-530-D said R >= 0.796 on its probes).
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

M_P = F(14249, 252252)
GP_MIN = {1: F(6,7), 2: F(66,91), 3: F(55,91), 4: F(1979,4004), 5: F(2243,5880),
          6: F(3029,10780), 7: F(45107,229320), 8: F(2479,17640), 9: F(10601,114660),
          10: F(14249,252252)}
def need_k(k):
    psz = 13 - k
    thr = F(0) if psz == 0 else 1 - GP_MIN[psz]
    return thr + M_P

def farey_pairs(n):
    a, b, c, d = 0, 1, 1, n
    while c <= n:
        yield (a, b, c, d)
        kk = (n + b) // d
        a, b, c, d = c, d, kk*c - a, kk*d - b

def mu_theta_AP_exact(n, theta: F):
    tot = F(0)
    for a, b, c, d in farey_pairs(n):
        w = F(1, b*d); y0, y1 = F(1, b), F(1, d)
        if y0 > y1: y0, y1 = y1, y0
        if y0 >= theta: tot += w
        elif y1 > theta: tot += w * (y1 - theta) / (y1 - y0)
    return tot

def grid(NG): return (np.arange(NG) + 0.5)/NG

def mu_batch(shapes, x, theta=1/7):
    """mu for a batch of equal-size shapes: shapes (B,k) int array, x (NG,)."""
    A = np.asarray(shapes, float)                      # (B,k)
    P = (A[:, None, :] * x[None, :, None]) % 1.0       # (B,NG,k)
    P.sort(axis=2)
    G = np.diff(P, axis=2)
    wrap = P[:, :, :1] + 1.0 - P[:, :, -1:]
    mg = np.maximum(G.max(axis=2), wrap[:, :, 0])
    return (mg > theta).mean(axis=1)                   # (B,)

def mu_one(E, x, theta=1/7):
    return float(mu_batch([E], x, theta)[0])

def primitive_normal(E):
    """affine-normalize: translate min->0, divide by gcd of diffs; canonical vs mirror."""
    E = sorted(E); E0 = [e - E[0] for e in E]
    g = reduce(gcd, E0[1:])
    E0 = tuple(e//g for e in E0)
    Em = tuple(sorted(E0[-1] - e for e in E0))
    return min(E0, Em)

if __name__ == "__main__":
    rng = np.random.default_rng(155155)
    TH = 1/7

    print("=== (a) EXACT per-leg tail-diameter crossings D_k* (bars need_k = thr_k + m_P) ===")
    Dstar = {}
    mu_table = {}
    for n in range(8, 100):
        mu_table[n] = mu_theta_AP_exact(n, F(1,7))
    for k in range(8, 14):
        nk = need_k(k)
        ns = [n for n in range(8, 100) if mu_table[n] >= nk]
        n_star = max(ns) if ns else None
        Dstar[k] = (n_star - 1) if n_star else None
        print(f"  k={k:>2}: need_k = {str(nk):>18} = {float(nk):.5f}   n_k* = {n_star}   D_k* = {Dstar[k]}"
              f"   [mu(AP_{n_star})={float(mu_table[n_star]):.5f} >= need > mu(AP_{n_star+1})={float(mu_table[n_star+1]):.5f}]")
    print(f"  (k=13 row should reproduce monad-S2's 75/76-77: D_13* = {Dstar[13]})")

    print("\n=== (b) FINITE BAND at k=8: exact-screen all affine-normalized shapes, diam in [D8*+1, D0] ===")
    D8 = Dstar[8]; D0 = 26
    xS = grid(4099); xM = grid(15013); xL = grid(46021)
    need8 = float(need_k(8))
    seen = set(); band_stats = {}
    worst_global = (2.0, None)
    for D in range(D8+1, D0+1):
        shapes = []
        for interior in combinations(range(1, D), 6):
            E = (0,) + interior + (D,)
            cn = primitive_normal(E)
            if cn in seen: continue
            seen.add(cn)
            if max(cn) != D: continue   # canonical rep sits at its own diameter
            shapes.append(cn)
        if not shapes:
            band_stats[D] = None; continue
        vals = np.concatenate([mu_batch(shapes[i:i+4096], xS) for i in range(0, len(shapes), 4096)])
        # refine the lowest 200 (or all below 0.75) on medium grid
        order = np.argsort(vals)
        cand_idx = [i for i in order[:200]] + [i for i in range(len(shapes)) if vals[i] < 0.75]
        cand_idx = sorted(set(cand_idx))
        best = (2.0, None)
        for i in cand_idx:
            v = mu_one(list(shapes[i]), xM)
            if v < best[0]: best = (v, shapes[i])
        vL = mu_one(list(best[1]), xL)
        band_stats[D] = (min(best[0], vL), best[1], len(shapes))
        if band_stats[D][0] < worst_global[0]:
            worst_global = (band_stats[D][0], best[1])
        print(f"  D={D:>2}: shapes={len(shapes):>6}  min mu = {band_stats[D][0]:.4f} at {best[1]}   margin over need_8: {band_stats[D][0]-need8:+.4f}")
    print(f"  BAND VERDICT: min over D in [{D8+1},{D0}] = {worst_global[0]:.4f} at {worst_global[1]}"
          f"  ({'ALL CLEAR of' if worst_global[0] >= need8 else 'BELOW'} need_8 = {need8:.5f})")

    print("\n=== (c) '7 tame + 1 far': effective crossing constant C_eff (Bis <= 1/7 + C_eff*T0/gap) ===")
    xF = grid(90007)
    for core_name, core in [("AP7 {0..6}", [0,1,2,3,4,5,6]),
                            ("tame7 {0,1,3,4,7,9,11}", [0,1,3,4,7,9,11]),
                            ("tame7 {0,2,3,7,8,12,15}", [0,2,3,7,8,12,15])]:
        T0 = max(core)
        print(f"  core={core_name} (T0={T0}):")
        for ratio in (3, 10, 40, 160):
            f = T0*ratio + 1
            E = sorted(core + [f])
            mu = mu_one(E, xF)
            bis = 1.0 - mu
            ceff = max(0.0, (bis - 1/7)) * (f - T0) / T0
            print(f"    far={f:>5} (gap/T0={ratio:>3}): mu = {mu:.4f}  Bis = {bis:.4f}  excess over 1/7 = {bis-1/7:+.4f}  C_eff = {ceff:.3f}")
    print("  (lemma shape: Bis <= Ind + C*T0/gap with Ind <= 1/7; crude provable C ~ 21+; measured C_eff below)")

    print("\n=== (d) RESIDUAL + the TRUE object at k=8: rho*, quasi-independence R over joint (P,E) ===")
    xR = grid(15013)
    def gp_mask(P, x):
        m = np.ones_like(x, bool)
        for p in P:
            m &= np.abs(((p*x) % 1.0) - 0.5) <= 0.5 - 1/14   # ||px|| >= 1/14
        return m
    def rho_R(P, E, x):
        g = gp_mask(P, x)
        A = np.asarray(E, float)
        Pts = np.sort((A[None, :] * x[:, None]) % 1.0, axis=1)
        G = np.diff(Pts, axis=1); wrap = Pts[:, :1] + 1.0 - Pts[:, -1:]
        mg = np.maximum(G.max(axis=1), wrap[:, 0])
        good = mg > 1/7
        mG = g.mean(); mu = good.mean(); rho = (g & good).mean()
        R = rho/(mG*mu) if mG*mu > 0 else np.inf
        return rho, mG, mu, R
    # joint adversarial: minimize rho (and R) over P (5-subsets of 1..13) and E (8-sets)
    worstRho = (2.0, None, None); worstR = (2.0, None, None)
    Plist = list(combinations(range(1, 14), 5))
    for trial in range(24):
        P = list(Plist[rng.integers(len(Plist))])
        H = int(rng.choice([10, 14, 20, 30]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        cur = rho_R(P, E, xR)[0]
        for step in range(50):
            if rng.random() < 0.35:
                P2 = list(Plist[rng.integers(len(Plist))]); E2 = E
            else:
                i = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 22, 34]))+1))
                if new in E: continue
                E2 = sorted(set(E) - {E[i]} | {new}); P2 = P
                if len(E2) != 8: continue
            c = rho_R(P2, E2, xR)[0]
            if c < cur - 1e-4: P, E, cur = P2, E2, c
        rho, mG, mu, Rr = rho_R(P, E, xR)
        if rho < worstRho[0]: worstRho = (rho, tuple(P), tuple(E))
        if Rr < worstR[0]: worstR = (Rr, tuple(P), tuple(E))
    print(f"  adversarial min rho* = {worstRho[0]:.4f} at P={worstRho[1]} E={worstRho[2]}  (bar m_P = {float(M_P):.4f})")
    rho, mG, mu, Rr = rho_R(list(worstRho[1]), list(worstRho[2]), grid(46021))
    print(f"    refined: rho*={rho:.4f} mG={mG:.4f} mu={mu:.4f} R={Rr:.3f}")
    print(f"  adversarial min R    = {worstR[0]:.3f} at P={worstR[1]} E={worstR[2]}  (THM-530-D probe floor was 0.796)")
    rho2, mG2, mu2, Rr2 = rho_R(list(worstR[1]), list(worstR[2]), grid(46021))
    print(f"    refined: rho*={rho2:.4f} mG={mG2:.4f} mu={mu2:.4f} R={Rr2:.3f}")
    print(f"  UNION-BOUND WASTE: union needs mu >= {need8:.4f}; if R >= R0 held uniformly, demand drops to")
    for R0 in (0.5, 0.75):
        print(f"    R0={R0}: mu >= m_P/(R0*mG_min) = {float(M_P)/(R0*float(GP_MIN[5])):.4f}")
