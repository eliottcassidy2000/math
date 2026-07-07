#!/usr/bin/env python3
"""
klein-2026-07-07-S154 -- THE BISECTION RECURSION for the density floor mu_theta.

EXACT IDENTITY (elementary; proof in the reflection):
  For any finite integer set E, any e_j in E, any theta in (0,1):
    N_theta(E) := {x : maxgap{frac(e x) : e in E} <= theta}
  decomposes as the DISJOINT union
    N_theta(E) = N_theta(E_j)  u  BIS_j(theta,E),   E_j := E \ {e_j},
  where BIS_j := {x : E_j has exactly one gap > theta (necessarily <= 2*theta; it is the
  union of the two E-gaps abutting frac(e_j x)), all other gaps <= theta, and frac(e_j x)
  lies in that gap's theta-middle window (length 2*theta - G)}.
  Proof: adding a point refines gaps, so N_theta(E_j) subset N_theta(E); conversely on
  N_theta(E) every E_j-gap except possibly the P_j-containing one is an E-gap (<= theta),
  and the P_j-one equals g_left+g_right <= 2*theta, split into two <= theta parts iff P_j
  lies in the middle window.  Both inclusions compose to the disjoint decomposition. QED

  COROLLARIES:
    (C1) mu_theta(E) = mu_theta(E_j) - Bis_j(theta,E)   for EVERY j  [exact recursion]
    (C2) mu_theta(E) >= max_j mu_{2 theta}(E_j)                      [containment floor]
    (C3) (A')-induction: AP-minimality at k reduces, given k-1, to
         exists j: mu(E_j) >= mu(AP_{k-1})  AND  Bis_j(E) <= Bis_end(AP_k)
         where Bis_end(AP_k) = mu_theta(AP_{k-1}) - mu_theta(AP_k)
         (at theta=1/7, k=13: = 13823/24255 - 477/1078 = 6181/48510 ~ 0.12742).
    (C4) floor-only: mu_theta(E) >= m_P needs only exists j: mu_theta(E_j) - Bis_j >= m_P.

THIS SCRIPT (numpy, two grids for stability; survey-numeric, not canon-exact):
  1. identity validation; 2. Bis structure at theta=1/7 k=13 (AP saturates bisection);
  3. Delta_j = Bis_j - Ind_j (independence excess) vs midpoint-relation count R3_j;
  4. C2 leave-one-out 2/7 floor: adversarial min of max_j mu_{2/7}(E_j) vs m_P;
  5. k=8 leg from the PROVED k=7 base: floor 1 - min_j Bis_j vs thr_8 (jump-move descent).
"""
import numpy as np
from math import gcd
from functools import reduce

M_P = 14249/252252
THR8 = 3637/5880

def grid(NG):
    return (np.arange(NG) + 0.5) / NG

def gapmat(E, x):
    """orbit points (NG, k) sorted, and gaps (NG, k) with wraparound."""
    P = np.sort(np.outer(x, np.asarray(E, float)) % 1.0, axis=1)
    G = np.diff(P, axis=1, append=(P[:, :1] + 1.0))
    return P, G

def mu_theta(E, theta, x):
    _, G = gapmat(E, x)
    return float(np.mean(G.max(axis=1) > theta))

def bis_parts(E, j_elem, theta, x):
    """returns (mu_E, mu_Ej, Bis_j, Ind_j, mu2_Ej) on grid x."""
    Ej = [e for e in E if e != j_elem]
    Pj, Gj = gapmat(Ej, x)
    mgj = Gj.max(axis=1)
    mu_Ej = float(np.mean(mgj > theta))
    mu2_Ej = float(np.mean(mgj > 2*theta))
    big = Gj > theta
    nbig = big.sum(axis=1)
    # unique big gap <= 2 theta
    sel = (nbig == 1) & (mgj <= 2*theta)
    idx = np.argmax(big, axis=1)
    rows = np.arange(len(x))
    gstart = Pj[rows, idx]           # start of the big gap
    glen = Gj[rows, idx]
    win_lo = (gstart + glen - theta) % 1.0
    win_len = 2*theta - glen
    pj = (np.asarray(j_elem, float) * x) % 1.0
    inwin = ((pj - win_lo) % 1.0) < win_len
    bis = float(np.mean(sel & inwin))
    ind = float(np.mean(np.where(sel, win_len, 0.0)))
    # mu_E directly
    _, GE = gapmat(E, x)
    mu_E = float(np.mean(GE.max(axis=1) > theta))
    return mu_E, mu_Ej, bis, ind, mu2_Ej

def R3_through(E, j_elem):
    S = set(E); c = 0
    for a in E:
        if a == j_elem: continue
        b = 2*j_elem - a
        if b in S and b != j_elem and b > a: c += 1
    return c

def primitive(E):
    g = reduce(gcd, E)
    return tuple(sorted(e//g for e in E))

BANK13 = {
    "AP {1..13}":                 list(range(1, 14)),
    "prim-sat 2{1..12}u{13}":     sorted([2*i for i in range(1, 13)] + [13]),
    "monad rec 2{1..11}u{11,13}": sorted([2*i for i in range(1, 12)] + [11, 13]),
    "dstar min":                  [1,3,5,6,7,8,9,10,11,13,15,20,29],
    "spread {6..58}":             [6,10,14,18,22,26,30,34,38,42,46,50,58],
    "primes<=41":                 [2,3,5,7,11,13,17,19,23,29,31,37,41],
    "fibonacci-ish":              [1,2,3,5,8,13,21,34,55,89,144,233,377],
}

if __name__ == "__main__":
    rng = np.random.default_rng(20260707)
    TH = 1/7
    XA, XB = grid(30011), grid(77003)   # two coprime-ish grids for stability

    print("=== 1. IDENTITY VALIDATION: mu(E) = mu(E_j) - Bis_j (every j) ===")
    worst = 0.0
    for nm, E in BANK13.items():
        for j in (E[0], E[len(E)//2], E[-1]):
            muE, muEj, bis, ind, _ = bis_parts(E, j, TH, XA)
            worst = max(worst, abs(muE - (muEj - bis)))
    print(f"  max |mu(E) - (mu(E_j)-Bis_j)| over bank x 3 j: {worst:.2e}  (grid-exact 0 expected)")

    print("\n=== 2. BISECTION STRUCTURE theta=1/7, k=13 (grids 30011 / 77003) ===")
    print(f"{'family':>28} {'mu(E)':>7} | best j: {'mu(Ej)':>7} {'Bis_j':>7} {'Ind_j':>7} {'Delta':>7} {'R3_j':>4} | {'C2 max_j mu2/7(Ej)':>18}")
    for nm, E in BANK13.items():
        best = None; c2 = 0.0
        for j in E:
            muE, muEj, bis, ind, mu2Ej = bis_parts(E, j, TH, XA)
            c2 = max(c2, mu2Ej)
            fl = muEj - bis
            if best is None or fl > best[0]:
                best = (fl, j, muE, muEj, bis, ind)
        fl, j, muE, muEj, bis, ind = best
        muB = mu_theta(E, TH, XB)
        print(f"{nm:>28} {muE:7.4f} |  j={j:>3}: {muEj:7.4f} {bis:7.4f} {ind:7.4f} {bis-ind:+7.4f} {R3_through(E,j):>4} | {c2:18.4f}   [gridB mu={muB:.4f}]")

    print("\n  -- AP: endpoint vs middle drop (alignment mechanism; exact Bis_end(AP13)=6181/48510=0.12742) --")
    E = list(range(1, 14))
    for j in (13, 7, 1, 12):
        muE, muEj, bis, ind, _ = bis_parts(E, j, TH, XA)
        print(f"   drop j={j:>2}: mu(E_j)={muEj:.4f}  Bis_j={bis:.4f}  Ind_j={ind:.4f}  Delta={bis-ind:+.4f}  R3_j={R3_through(E,j)}")

    print("\n=== 3. Delta_j vs midpoint count R3_j (mechanism), 60 random (E,j) ===")
    pairs = []
    for _ in range(60):
        E = sorted(rng.choice(np.arange(1, 45), size=13, replace=False).tolist())
        j = int(rng.choice(E))
        _, _, bis, ind, _ = bis_parts(E, j, TH, XA[:15005])
        pairs.append((R3_through(E, j), bis - ind))
    r = np.array([p[0] for p in pairs], float); d = np.array([p[1] for p in pairs])
    corr = float(np.corrcoef(r, d)[0, 1])
    print(f"  corr(R3_j, Delta_j) = {corr:+.3f}")
    for band, lab in [((r == 0), "R3=0"), ((r >= 1) & (r <= 2), "R3=1-2"), ((r >= 3), "R3>=3")]:
        if band.any():
            print(f"    mean Delta | {lab:>6}: {d[band].mean():+.4f}   (n={band.sum()})")

    print("\n=== 4. C2 LEAVE-ONE-OUT 2/7 FLOOR: minimize max_j mu_{2/7}(E_j) (13-sets) ===")
    print(f"   target: stay >= m_P = {M_P:.4f} (then C2 alone gives the k=13 witness floor)")
    def loo27(E, x):
        return max(mu_theta([e for e in E if e != j], 2/7, x) for j in E)
    for nm, E in BANK13.items():
        print(f"   {nm:>28}: {loo27(E, XA[:15005]):.4f}")
    XS = grid(6007)
    gbest = (2.0, None)
    for trial in range(30):
        H = int(rng.choice([16, 20, 28, 40]))
        E = sorted(rng.choice(np.arange(1, H+1), size=13, replace=False).tolist())
        cur = loo27(E, XS)
        for step in range(60):
            i = int(rng.integers(13)); new = int(rng.integers(1, int(rng.choice([16, 30, 44]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 13: continue
            c = loo27(cand, XS)
            if c < cur - 1e-3: E, cur = cand, c
        v = loo27(E, XA[:15005])
        if v < gbest[0]: gbest = (v, primitive(E))
    vb = loo27(list(gbest[1]), XB[:38501])
    print(f"   ADVERSARIAL MIN of max_j mu_2/7(E_j): {gbest[0]:.4f} (gridB {vb:.4f}) at {gbest[1]}")
    print(f"   -> C2-alone k=13 floor {'HOLDS' if min(gbest[0],vb) >= M_P else 'FAILS'} on this bank (m_P={M_P:.4f})")

    print("\n=== 5. k=8 LEG from PROVED k=7 base: mu_1/7(E_8) = 1 - Bis_j, need >= thr_8 ===")
    print(f"   thr_8 = {THR8:.4f}  <=>  exists j with Bis_j <= {1-THR8:.4f}")
    def minbis8(E, x):
        best = None
        for j in E:
            _, _, bis, _, _ = bis_parts(E, j, TH, x)
            if best is None or bis < best: best = bis
        return best
    bank8 = {
        "consec {1..8}":        list(range(1, 9)),
        "perforated {1,3..9}":  [1,3,4,5,6,7,8,9],
        "2{1..7}u{9}":          [2,4,6,8,9,10,12,14],
        "spread":               [1,6,12,18,27,34,42,51],
    }
    for nm, E in bank8.items():
        mb = minbis8(E, XA[:15005])
        muE = mu_theta(E, TH, XA[:15005])
        print(f"   {nm:>22}: min_j Bis_j = {mb:.4f} -> floor {1-mb:.4f}  (true mu={muE:.4f})  [{'OK' if 1-mb >= THR8 else 'below'}]")
    gworst = (0.0, None)
    for trial in range(30):
        H = int(rng.choice([9, 12, 16, 24]))
        E = sorted(rng.choice(np.arange(1, H+1), size=8, replace=False).tolist())
        cur = minbis8(E, XS)
        for step in range(60):
            i = int(rng.integers(8)); new = int(rng.integers(1, int(rng.choice([12, 20, 30]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 8: continue
            c = minbis8(cand, XS)
            if c > cur + 1e-3: E, cur = cand, c
        v = minbis8(E, XA[:15005])
        if v > gworst[0]: gworst = (v, primitive(E))
    vb = minbis8(list(gworst[1]), XB[:38501])
    print(f"   ADVERSARIAL MAX of min_j Bis_j (8-sets): {gworst[0]:.4f} (gridB {vb:.4f}) at {gworst[1]}")
    print(f"   -> uniform k=8 recursion floor 1 - maxminBis = {1-max(gworst[0],vb):.4f} vs thr_8 {THR8:.4f}"
          f"  ({'HOLDS' if 1-max(gworst[0],vb) >= THR8 else 'FAILS'} on this bank)")
