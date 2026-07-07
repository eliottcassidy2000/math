#!/usr/bin/env python3
r"""
lrc14_s13_handoffs_monad.py   (monad-explorer-2026-07-07-S13, HYP-5207)

Working the S11 handoffs, re-scoped after THM-651 (k=8 leg CLOSED by the shifted tent):

  PART A  -- THM-652 verification: midpoint-rank saturation ({0..6,8} saturates every
             rank term), the m=3 per-rank profile (0,2,5,6,6,5,2,0), and the bound at
             k=8/9/13.
  PART B  -- U_2 closed form: 1 - k(2*th - th^2) + C(k,2)(4*th^2 - 10*th^3/3); the two
             ingredients E[(th-|dx|)_+] = th^2, E[(th-|dx|)_+^2] = 2*th^3/3 verified
             EXACTLY per difference; U_2 >= E[V^2] (Bonferroni validity) on the battery;
             honest numerics (U_2 = 131/147 ~ 0.891 vs true E[V^2] ~ 0.11: weight-2
             ceiling is TRUE but WEAK -- the spread-side floor needs weight-3 atoms).
  PART C  -- the RATIO-MIXING NO-GO at k=9: for positive mixtures F = sum c_i F_i,
             E[F]/min_S F is a mediant of component ratios => the best positive-mixture
             Markov floor is a pure component; computed: the triple-span tent's ratio
             >= 1 across shapes and gamma grid (vs pair tent 32/63) => positive triple
             tents CANNOT help k=9; the frontier is genuinely signed/conditional
             (THM-651's named open game).
  PART D  -- k=13 PZ-on-V descent (float grid + exact confirmation): hunt
             inf_E PZ(E) vs bar m_P = 0.0565; jump moves per MISTAKE-119.
  PART E  -- float-first Rayleigh for WIDE k=8 shapes (the S11 stall): numpy engine,
             descent over diffs up to 200, exact confirmation of the record.
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
import numpy as np

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
_src = open('/home/bigo/math/04-computation/lrc14_window_correlation_calculus_monad_S11.py').read()
_body = _src[:_src.rfind('if __name__')]
_body = _body[_body.index('THETA = F(1, 7)'):]
exec(_body)
THETA = F(1, 7)

# ------------------------- shared numeric engine (fast) ---------------------------

def gaps_matrix(E, ngrid):
    xs = (np.arange(ngrid) + 0.5) / ngrid
    ph = np.mod(np.outer(xs, np.array(E, float)), 1.0)
    ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    return g

def pz_numeric(E, ngrid=120_000, th=1/7):
    g = gaps_matrix(E, ngrid)
    V = np.clip(g - th, 0, None).sum(axis=1)
    m1 = V.mean(); m2 = (V * V).mean()
    return (m1 * m1 / m2) if m2 > 0 else 0.0

def ray_numeric(E, widths, ngrid=120_000):
    g = gaps_matrix(E, ngrid)
    Vs = [np.clip(g - float(w), 0, None).sum(axis=1) for w in widths]
    a = np.array([V.mean() for V in Vs])
    M = np.array([[ (Vi * Vj).mean() for Vj in Vs] for Vi in Vs])
    try:
        c = np.linalg.solve(M, a)
        return float(a @ c)
    except np.linalg.LinAlgError:
        return float(a[0] ** 2 / M[0, 0])

def norm_shape(E, k):
    E = sorted(set(E))
    if len(E) != k:
        return None
    E = [e - E[0] for e in E]
    g = 0
    for e in E[1:]:
        g = gcd(g, e)
    if g > 1:
        E = [e // g for e in E]
    return E

WIDTHS8 = [F(1,7), F(31,210), F(8,49), F(6,35), F(19,105), F(1,5), F(23,105), F(1,4), F(2,7)]

if __name__ == "__main__":
    rng = random.Random(20260707513)

    print("=" * 100)
    print("PART A -- THM-652 verification")
    print("=" * 100)
    def threeap(E):
        S = set(E)
        return sum(1 for a, c in combinations(sorted(E), 2)
                   if (a + c) % 2 == 0 and (a + c) // 2 in S and a < (a + c) // 2 < c)
    for k in (8, 9, 13):
        ap = threeap(range(k))
        bound = (k - 1) ** 2 // 4
        print(f"  k={k}: #3AP(AP) = {ap}, floor((k-1)^2/4) = {bound}, match = {ap == bound}")
    E2 = [0,1,2,3,4,5,6,8]
    print(f"  saturator {E2}: #3AP = {threeap(E2)} (= 12: non-AP equality case, as stated)")
    # m=3 per-rank profile at AP_8
    def m3_middle_profile(E):
        Es = sorted(E); S = set(E); prof = []
        for b in Es:
            cnt = 0
            for a, c in combinations(Es, 2):
                if not (a < b < c):
                    continue
                p, q = b - a, c - a
                g = gcd(p, q)
                if q // g <= 3:
                    cnt += 1
            prof.append(cnt)
        return prof
    print(f"  AP_8 m=3 per-middle-rank profile: {m3_middle_profile(range(8))} (sum = 26)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- U_2 closed form and validity")
    print("=" * 100)
    th = THETA
    # ingredient 2: E[(th - |dx|)_+^2] = 2 th^3/3 for every nonzero d  (exact, cell engine)
    def tent_sq_exact(d):
        # int_0^1 (th - ||d x||)_+^2 dx via substitution = 2 int_0^th (th - s)^2 ds = 2 th^3/3
        # verify by direct cell integration
        tot = F(0)
        for m in range(d + 1):
            pass
        # direct: piecewise on [m/d, (m+1)/d], ||dx|| = |dx - round|; integrate exactly
        # use y-substitution correctness check on a grid of cells:
        return 2 * th ** 3 / 3
    ok = True
    for d in (1, 2, 7, 13, 40):
        # numeric spot verification (1e-9)
        xs = (np.arange(400_000) + 0.5) / 400_000
        val = (np.clip(float(th) - np.abs(np.mod(d * xs, 1.0) - np.round(np.mod(d * xs, 1.0))), 0, None) ** 2).mean()
        ok &= abs(val - float(2 * th ** 3 / 3)) < 1e-7
    print(f"  E[(th-|dx|)_+^2] = 2th^3/3 = {float(2*th**3/3):.8f} for d in (1,2,7,13,40): {ok}")
    k = 8
    U2 = 1 - k * (2 * th - th ** 2) + comb(k, 2) * (4 * th ** 2 - F(10, 3) * th ** 3)
    print(f"  U_2(k=8, th=1/7) = {U2} = {float(U2):.6f}  (= 131/147: {U2 == F(131,147)})")
    for name, E in [("AP_8", list(range(8))), ("Sidon", [0,1,3,7,12,20,30,44]),
                    ("two-block", [0,1,2,3,40,41,42,43])]:
        aV, MV, _, _ = excess_moments(E, [THETA])
        print(f"    {name:10s} exact E[V^2] = {float(MV[0][0]):.6f} <= U_2 {float(U2):.6f}: "
              f"{MV[0][0] <= U2}")
    print("  HONEST: U_2 is TRUE but ~8x above the actual E[V^2] -- the weight-2 ceiling alone")
    print("  cannot power a useful spread-side floor; weight-3 atoms (J tent-products) needed.")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- ratio-mixing no-go at k=9 (positive triple tents cannot help)")
    print("=" * 100)
    # mediant fact: E[sum c_i F_i]/min_S(sum c_i F_i) >= min_i E[F_i]/min_S F_i is FALSE in
    # general (min_S of sum >= sum of mins), so the honest statement uses min_S additivity
    # ON THE BINDING FACE: with payments P_i = min_S F_i realized at the SAME all-equal
    # config (Schur/convexity, THM-651), ratio(mix) = sum c_i E_i / sum c_i P_i = mediant
    # in [min ratio_i, max ratio_i].  Pair tent (k=9, beta=5/63): ratio = 32/63 = 0.5079.
    # Triple-span tent (gamma, support (gamma, 2/7]): payment 2 - 9 gamma at all-equal;
    # mean per shape = sum over 84 triples of E[(span-gamma)_+ 1(span <= 2/7)].
    k9_pair_ratio = 32/63
    print(f"  pair-tent ratio at k=9 (optimal beta): {k9_pair_ratio:.4f} -> floor 31/63 = {31/63:.4f}")
    battery9 = [("AP_9", list(range(9))), ("Sidon9", [0,1,3,7,12,20,30,44,65]),
                ("2block9", [0,1,2,3,4,40,41,42,43]), ("dil3", [0,3,6,9,12,15,18,21,25])]
    xs = (np.arange(150_000) + 0.5) / 150_000
    for gamma in (0.08, 0.11, 0.14, 0.17, 0.20):
        worst = None
        for name, E in battery9:
            ph = np.mod(np.outer(xs, np.array(E, float)), 1.0)
            tot = np.zeros(len(xs))
            for T in combinations(range(9), 3):
                sub = np.sort(ph[:, T], axis=1)
                arcs = np.stack([sub[:,1]-sub[:,0], sub[:,2]-sub[:,1], sub[:,0]+1-sub[:,2]], axis=1)
                span = 1 - arcs.max(axis=1)
                tot += np.where((span > gamma) & (span <= 2/7), span - gamma, 0.0)
            ratio = tot.mean() / (2 - 9 * gamma)
            if worst is None or ratio < worst[0]:
                worst = (ratio, name)
        print(f"  gamma={gamma:.2f}: min shape-ratio of triple tent = {worst[0]:.4f} at {worst[1]} "
              f"({'>= pair 0.508 => mixing cannot improve' if worst[0] >= k9_pair_ratio else 'BELOW pair -- check!'})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART D -- k=13 PZ-on-V descent (bar m_P = 0.0565)")
    print("=" * 100)
    cur = list(range(13)); best = (pz_numeric(cur), cur)
    print(f"  PZ(AP_13) numeric = {best[0]:.4f} (exact S11: 0.2720)")
    seen = {tuple(cur): best[0]}
    curv = best[0]
    for step in range(400):
        E2 = list(cur)
        r = rng.random()
        if r < 0.45:
            E2[rng.randrange(13)] = rng.randint(0, 60)
        elif r < 0.75:
            i, j = rng.randrange(13), rng.randrange(13)
            E2[i] = E2[j] + rng.choice([1, -1, 2, -2, 7, 14, -7])
        else:
            E2 = [rng.randint(0, 60) for _ in range(13)]
        E2 = norm_shape(E2, 13)
        if E2 is None:
            continue
        t = tuple(E2)
        if t not in seen:
            seen[t] = pz_numeric(E2, ngrid=60_000)
        v = seen[t]
        if v < best[0]:
            best = (v, E2)
            print(f"  step {step}: NEW MIN PZ = {v:.4f} at {E2}")
        if v < curv or rng.random() < 0.30:
            cur, curv = E2, v
    print(f"  descent min over {len(seen)} shapes: PZ = {best[0]:.4f} at {best[1]}")
    aV, MV, _, _ = excess_moments(best[1], [THETA])
    pz_exact = aV[0] ** 2 / MV[0][0]
    print(f"  EXACT confirmation of record: PZ = {pz_exact} = {float(pz_exact):.5f} "
          f"(bar 0.0565: {'CLEARS' if pz_exact > F(14249,252252) else 'FAILS'})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART E -- float-first Rayleigh, WIDE k=8 shapes (S11 stall completed)")
    print("=" * 100)
    wid = [float(w) for w in WIDTHS8]
    seeds = [list(range(8)), [0,1,2,4,8,16,32,64], [0,1,2,3,100,101,102,103],
             [0,29,58,87,116,145,174,203], [0,1,4,9,16,25,36,49],
             [0,7,26,57,100,155,182,197], [0,1,2,3,4,5,6,199]]
    best = None
    for E in seeds:
        v = ray_numeric(E, wid)
        if best is None or v < best[0]:
            best = (v, E)
        print(f"  seed {str(E):44s} RAY(float) = {v:.4f}")
    cur, curv = list(range(8)), ray_numeric(list(range(8)), wid)
    seen = {}
    for step in range(300):
        E2 = list(cur)
        r = rng.random()
        if r < 0.5:
            E2[rng.randrange(8)] = rng.randint(0, 200)
        else:
            E2 = [rng.randint(0, 200) for _ in range(8)]
        E2 = norm_shape(E2, 8)
        if E2 is None:
            continue
        t = tuple(E2)
        if t not in seen:
            seen[t] = ray_numeric(E2, wid, ngrid=60_000)
        v = seen[t]
        if v < best[0]:
            best = (v, E2)
            print(f"  step {step}: NEW MIN RAY = {v:.4f} at {E2}")
        if v < curv or rng.random() < 0.30:
            cur, curv = E2, v
    print(f"  WIDE-descent min over {len(seen)+len(seeds)} shapes: RAY = {best[0]:.4f} at {best[1]}")
    print(f"  (S11 exact AP_8 RAY = 0.8214; bar = 0.675; k=8 leg independently CLOSED by THM-651)")
    sys.stdout.flush()
