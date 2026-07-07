"""
lrc_2anchor_handoffs_opus_S139.py

THE FOUR 2-ANCHOR HANDOFFS (owner; continuing opus-S138 / boxeph HYP-4801 / klein-S160):

 (A1) EXACT NEAR-LIMIT SWEEP: exact PA_2 for spread APs {a+dj} over a (a,d) grid at k=8
      (+ spot k=9,10) and one-entry perturbations of the best — pins the finite correction
      below the T^2 limit (S138 saw (2,5): 0.79524 < limit 0.80169).
 (A2) CERTIFIED 2-ANCHOR LIMIT TABLE: E_y[meas(S(y) ∪ (S(y)-1/2))] on a y-grid with a
      rigorous Lipschitz error bound.  Each config point frac(jy) moves at speed j <= K, so
      each gap endpoint at speed <= K and |d/dy meas(S ∪ S-1/2)| <= 4K(K+1) crudely (each of
      <= 2(K+1) interval endpoints moving at <= 2K).  Grid step h => |table - truth| <=
      C_K h / 2 with C_K = 4K(K+1).  Reported as certified intervals.
 (A3) THE DOUBLE-COVER RIGIDITY CHAIN (rigorous, klein-S160's identity direction):
      doubling z -> 2z maps config(E,x) onto config(E,2x) 2-to-1; a point within delta of 0
      OR 1/2 downstairs lands within 2*delta of 0 upstairs, so
          gap@0(E, 2x) <= 2 * min(gap@0, gap@1/2)(E, x).
      Hence P(min(gap@0,gap@1/2) > 1/7) >= P(gap@0(E, 2x) > 2/7) = P(gap@0 > 2/7)  (x->2x
      measure-preserving), and PA_2 >= P(min > 1/7) >= the 2/7-tail of the ORIGIN gap.
      In the spread-AP limit the RHS = E_y[Sigma_{g>2/7} g] — EXACT via the S138 long-gap
      mass at theta = 2/7.  Verified pointwise + tabulated (honest: lossy; how lossy?).
 (A4) OBSERVER-RANK LEDGER CLASSIFIER: E[PS_0]-rank of the observer vs "shape in the proved
      ledger zone" over structured + random families — contingency table (cheap classifier?).
"""
from fractions import Fraction as F
from math import gcd
import random, sys, time, bisect

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import order_cells, cell_gap_affines
from lrc_PA2_shifted_roof_opus_S138 import PA2_exact_finite, EL_exact
from lrc_personal_space_tournament_opus_S138 import EPS_exact

THETA = F(1, 7)
TK = {8: 0.6185, 9: 0.5057, 10: 0.3956, 11: 0.2747, 12: 0.1429, 13: 0.0565}
LIMIT2 = {8: 0.80169, 9: 0.69274, 10: 0.59742, 11: 0.49776, 12: 0.42866, 13: 0.36647}

def PA2_exact_family(E, theta=THETA):
    """Exact PA_2 for an arbitrary integer family via order cells (anchored gaps at 0, 1/2)."""
    tot = F(0)
    for c0, c1 in order_cells(E):
        gaps = cell_gap_affines(E, c0, c1)
        subbp = set([c0, c1])
        for i in range(len(gaps)):
            ci, bi = gaps[i]
            for j in range(i + 1, len(gaps)):
                cj, bj = gaps[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if c0 < xc < c1: subbp.add(xc)
        for u, v in zip(sorted(subbp), sorted(subbp)[1:]):
            mid = (u + v) / 2
            pts = sorted(((e * mid) % 1, e) for e in E)
            phis = [p for p, _ in pts]
            def gform(anchor):
                i = bisect.bisect_left(phis, anchor)
                lo_e = pts[i - 1][1] if i > 0 else pts[-1][1]
                hi_e = pts[i][1] if i < len(pts) else pts[0][1]
                c = F(hi_e - lo_e)
                val = (((hi_e * mid) % 1) - ((lo_e * mid) % 1)) % 1
                return c, val - c * mid
            forms = [gform(F(1, 10**9)), gform(F(1, 2))]
            pts2 = set([u, v])
            for cc, bb in forms:
                if cc != 0:
                    xc = (theta - bb) / cc
                    if u < xc < v: pts2.add(xc)
            for uu, vv in zip(sorted(pts2), sorted(pts2)[1:]):
                mm = (uu + vv) / 2
                if max(cc * mm + bb for cc, bb in forms) > theta:
                    tot += vv - uu
    return tot

def part_A1():
    print("=" * 100)
    print("(A1) EXACT NEAR-LIMIT SWEEP (spread APs + perturbations)")
    print("=" * 100)
    for k in (8, 9, 10):
        best = None
        cap = 24 if k == 8 else 16
        for d in range(2, cap):
            for a in range(1, d):
                if gcd(a, d) != 1: continue
                v = PA2_exact_finite(a, d, k)
                if best is None or v < best[0]: best = (v, (a, d))
        v, (a, d) = best
        print(f"  k={k}: min exact PA_2 over (a,d) grid = {float(v):.6f} at (a,d)=({a},{d})"
              f"   T^2 limit {LIMIT2[k]:.5f}   T_k={TK[k]:.4f}   margin {float(v)-TK[k]:+.4f}")
        if k == 8:
            E0 = [a + d * j for j in range(k)]
            pbest = None
            for idx in range(k):
                for eps in (-2, -1, 1, 2):
                    E = sorted(set(E0[:idx] + [E0[idx] + eps] + E0[idx+1:]))
                    if len(E) != k or min(E) < 1: continue
                    v2 = PA2_exact_family(E)
                    if pbest is None or v2 < pbest[0]: pbest = (v2, E)
            print(f"       one-entry perturbations of the min: best = {float(pbest[0]):.6f} at {pbest[1]}"
                  f"   ({'DEEPER dip' if pbest[0] < v else 'no deeper'})")

def part_A2():
    print()
    print("=" * 100)
    print("(A2) CERTIFIED 2-ANCHOR LIMIT TABLE (grid + Lipschitz certificate)")
    print("=" * 100)
    G = 400_000
    for k in range(8, 14):
        K = k - 1
        CK = 4 * K * (K + 1)
        tot = 0.0
        for i in range(G):
            y = (i + 0.5) / G
            pts = sorted((j * y) % 1.0 for j in range(0, K + 1))
            ivs = []
            n = len(pts)
            for s in range(n):
                lo = pts[s]
                hi = pts[s + 1] if s < n - 1 else pts[0] + 1
                if hi - lo > 1.0 / 7:
                    ivs.append((lo, hi))
            flat = []
            for lo, hi in ivs:
                for sh in (0.0, 0.5):
                    l2, h2 = (lo - sh) % 1.0, 0.0
                    h2 = l2 + (hi - lo)
                    if h2 <= 1.0: flat.append((l2, h2))
                    else: flat.append((l2, 1.0)); flat.append((0.0, h2 - 1.0))
            flat.sort()
            m = 0.0; cur = None; hi0 = None
            for lo, hi in flat:
                if cur is None: cur, hi0 = lo, hi
                elif lo > hi0:
                    m += hi0 - cur; cur, hi0 = lo, hi
                else: hi0 = max(hi0, hi)
            if cur is not None: m += hi0 - cur
            tot += m
        est = tot / G
        err = CK / (2 * G)
        print(f"  k={k:2d}: 2-anchor limit = {est:.6f} ± {err:.6f} (certified)   T_k={TK[k]:.4f}"
              f"   certified margin >= {est - err - TK[k]:+.5f}")

def part_A3():
    print()
    print("=" * 100)
    print("(A3) DOUBLE-COVER RIGIDITY CHAIN: PA_2 >= P(min>1/7) >= P(gap@0 > 2/7); exact 2/7 mass")
    print("=" * 100)
    # pointwise verification of gap@0(E,2x) <= 2*min(gap@0, gap@1/2)(E,x)
    rng = random.Random(139)
    bad = 0
    for _ in range(20000):
        k = rng.choice([8, 13])
        E = sorted(rng.sample(range(1, 60), k))
        x = rng.random()
        pts = sorted((e * x) % 1.0 for e in E)
        def anch_gap(pts_, a):
            i = bisect.bisect_left(pts_, a)
            lo = pts_[i - 1] if i > 0 else pts_[-1] - 1
            hi = pts_[i] if i < len(pts_) else pts_[0] + 1
            return hi - lo
        g0 = anch_gap(pts, 1e-9); gh = anch_gap(pts, 0.5)
        pts2 = sorted((e * 2 * x) % 1.0 for e in E)
        g0up = anch_gap(pts2, 1e-9)
        if g0up > 2 * min(g0, gh) + 1e-9: bad += 1
    print(f"  pointwise identity check (20000 random): violations = {bad}")
    print("  exact 2/7-origin-gap mass (spread-AP limit lower bound for P(min>1/7)):")
    for k in range(8, 14):
        q = EL_exact(k, theta=F(2, 7))
        print(f"    k={k:2d}: E_y[Sigma_(g>2/7) g] = {q} = {float(q):.6f}   vs T_k={TK[k]:.4f}"
              f"   {'CLEARS (chain suffices on the limit!)' if float(q) > TK[k] else 'lossy (chain alone insufficient)'}")

def part_A4():
    print()
    print("=" * 100)
    print("(A4) OBSERVER-RANK vs LEDGER ZONE (kps-S60 intersected bites: k=13 diam<=75)")
    print("=" * 100)
    rng = random.Random(1394)
    BITE13 = 75
    rows = {(0, True): 0, (0, False): 0, (1, True): 0, (1, False): 0, ("2+", True): 0, ("2+", False): 0}
    for trial in range(120):
        if trial % 3 == 0:
            E = sorted(rng.sample(range(1, 40), 13))          # compressed-ish
        elif trial % 3 == 1:
            E = sorted(rng.sample(range(1, 200), 13))         # spread
        else:
            d = rng.choice([2, 3]); a = rng.randrange(1, d + 6)
            E = sorted(set(a + d * j + rng.choice([0, 0, 1]) for j in range(13)))
            while len(E) < 13: E.append(max(E) + rng.randrange(1, 5)); E = sorted(set(E))
        V = [0] + E
        eps = EPS_exact(V)
        order = sorted(eps.items(), key=lambda kv: (-kv[1], kv[0]))
        rank = [v for v, _ in order].index(0)
        rk = rank if rank <= 1 else "2+"
        inzone = (max(E) - min(E)) <= BITE13
        rows[(rk, inzone)] = rows.get((rk, inzone), 0) + 1
    print("   rank \\ in-ledger-zone:   True   False")
    for rk in (0, 1, "2+"):
        print(f"     {str(rk):>3}                {rows[(rk, True)]:5d}  {rows[(rk, False)]:5d}")
    print("   (classifier value = how strongly rank 0 predicts in-zone / rank>=1 predicts spread)")

if __name__ == "__main__":
    t0 = time.time()
    part_A1()
    part_A2()
    part_A3()
    part_A4()
    print(f"\n[total {time.time()-t0:.0f}s]")
