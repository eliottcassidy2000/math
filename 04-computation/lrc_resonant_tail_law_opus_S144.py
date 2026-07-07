"""
lrc_resonant_tail_law_opus_S144.py   (opus-2026-07-07-S144, HYP-5137, part 3)

THE RESONANT-TAIL LAW (new conjecture, born from the part-2 V_r data):

  V_r(E) = meas{x : exactly r circular gaps of {frac(ex)} exceed 1/7}   (affine-invariant)
  T_j(E) = sum_{r >= j} V_r(E)                                          (tail sums)

  (P1)  Among k-element integer families, the AP maximizes T_j for j >= 4
        (high-resonance mass), while minimizing the total mu = T_1 (R2/(A)).
        The mu-minimality is the trade-off: perturbations destroy resonant mass
        faster than... no -- they CREATE bulk (V_1..V_3) mass faster than they
        destroy resonant (V_4..V_6) mass.  P1 is the destruction half.

(1) EXACT V_r(AP_13) from the roof + cluster rates.  Near node p/q (q <= 6) the phase
    vector is EXACTLY linear in delta = x - p/q: cluster residue-r position = r*inv + e*delta.
    The q inter-cluster gaps close at integer rates (min of next cluster - max of current,
    per side); gap i survives > 1/7 until delta_i = (1/q - 1/7)/rate_i; sorting the rates
    slices the window into r-constant pieces.  Intra-cluster sub-gaps stay < 1/7 inside
    all windows for the AP (max spread 2*delta < 1/7).  Sum with Fractions; verify vs
    numerics; report exact T_4, T_5, T_6 -- the law's RHS.

(2) ADVERSARIAL SEARCH on P1: structured zoo + random + hill-climb from AP and from
    random starts, objective = maximize T_5 (and T_4, T_6), k = 13 and k = 8.
    (MISTAKE-119 discipline: structured families + multi-start local search, not just
    random sampling.)
"""
from fractions import Fraction as F
from math import gcd
import numpy as np
import itertools, random, time

THETA = F(1, 7)

# ---------------------------------------------------------------- (1) exact AP V_r
def exact_Vr_AP(k):
    """Exact V_r profile of AP_k = {1..k} at theta = 1/7, for k >= 12 (disjoint windows,
       one q<=6 node per Farey cell boundary).  Returns dict r -> Fraction."""
    V = {r: F(0) for r in range(1, 7)}
    # q = 1 (origin + 1): wrap gap 1 - (k-1)delta > 1/7 while delta < (6/7)/(k-1); r = 1.
    V[1] += 2 * (F(6, 7) / (k - 1))
    E = list(range(1, k + 1))
    for q in range(2, 7):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            # clusters: residue class c (position c/q) holds elements e with e*p == c (mod q)
            clusters = {}
            for e in E:
                clusters.setdefault((e * p) % q, []).append(e)
            order = sorted(clusters)  # positions c/q ascending; all q residues occupied for AP
            assert len(order) == q
            for side in (+1, -1):
                rates = []
                for i, c in enumerate(order):
                    cn = order[(i + 1) % len(order)]
                    if side == +1:
                        rate = max(clusters[c]) - min(clusters[cn])   # gap shrinks at this rate
                    else:
                        rate = max(clusters[cn]) - min(clusters[c])
                    assert rate > 0, (q, p, side, c)
                    rates.append(rate)
                rates.sort()
                Wd = F(1, q) - THETA
                # gap with rate rates[i] survives until Wd/rates[i]; slice
                deltas = sorted({Wd / r for r in rates}, reverse=True)
                # count of gaps surviving in (d_{next}, d): #rates with Wd/rate >= d
                prev = F(0)
                for d in sorted(deltas):
                    cnt = sum(1 for r in rates if Wd / r >= d)
                    # on (prev, d): gaps surviving = #{i: Wd/rates[i] > prev} evaluated inside
                    pass
                # simpler: walk pieces between consecutive breakpoints ascending
                bps = [F(0)] + sorted(set(Wd / r for r in rates))
                for i in range(len(bps) - 1):
                    lo, hi = bps[i], bps[i + 1]
                    mid_cnt = sum(1 for r in rates if Wd / r > lo)  # survive past lo
                    V[mid_cnt] += (hi - lo)
    return V

# ---------------------------------------------------------------- numeric V_r
def Vr_numeric(E, res=800_000, theta=1 / 7, chunk=200_000):
    E = np.asarray(sorted(E), dtype=np.float64)
    Vc = np.zeros(8, dtype=np.int64)
    for lo in range(0, res, chunk):
        x = (np.arange(lo, min(lo + chunk, res)) + 0.5) / res
        ph = np.sort((x[:, None] * E[None, :]) % 1.0, axis=1)
        gaps = np.diff(ph, axis=1)
        wrap = (ph[:, 0] + 1 - ph[:, -1])[:, None]
        allg = np.concatenate([gaps, wrap], axis=1)
        rcount = (allg > theta).sum(axis=1)
        for r in range(1, 7):
            Vc[r] += int((rcount == r).sum())
    return {r: Vc[r] / res for r in range(1, 7)}

def tails(V):
    return {j: sum(V[r] for r in range(j, 7)) for j in range(1, 7)}

# ---------------------------------------------------------------- (2) adversarial search
def stress(k, TAP, nrand=200, seed=1):
    """Return worst (max) T_4, T_5, T_6 found over zoo + random + hill-climbs."""
    rng = random.Random(seed)
    zoo = {
        "AP": list(range(1, k + 1)),
        "GW-like": (list(range(1, k - 1)) + [k, 2 * (k - 1)])[:k],
        "AP+tail": list(range(1, k)) + [k + 1],
        "two-block": list(range(1, k // 2 + 1)) + list(range(3 * k, 3 * k + (k - k // 2))),
        "near-AP swap": list(range(1, k - 1)) + [k - 1 + 2, k + 2],
        "mod6-perfect-wide": [6 * j + (j % 6) for j in range(1, k + 1)],
        "mod7-perfect": [7 * j + (j % 7) for j in range(1, k + 1)],
        "tight+far": list(range(1, k)) + [14 * (k - 1)],
    }
    for t in range(nrand):
        d = rng.randint(20, 300)
        zoo[f"rnd{t}"] = sorted(rng.sample(range(1, d), k - 1)) + [d]
    best = {4: (TAP[4], "AP"), 5: (TAP[5], "AP"), 6: (TAP[6], "AP")}
    for nm, E in zoo.items():
        E = sorted(set(E))
        if len(E) != k:
            continue
        V = Vr_numeric(E, res=300_000)
        T = tails(V)
        for j in (4, 5, 6):
            if T[j] > best[j][0] + 1e-9:
                best[j] = (T[j], nm if nm.startswith(("AP", "GW")) else f"{nm}:{E}")
    # hill-climb on T_5 from AP and 3 random starts
    def T5(E):
        return tails(Vr_numeric(E, res=200_000))[5]
    for start in range(4):
        if start == 0:
            E = list(range(1, k + 1))
        else:
            E = sorted(rng.sample(range(1, 60), k))
        cur = T5(E)
        for it in range(60):
            improved = False
            for i in range(k):
                for step in (-2, -1, 1, 2, 5, -5):
                    E2 = sorted(set(E) - {E[i]} | {max(1, E[i] + step)})
                    if len(E2) != k:
                        continue
                    v = T5(E2)
                    if v > cur + 1e-6:
                        E, cur, improved = E2, v, True
                        break
                if improved:
                    break
            if not improved:
                break
        if cur > best[5][0] + 2e-3:
            best[5] = (cur, f"hillclimb:{E}")
    return best

def main():
    t0 = time.time()
    print("=" * 100)
    print("(1) EXACT V_r(AP_13) from roof + cluster rates, vs numerics")
    print("=" * 100)
    Vx = exact_Vr_AP(13)
    Vn = Vr_numeric(list(range(1, 14)), res=2_000_000)
    print("    r :   exact              float      numeric     |diff|")
    for r in range(1, 7):
        print(f"    {r} :   {str(Vx[r]):>14}   {float(Vx[r]):.6f}   {Vn[r]:.6f}   {abs(float(Vx[r])-Vn[r]):.6f}")
    s = sum(Vx.values())
    print(f"    sum = {s} = {float(s):.6f}  vs mu(AP_13) = 477/1078 = {477/1078:.6f}"
          f"   {'MATCH' if s == F(477,1078) else '*** MISMATCH ***'}")
    Tx = {j: sum(Vx[r] for r in range(j, 7)) for j in range(1, 7)}
    print("    EXACT TAILS: " + "  ".join(f"T_{j} = {Tx[j]} = {float(Tx[j]):.5f}" for j in (4, 5, 6)))

    print()
    print("=" * 100)
    print("(2) ADVERSARIAL SEARCH: can anything beat the AP's resonant tail?  k=13")
    print("=" * 100)
    TAP13 = tails(Vr_numeric(list(range(1, 14)), res=800_000))
    best13 = stress(13, TAP13)
    for j in (4, 5, 6):
        val, who = best13[j]
        print(f"    max T_{j}: {val:.5f} at {who}   (AP: {TAP13[j]:.5f})"
              f"{'   AP-MAX SURVIVES' if who == 'AP' else '   *** BEATEN ***'}")

    print()
    print("    k=8:")
    TAP8 = tails(Vr_numeric(list(range(1, 9)), res=800_000))
    best8 = stress(8, TAP8, seed=2)
    for j in (4, 5, 6):
        val, who = best8[j]
        print(f"    max T_{j}: {val:.5f} at {who}   (AP: {TAP8[j]:.5f})"
              f"{'   AP-MAX SURVIVES' if who == 'AP' else '   *** BEATEN ***'}")

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
