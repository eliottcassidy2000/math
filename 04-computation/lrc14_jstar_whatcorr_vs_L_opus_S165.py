"""
lrc14_jstar_whatcorr_vs_L_opus_S165.py   (opus-2026-07-08-S165)

UNIFY the capstone routes: show kps-S90's W-hat partial-sum correction (HYP-5507) is BOUNDED for
small longest-AP and LARGE for near-AP -- i.e. the small-L / near-AP split of the good period j* is
EXACTLY the opus-S154 L^1-converges / L^1-diverges (resonance) dichotomy.

kps route: j* <= N  <=  S_N := sum_{j=1}^{N} W(j/Vmax) > 0, and S_N = N*(6/7)^k + Corr_N,
  W(x) = uncovered measure of the cluster = sum_i (gap_i - 1/7)_+  (>0 iff x good, i.e. maxgap>1/7),
  j=0 anchor: W(0)=6/7 (all phases at 0).  So a good j<=N exists once N*(6/7)^k > |Corr_N|.
Corr_N is the resonance correction = the SAME W-hat sum as the density-floor tail (S154/LEM-011).

CLAIM (this file): |Corr_N| / (N*(6/7)^k) is SMALL for small longest-AP (route wins early, j* small)
and LARGE for near-AP (resonant, route needs N ~ ceil(7(k-1)/6) or the embedded-AP argument).  This
is the S154 mechanism (L^1 resonance sum converges for dissociated sets, diverges for AP-structured).
"""
import sys, random
from math import gcd


def W_uncovered(E, Vmax, j):
    """uncovered measure (in units of Vmax) of {e_i j mod Vmax}: sum of (gap - Vmax/7)_+ ."""
    k = len(E)
    ph = sorted((e * j) % Vmax for e in E)
    thr = Vmax / 7.0
    tot = 0.0
    prev = ph[-1] - Vmax
    for p in ph:
        g = p - prev
        if g > thr:
            tot += g - thr
        prev = p
    return tot / Vmax  # normalize to [0,1) circle


def jstar(E, Vmax):
    for j in range(1, Vmax):
        if W_uncovered(E, Vmax, j) > 1e-12:
            return j
    return None


def longest_ap(E):
    S = set(E); E = sorted(E); best = 1
    for i, a in enumerate(E):
        for b in E[i + 1:]:
            d = b - a
            if a - d in S: continue
            L = 2; x = b + d
            while x in S: L += 1; x += d
            best = max(best, L)
    return best


def corr_ratio(E, Vmax, N):
    """|Corr_N| / (N*(6/7)^k)  where S_N = sum_{j=1}^N W = N*(6/7)^k + Corr_N."""
    k = len(E)
    base = (6.0 / 7.0) ** k
    SN = sum(W_uncovered(E, Vmax, j) for j in range(1, N + 1))
    corr = SN - N * base
    return abs(corr) / (N * base), SN


def main():
    print("=" * 96)
    print("W-hat CORRECTION vs longest-AP: small-L bounded (route wins early) / near-AP large (resonant)")
    print("  = the opus-S154 L^1-converges/diverges dichotomy, driving kps-S90's W-hat route (HYP-5507)")
    print("=" * 96)
    r = random.Random(165)
    for k in (11, 13):
        print(f"\n  k={k}:  (base (6/7)^k = {(6/7)**k:.5f})")
        # bucket by longest-AP: max j*, and |Corr|/(N base) at N=ceil(7(k-1)/6)
        Nk = -(-7 * (k - 1) // 6)
        buckets = {}
        cands = []
        for _ in range(6000):
            spread = r.randint(k, 45)
            E = sorted(set([0] + r.sample(range(1, spread), k - 2) + [spread]))
            if len(E) == k: cands.append(E)
        for d in range(1, 6):
            base = [i * d for i in range(k)]
            cands.append(base)
            for jd in range(1, k):
                for dl in (-1, 1):
                    E = sorted(set(x + (dl if idx == jd else 0) for idx, x in enumerate(base)))
                    if len(E) == k and min(E) == 0: cands.append(E)
        for E in cands:
            spread = max(E)
            for Vmax in range(spread + 1, (7 * spread) // 6 + 1):
                js = jstar(E, Vmax)
                if js is None:
                    continue
                L = longest_ap(E)
                cr, SN = corr_ratio(E, Vmax, min(Nk, Vmax - 1))
                b = buckets.setdefault(L, [0, 0.0, 0])  # [max j*, max corr-ratio, count]
                b[0] = max(b[0], js); b[1] = max(b[1], cr); b[2] += 1
        print(f"    longest-AP L | max j* | max |Corr|/(N base) at N={Nk} | (count)")
        for L in sorted(buckets):
            b = buckets[L]
            flag = "  <- small-L: bounded" if b[1] < 1.0 else ("  <- near-AP: LARGE" if L >= k - 3 else "")
            print(f"       L={L:2d}      |  {b[0]:3d}   |   {b[1]:8.3f}                | ({b[2]}){flag}")
    print()
    print("  READING: small L -> |Corr|/(N base) < 1 (S_N>0 at small N => j* small, kps route wins);")
    print("  near-AP (L~k) -> ratio >> 1 (resonant, S154 L^1 divergence) => route needs N~ceil(7(k-1)/6)")
    print("  or mac-mini's embedded-AP Dirichlet.  UNIFIES kps-S90 + mac-mini + opus-S154/S164.")


if __name__ == "__main__":
    main()
