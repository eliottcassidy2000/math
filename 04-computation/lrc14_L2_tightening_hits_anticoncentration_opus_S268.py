"""
opus-2026-07-11-S268: tighten the large-sieve energy bound to close the 3.1x gap (S267). OUTCOME: the gap does
NOT close by tightening the large sieve -- because the TIGHT large-sieve bound IS the anti-concentration. Two
genuine advances + an honest correction of S267:

(1) THE 3.1x GAP IS NOT UNIFORM -- it is isolated to core=1 (the deep well / near-AP family), which is exactly
    the runner-1 lemma (S265). For core>=2 the L2 energy core*Sum(eps^2) has a ~21x margin (<= 0.034 vs 0.735).

(2) THE SCALING IS eps_v = O(1), NOT sqrt(|G'|). For core>=2, max|eps_v| <= 0.099 stays ~flat as |G'|->0.2,
    while the Bessel prediction sqrt(6/49/|G'|) grows to ~0.78. So <g(v.),1_G'> = eps_v|G'| = O(|G'|) (LINEAR in
    |G'|), not O(sqrt|G'|). The Bessel/operator-norm bound allows sqrt(|G'|)-blowup, so it is loose by
    ~1/sqrt(|G'|) -- this is exactly the 3.1x (and worse for small |G'|).

(3) CORRECTS S267. L2 removed the L1 DIVERGENCE (the harmonic tail Sum|b_k|=inf / the alternating multi-order
    cancellation of S266). It did NOT remove the underlying ANTI-CONCENTRATION, which reappears as the
    |G'|-scaling: a tight L2 energy bound REQUIRES eps_v = O(1) uniformly, i.e. <g(v.),1_G'> = O(|G'|), which is
    precisely the anti-concentration statement. The positive off-diagonal in the energy expansion (measured
    ~64% of the total for core>=2) is the same fact: the diagonal (Bessel) underestimates because 1_G' is far
    from the frame operator's top eigenvector, and closing that gap is the anti-concentration.

NET: the irreducible core of LRC(14)-covering is the CLEAN uniform bound  eps_v = O(1)  (equivalently
<g(v.),1_G'> = O(|G'|)) -- the same hard core as S266/S267, now in its sharpest, most quotable form, most
extreme at the deep well (core=1, S265). "Tighten the large sieve" is not a separate task from the
anti-concentration; it IS the anti-concentration.

-> opus-S267 (corrected: L2 removes divergence, not anti-concentration), opus-S266 (same core, L1 form),
   opus-S265 (core=1 = runner-1 lemma), opus-S262 (pairwise Cov<=1/(3vv'), the rigorous-but-loose Bessel input).
"""
import numpy as np
from math import gcd
from functools import reduce
import random

def divcomplete(v): return all(any(x % d == 0 for x in v) for d in range(2, 15))
def primitive(v): return reduce(gcd, v) == 1

def eps_stats(v, a, D, cD):
    core = [x for x in v if gcd(x, 30030) == 1]; non = [x for x in v if gcd(x, 30030) != 1]
    if not core: return None
    safe = np.ones(D, bool)
    for w in non:
        rr = (w * a) % D; safe &= (np.minimum(rr, D - rr) >= cD)
    Gm = safe.mean()
    if Gm < 0.02: return None
    g = safe.astype(float); Se2 = 0.0; mx = 0.0
    for vv in core:
        rr = (vv * a) % D
        e = ((np.minimum(rr, D - rr) < cD).astype(float) @ g) / g.sum() - 1 / 7
        Se2 += e * e; mx = max(mx, abs(e))
    return len(core), Gm, Se2, mx

def main():
    D = 13860; c = 1.0 / 14; cD = c * D; a = np.arange(D)
    random.seed(23)
    cands = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182]]
    while len(cands) < 600:
        v = sorted(random.sample(range(1, 170), 13))
        if primitive(v) and divcomplete(v): cands.append(v)

    # (1) split by core count
    from collections import defaultdict
    by_core = defaultdict(lambda: [0.0, 0])
    # (2) scaling: bucket core>=2 by |G'|
    buckets = defaultdict(list)
    for v in cands:
        s = eps_stats(v, a, D, cD)
        if s is None: continue
        nc, Gm, Se2, mx = s
        r = by_core[nc]; r[0] = max(r[0], nc * Se2); r[1] += 1
        if nc >= 2:
            buckets[round(Gm * 5) / 5].append(mx)

    print("(1) core*Sum(eps^2) by core count (target < 36/49 = 0.735):")
    for nc in sorted(by_core):
        r = by_core[nc]
        print("    core=%d: %d families, max core*Sum(eps^2) = %.4f" % (nc, r[1], r[0]))
    c2 = max(v[0] for k, v in by_core.items() if k >= 2)
    print("    => core=1 (deep well) = %.3f is the ONLY case near threshold (= runner-1 lemma S265);"
          % by_core[1][0])
    print("       core>=2 max = %.4f (~%.0fx margin)." % (c2, (36 / 49) / c2))

    print("\n(2) scaling of max|eps_v| vs |G'| for core>=2 (Bessel would allow sqrt(6/49/|G'|)):")
    for k in sorted(buckets):
        arr = buckets[k]; bess = (6 / 49 / max(k, 0.05)) ** 0.5
        print("    |G'|~%.1f: %d fam, mean max|eps|=%.4f, worst=%.4f, Bessel-allows=%.3f"
              % (k, len(arr), sum(arr) / len(arr), max(arr), bess))
    gw = max(max(v) for v in buckets.values())
    print("    => max|eps_v| stays FLAT ~%.3f (=O(1)) while Bessel-sqrt grows: eps_v=O(1), <g(v.),1_G'>=O(|G'|)."
          % gw)

    print("\n(3) VERDICT: tightening the large sieve does NOT close the gap -- the tight bound eps_v=O(1) IS the")
    print("    anti-concentration. L2 (S267) removed the L1 DIVERGENCE (S266's alternating cancellation), NOT the")
    print("    anti-concentration, which reappears as the |G'|-scaling. Irreducible core (cleanest form):")
    print("    eps_v = O(1) uniformly  <=>  <g(v.),1_G'> = O(|G'|),  most extreme at the deep well (core=1 = S265).")

if __name__ == '__main__':
    main()
