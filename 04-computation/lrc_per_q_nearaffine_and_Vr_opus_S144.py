"""
lrc_per_q_nearaffine_and_Vr_opus_S144.py   (opus-2026-07-07-S144, HYP-5137, part 2)

TWO QUESTIONS left open by the part-1 audit (lrc_per_q_audit_opus_S144.py):

(1) NEAR-AFFINE ADVERSARIES vs the per-q (Voronoi) minimality claim of HYP-5117.
    The affine image 2AP-1 (all-odd) undercuts W_1, W_3, W_5 with EQUAL mu -- kps calls
    these rows bookkeeping artifacts since mu is affine-invariant and per-q is not.
    But a PRIMITIVE, NON-AFFINE small perturbation of the spread AP inherits the shifted
    per-q profile while being a legitimate comparison family.  If it undercuts any
    W_q(AP_13) (q = 2..6), the per-q minimality LEMMA AS STATED is refuted in the
    primitive normalization, and the program needs a dilation-covariant restatement.

(2) THE INTRINSIC DECOMPOSITION V_r(E) = meas{x : exactly r circular gaps > 1/7}.
    r(x) is a function of the phase MULTISET, so V_r is fully affine-invariant
    (dilation + translation + reflection) -- the dilation bookkeeping problem
    disappears by construction.  mu = sum_{r>=1} V_r.  Questions:
      (a) verify affine invariance numerically (spread AP == AP profile);
      (b) does the AP minimize each V_r?  or the tail sums sum_{r>=j} V_r
          (a majorization/fragmentation statement)?
      (c) closed forms for V_r(AP_13) from the roof windows (recorded in comments).
"""
from math import gcd
import numpy as np
import time

def profiles(E, res=2_000_000, theta=1 / 7, Q=6, chunk=200_000):
    """Returns mu, W[1..6] (nearest-rational Voronoi attribution), V[1..6] (#gaps>theta)."""
    E = np.asarray(sorted(E), dtype=np.float64)
    nodes = []
    for q in range(1, Q + 1):
        for p in range(0, q + 1):
            if gcd(p, q) == 1:
                nodes.append((p / q, q))
    nv = np.array([n[0] for n in nodes])
    nq = np.array([n[1] for n in nodes])
    Wc = np.zeros(Q + 1, dtype=np.int64)
    Vc = np.zeros(8, dtype=np.int64)
    good_total = 0
    for lo in range(0, res, chunk):
        x = (np.arange(lo, min(lo + chunk, res)) + 0.5) / res
        ph = np.sort((x[:, None] * E[None, :]) % 1.0, axis=1)
        gaps = np.diff(ph, axis=1)
        wrap = (ph[:, 0] + 1 - ph[:, -1])[:, None]
        allg = np.concatenate([gaps, wrap], axis=1)
        rcount = (allg > theta).sum(axis=1)
        good = rcount > 0
        good_total += int(good.sum())
        for r in range(1, 7):
            Vc[r] += int((rcount == r).sum())
        if good.any():
            xg = x[good]
            d = np.abs(xg[:, None] - nv[None, :])
            d = np.minimum(d, 1 - d)
            qq = nq[np.argmin(d, axis=1)]
            for q in range(1, Q + 1):
                Wc[q] += int((qq == q).sum())
    mu = good_total / res
    return (mu, {q: Wc[q] / res for q in range(1, Q + 1)},
            {r: Vc[r] / res for r in range(1, 7)})

def main():
    t0 = time.time()
    k = 13
    AP = list(range(1, k + 1))
    muAP, WAP, VAP = profiles(AP)
    spread = [2 * j - 1 for j in range(1, k + 1)]          # 2AP - 1 (affine image)

    fams = {
        "AP {1..13}": AP,
        "spread 2AP-1 {1,3,..,25}": spread,
        "spread3 3AP-2 {1,4,..,37}": [3 * j - 2 for j in range(1, k + 1)],
        # primitive near-affine perturbations of the spread AP (NOT affine images):
        "spread+bump {1,3,..,23,26}": [2 * j - 1 for j in range(1, k)] + [26],
        "spread+bump2 {1,3,..,23,27}": [2 * j - 1 for j in range(1, k)] + [27],
        "spread+mid {1,3,..,11,14,15,..,25odd}": sorted(set(spread) - {13} | {14}),
        "spread+head {2,3,5,7,..,25}": sorted(set(spread) - {1} | {2}),
        "spread3+bump {1,4,..,34,38}": [3 * j - 2 for j in range(1, k)] + [38],
    }
    print("=" * 108)
    print("(1) PER-q VORONOI PROFILE W_q  --  near-affine primitive adversaries vs AP (k=13, res 2e6)")
    print("=" * 108)
    print(f"    W_q(AP_13) exact: 1/7, 5/77, 3/35, 8/147, 23/294, 4/245")
    print(f"    {'family':38s} {'mu':>7} {'dmu':>8}  {'W_1..W_6':>44}")
    rows = {}
    for nm, E in fams.items():
        g = 0
        diffs = sorted(E)
        for i in range(len(diffs) - 1):
            g = gcd(g, diffs[i + 1] - diffs[i])
        mu, W, V = profiles(E)
        rows[nm] = (mu, W, V)
        flags = [f"q{q}<AP" for q in range(1, 7) if W[q] < WAP[q] - 0.003]
        wstr = " ".join(f"{W[q]:.4f}" for q in range(1, 7))
        print(f"    {nm:38s} {mu:7.4f} {mu-muAP:+8.4f}  {wstr}  gcd-diff={g}"
              f"{('  *** BELOW AP at: ' + ','.join(flags) + ' ***') if flags else ''}")

    print()
    print("=" * 108)
    print("(2) INTRINSIC PROFILE V_r = meas{x : exactly r gaps > 1/7}  (affine-invariant by construction)")
    print("=" * 108)
    print(f"    {'family':38s} {'mu':>7}  {'V_1..V_6':>44}   tail sums T_j = sum_(r>=j) V_r")
    for nm, (mu, W, V) in rows.items():
        vstr = " ".join(f"{V[r]:.4f}" for r in range(1, 7))
        tails = [sum(V[r] for r in range(j, 7)) for j in range(1, 7)]
        tailsAP = [sum(VAP[r] for r in range(j, 7)) for j in range(1, 7)]
        flags = [f"V{r}<AP" for r in range(1, 7) if V[r] < VAP[r] - 0.003]
        tflags = [f"T{j}<AP" for j in range(1, 7) if tails[j - 1] < tailsAP[j - 1] - 0.003]
        print(f"    {nm:38s} {mu:7.4f}  {vstr}"
              f"{('   [' + ','.join(flags) + ']') if flags else '   [per-r >= AP]'}"
              f"{('  [' + ','.join(tflags) + ' TAILS]') if tflags else '  [tails >= AP]'}")

    print()
    print("    V_r(AP_13) roof-derived structure (for the record):")
    print("      q=1 window: r=1 throughout       -> V_1 gains 1/7")
    print("      q=2 window: r=2 throughout       -> V_2 gains 5/77")
    print("      q=3 windows: r=3 then r=1 slivers; q=4: r=4->2->1...; etc.")
    print(f"    measured V_r(AP_13): " + " ".join(f"{VAP[r]:.4f}" for r in range(1, 7)))

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
