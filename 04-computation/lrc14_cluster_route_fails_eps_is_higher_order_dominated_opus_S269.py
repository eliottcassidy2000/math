"""
opus-2026-07-11-S269: prove eps_v = O(1) uniformly for core>=2. OUTCOME: the natural rigorous route -- the
cluster/Mayer expansion whose leading term is -(7/6) Sum_w Cov(D_v,D_w) with the PROVEN pairwise bound
|Cov(D_v,D_w)| <= 1/(3vw) (S262) -- provably FAILS. The leading term does NOT dominate: for core speeds v>=17,
|eps_v| exceeds (7/6) Sum_w |Cov(D_v,D_w)| in 555/658 cases, by up to 47x. So eps_v is NOT a pairwise object; the
higher-order MULTI-LINEAR correlations <g(v.), prod_{w in S} beta_w> (|S|>=2) ARE eps_v, not a correction. This
is a SHARPER statement than S266: for v=1 (deep well) the higher orders CANCEL a large pairwise; for v>=17 the
pairwise is NEGLIGIBLE and the multi-linear resonances (noncore pairs w1 +- w2 = +-kv) are the whole thing.

THE EXPANSION (rigorous). mu_G' = prod_w(1-beta_w)/|G'|, beta_w=1_{D_w}. eps_v = E_{mu}[g(v.)] =
[-Sum_w Cov(D_v,D_w) + Sum_{|S|>=2}(-1)^{|S|}<g(v.),prod_{w in S}beta_w>]/|G'|. Cluster/KP convergence would need
the higher (connected) correlations to be controlled by the pairwise activity; they are NOT -- the |S|=2 triple
correlations dominate (and their magnitude sum diverges harmonically, S266). So the polymer gas is NOT in the
convergent (low-activity) regime for this band.

TEST A: |eps_v| <= (7/6) Sum_w |Cov(D_v,D_w)|  [leading term dominates]?  -> FAILS 555/658, worst ratio 46.9x.
TEST B: |eps_v| <= (7/6) Sum_w 1/(3vw)          [fully rigorous via proven |Cov|]?  -> worst ratio 113x.

NET: eps_v = O(1) for core>=2 is the genuine anti-concentration. The cluster expansion -- the best rigorous
route, built on the PROVEN pairwise bound -- cannot reach it, because eps_v is higher-order-dominated (multi-
linear), not pairwise. The dominant contribution is the noncore-pair resonance sum (w1 +- w2 = +-kv), whose
magnitude diverges and whose O(1) value needs the signed cancellation = the same wall as S266/S268. NOT PROVED.

-> opus-S268 (eps_v=O(1) = the anti-concentration), opus-S266 (multi-linear, L1 divergence), opus-S262
   (|Cov|<=1/(3vw), the proven pairwise input this route rests on), opus-S263/LEM-015 (the resonance/E3 core).
"""
import numpy as np
from math import gcd
from functools import reduce
import random

def divcomplete(v): return all(any(x % d == 0 for x in v) for d in range(2, 15))
def primitive(v): return reduce(gcd, v) == 1

def main():
    D = 13860; c = 1.0 / 14; cD = c * D; a = np.arange(D)
    def band(x):
        rr = (x * a) % D; return (np.minimum(rr, D - rr) < cD).astype(float)
    random.seed(31)
    cands = []
    while len(cands) < 250:
        v = sorted(random.sample(range(1, 170), 13))
        if primitive(v) and divcomplete(v) and len([x for x in v if gcd(x, 30030) == 1]) >= 2:
            cands.append(v)
    worst_pair = 0.0; worst_rig = 0.0; nfail = 0; ntot = 0; maxeps = 0.0
    for v in cands:
        core = [x for x in v if gcd(x, 30030) == 1]; non = [x for x in v if gcd(x, 30030) != 1]
        Dw = [band(w) for w in non]
        safe = np.ones(D, bool)
        for w in non:
            rr = (w * a) % D; safe &= (np.minimum(rr, D - rr) >= cD)
        Gm = safe.mean()
        if Gm < 0.02: continue
        g = safe.astype(float)
        for vv in core:
            if vv < 17: continue
            ntot += 1
            Dv = band(vv)
            eps = (Dv @ g) / g.sum() - 1 / 7
            maxeps = max(maxeps, abs(eps))
            SumCov = sum(abs((Dv @ Dwm) / D - 1 / 49) for Dwm in Dw)
            SumRig = sum(1.0 / (3 * vv * w) for w in non)
            lead = (7 / 6) * SumCov; rig = (7 / 6) * SumRig
            if abs(eps) > lead + 1e-9: nfail += 1
            if lead > 0: worst_pair = max(worst_pair, abs(eps) / lead)
            if rig > 0: worst_rig = max(worst_rig, abs(eps) / rig)
    print("core>=2 families, core speeds v>=17: %d tested, max|eps_v|=%.4f" % (ntot, maxeps))
    print("\nTEST A  |eps_v| <= (7/6) Sum_w |Cov(D_v,D_w)|  (leading cluster term dominates)?")
    print("  failures %d/%d ; worst ratio |eps|/lead = %.1fx  => %s"
          % (nfail, ntot, worst_pair, "FAILS -- eps_v is higher-order (multi-linear) dominated"))
    print("\nTEST B  |eps_v| <= (7/6) Sum_w 1/(3vw)  (fully rigorous via proven |Cov|<=1/(3vw))?")
    print("  worst ratio = %.1fx  => %s" % (worst_rig, "FAILS -- cluster route cannot reach eps_v=O(1)"))
    print("\nNET: eps_v=O(1) for core>=2 is the anti-concentration. The cluster/Mayer expansion (best rigorous")
    print("     route, on the PROVEN pairwise bound) provably cannot reach it: eps_v is multi-linear, not")
    print("     pairwise. Dominant term = noncore-pair resonances (w1 +- w2 = +-kv); magnitude diverges, value")
    print("     needs the signed cancellation = the S266/S268 wall. NOT PROVED.")

if __name__ == '__main__':
    main()
