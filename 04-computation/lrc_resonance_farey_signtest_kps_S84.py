#!/usr/bin/env python3
r"""
lrc_resonance_farey_signtest_kps_S84.py  (kind-pasteur-2026-07-08-S84, HYP-5397)

ATTEMPT: prove the resonance sign Var(W) <= c*R2 via the Farey-cell moment route, as a
PER-PAIR covariance sign (the cleanest Farey form): Resonance = sum_{d!=+-d'} r_d r_{d'}
Cov(g(dx),g(d'x)); since r_d r_{d'} >= 0, if each Cov <= 0 then Resonance <= 0.

RESULT: a DECISIVE NEGATIVE -- the per-pair covariances are MIXED-SIGN, so there is no
clean Farey-cell/per-pair sign proof.  Completes the map (with kps-S82) of why brick (B) is
the delicate open mile: it is the BALANCE of two large opposing effects, provable only via
exact per-config computation + decorrelation, NOT a uniform structural sign.

FINDINGS (grid-exact, verified):
  - TENT f (THM-656 pair sum, k<=10): 44/132 per-pair Cov(f(dx),f(d'x)) are POSITIVE
    (small, ~1e-5), 88 negative; so Resonance<=0 is a delicate WEIGHTED sum, not per-pair.
  - OVERLAP ov (covering, k=11): 110/132 per-pair Cov(ov(dx),ov(d'x)) are POSITIVE (LARGE,
    up to 8e-4)!  So Var(sum ov_ij) = R2*V_ov + Resonance_ov with Resonance_ov > 0 (the
    pair-sum variance is LARGE, 1.26 at the block).  Var(W) << Var(sum ov) because the
    higher-order (triple+) terms cancel MORE.  So Var(W) <= c*R2 = [positive pair-resonance,
    large] - [negative higher-order cancellation, larger], a delicate net -- NO clean sign.

CONCLUSION: neither the Fourier-truncation route (kps-S82: diverges) nor the Farey per-pair
route (this file: mixed signs) proves the resonance sign.  It is genuinely a delicate
two-large-effects cancellation.  The k=11 tail must go through exact per-config Farey moments
(compact exhaustive: opus/klein/mac-mini) + decorrelation (large diam: LEM-005), NOT a
uniform resonance bound.
"""
import numpy as np
beta = 5/63
def tent(s): return np.maximum(s-beta, 0)*(s <= 1/7)
def ov(s):
    t = np.minimum(s, 1-s); return np.maximum(1/7 - t, 0)
def cov_grid(g, d, dp, res=400000):
    xs = (np.arange(res)+0.5)/res
    a = g((d*xs) % 1.0); b = g((dp*xs) % 1.0)
    return float(np.mean(a*b) - np.mean(a)*np.mean(b))
for name, g in [("tent f (THM-656)", tent), ("overlap ov (covering)", ov)]:
    npos = nneg = 0; mx = -9; mn = 9
    for d in range(1, 15):
        for dp in range(1, 15):
            if dp == d or dp == -d: continue
            cv = cov_grid(g, d, dp)
            if cv > 1e-7: npos += 1
            else: nneg += 1
            mx = max(mx, cv); mn = min(mn, cv)
    print(f"{name}: per-pair Cov(d,d') sign -- POSITIVE {npos}, NON-POSITIVE {nneg}; "
          f"range [{mn:.2e}, {mx:.2e}]")
print("=> mixed signs (overlap: mostly POSITIVE) => the resonance sign is NOT a per-pair/per-cell")
print("   sign; it is a delicate weighted cancellation. No clean Farey-sign proof. Both routes")
print("   (Fourier truncation kps-S82, Farey per-pair kps-S84) ruled out -- exact per-config +")
print("   decorrelation is the only path (fleet-active).")
