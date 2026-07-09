"""
lrc14_riesz_push_below_one_opus_S174.py   (opus-2026-07-09-S174)

Can a PRINCIPLED / globally-optimized Riesz product break the ratio int(M*R)/int(R) below 1 for a hard
LOOSE 13-set?  (<1 => positive-measure lonely set => S loose => a constructive looseness certificate;
if it generalizes uniformly => inf L>0.)  S173 naive coordinate-descent stalled at 1.07 on D=speeds.
Here: scipy differential_evolution (global) over several dissociated / speed-based D of increasing size,
tracking the best ratio.  Objective ratio(a)=int(M*R)/int(R), R=prod(1+a_m cos2pi m tau)>=0 (|a_m|<=1).
"""
import numpy as np
from scipy.optimize import differential_evolution

NG = 60013                      # coarse grid for optimization (prime)
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14
NGF = 400013                    # fine grid for final verification
TAUF = (np.arange(NGF) + 0.5) / NGF


def Mmult(S, tau):
    M = np.zeros(len(tau))
    for v in S:
        d = np.abs(((v * tau + 0.5) % 1.0) - 0.5)
        M += (d <= H).astype(float)
    return M


def make_ratio(S, D, tau):
    M = Mmult(S, tau)
    cb = np.array([np.cos(2 * np.pi * m * tau) for m in D])   # |D| x NG

    def ratio(a):
        R = np.ones(len(tau))
        for i in range(len(D)):
            R = R * (1 + a[i] * cb[i])
        r = R.mean()
        if r < 1e-6:
            return 5.0
        return float((M * R).mean() / r)
    return ratio


def best_ratio(S, D, maxiter=60, seed=0):
    ratio = make_ratio(S, D, TAU)
    bounds = [(-0.999, 0.999)] * len(D)
    res = differential_evolution(ratio, bounds, maxiter=maxiter, popsize=15, tol=1e-7,
                                 seed=seed, polish=True, init='sobol')
    # verify best on fine grid
    ratiof = make_ratio(S, D, TAUF)
    return ratiof(res.x), res.x


print("=" * 92)
print("RIESZ push below 1: global-opt ratio int(M*R)/int(R) for hard LOOSE sets.  main = 13/7 = 1.857")
print("=" * 92)
sets = {
    "{1..13}\\{6} U {56} [LOOSE .0056]": sorted(set(range(1, 14)) - {6} | {56}),
    "{1..12} U {182}     [LOOSE .024] ": list(range(1, 13)) + [182],
    "7*{1..12} U {13}    [LOOSE .029] ": sorted(set([7 * i for i in range(1, 13)] + [13])),
    "{1..13}             [TIGHT lm=0] ": list(range(1, 14)),                 # validity anchor
    "2*{1..13}           [TIGHT lm=0] ": [2 * i for i in range(1, 14)],      # dilate, still tight
}
for name, S in sets.items():
    print(f"\n  S = {name}")
    Ds = {
        "speeds13": sorted(S),
        "speeds+lac": sorted(S) + [128, 256, 512],
        "lac-dense": [1, 2, 3, 5, 8, 13, 21, 34, 55, 89],       # Fibonacci-ish (near-dissociated)
    }
    best_overall = 9.9
    for label, D in Ds.items():
        r, a = best_ratio(S, D, maxiter=50, seed=1)
        best_overall = min(best_overall, r)
        flag = "  <<< BELOW 1 (CERTIFIES LOOSE)" if r < 1.0 - 1e-4 else ""
        print(f"    D={label:>12} (|D|={len(D):2d}): best ratio = {r:.4f}{flag}")
    print(f"    => best over D: {best_overall:.4f}  ({'CERTIFICATE FIRES' if best_overall < 1 else 'still >= 1'})")
print()
print("  READING: ratio<1 on a LOOSE set = the Riesz method firing (constructive looseness). The TIGHT")
print("  set must stay >=1 (validity). If global-opt also stalls >1, the <1 bar genuinely needs the")
print("  structured dissociated+hypercontractive construction (Bedert 2025), not generic optimization.")
print("=" * 92)
