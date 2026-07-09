"""
lrc14_riesz_product_certificate_opus_S173.py   (opus-2026-07-09-S173)

Push the RIESZ-PRODUCT certificate for L(S)>0 (S loose) -- THM-515 part C / HYP-2540, the Bedert-2025
method (arXiv:2511.16636).  Certificate: find a Riesz product R(tau)=prod_{m in D}(1+a_m cos 2pi m tau)
>= 0 with int R = 1 and int M*R < 1, where M(tau)=sum_v 1_{||v tau||<=1/14} is the covering multiplicity.
Since a TIGHT S has M>=1 a.e. => int M*R >= int R = 1, the strict inequality int M*R < 1 CERTIFIES S is
loose (the lonely set has positive measure).  Main term int M*R = sum_v s(0) = 13/7 = 1.857; the signed
Riesz corrections on the relation frequencies pull it down.  Hand-built R reached 1.41; target < 1.

This OPTIMIZES the Riesz product (dissociated D + coefficients a) to MINIMIZE int M*R for hard 13-sets
(the L-extremizers), by coordinate descent on a, over several candidate dissociated D tuned to S's
relation lattice.  Reports the minimum int M*R reached and whether the certificate (<1) fires.
"""
import numpy as np
from itertools import combinations

NG = 200003          # fine tau-grid (prime, avoids aliasing with integer speeds)
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14         # danger radius


def Mmult(S):
    """covering multiplicity M(tau) = #{v : ||v tau|| <= 1/14} on the grid."""
    M = np.zeros(NG)
    for v in S:
        d = np.abs(((v * TAU + 0.5) % 1.0) - 0.5)   # ||v tau||
        M += (d <= H).astype(float)
    return M


def cosbasis(D):
    """precompute cos(2 pi m tau) for m in D."""
    return {m: np.cos(2 * np.pi * m * TAU) for m in D}


def riesz(D, a, cb):
    R = np.ones(NG)
    for m, am in zip(D, a):
        R = R * (1 + am * cb[m])
    return R


def ratioMR(M, D, a, cb):
    """the VALID certificate objective: int(M*R)/int(R) for R>=0.  <1 => lonely set has positive
    measure => S loose.  (Tight S has M>=1 a.e. => ratio>=1, so no false positive.)"""
    R = riesz(D, a, cb)
    r = float(R.mean())
    return float((M * R).mean()) / r, r


def optimize(S, D, iters=40):
    """coordinate descent on a in [-1,1]^|D| to minimize the RATIO int(M*R)/int(R) (R>=0 via |a_m|<=1)."""
    M = Mmult(S)
    cb = cosbasis(D)
    a = np.zeros(len(D))
    base = float(M.mean())        # = int M = 13/7
    for _ in range(iters):
        improved = False
        for i in range(len(D)):
            best_ai = a[i]; best_val = ratioMR(M, D, a, cb)[0]
            for cand in np.linspace(-0.999, 0.999, 61):
                a[i] = cand
                val, _ = ratioMR(M, D, a, cb)
                if val < best_val - 1e-9:
                    best_val = val; best_ai = cand; improved = True
            a[i] = best_ai
        if not improved:
            break
    val, rnorm = ratioMR(M, D, a, cb)
    return val, rnorm, a, base


def relations_freqs(S, kmax=2):
    """small relation frequencies: differences and near-relations, as dissociated-set candidates."""
    S = sorted(S)
    diffs = sorted(set(abs(a - b) for a, b in combinations(S, 2) if a != b))
    return diffs


print("=" * 96)
print("RIESZ-PRODUCT certificate: minimize int M*R (target <1 => S loose).  main int M = 13/7 = 1.857")
print("=" * 96)
extremizers = {
    "{1..13}\\{6} U {56} [LOOSE L~.0056]": (sorted(set(range(1, 14)) - {6} | {56}), "loose"),
    "{1..12} U {182}     [TIGHT L=0]    ": (list(range(1, 13)) + [182], "tight"),
    "7*{1..12} U {13}    [LOOSE L~.029] ": (sorted(set([7 * i for i in range(1, 13)] + [13])), "loose"),
    "{1..13}             [ref]          ": (list(range(1, 14)), "?"),
}
print(f"  {'S':>38} {'best ratio int(MR)/int(R)':>26} {'bestD':>9} {'verdict':>16}")
for name, (S, kind) in extremizers.items():
    diffs = relations_freqs(S)
    cands = [("speeds", sorted(S)[:7]), ("lacunary", [1, 2, 4, 8, 16, 32, 64]),
             ("diffs", diffs[:7]), ("speeds+lac", sorted(S)[:4] + [16, 32, 64])]
    best = None
    for label, D in cands:
        val, rnorm, a, base = optimize(S, D, iters=25)
        if best is None or val < best[0]:
            best = (val, rnorm, label, D, a)
    val, rnorm, label, D, a = best
    verdict = "CERTIFIES LOOSE" if val < 1.0 - 1e-6 else ("(tight: >=1 OK)" if kind == "tight" else f"only {val:.3f}")
    print(f"  {name:>38} {val:>26.4f} {label:>9} {verdict:>16}")
print()
print("  READING: ratio<1 => positive-measure lonely set => S LOOSE (a valid certificate). TIGHT set")
print("  {1..12}U{182} must give ratio>=1 (VALIDITY check -- no false positive; M>=1 a.e.). How far")
print("  coordinate-descent Riesz reaches on LOOSE extremizers = the method's strength vs 2025 optimal.")
print("=" * 96)
