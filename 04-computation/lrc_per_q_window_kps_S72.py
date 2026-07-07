#!/usr/bin/env python3
r"""
lrc_per_q_window_kps_S72.py   (kind-pasteur-2026-07-07-S72, HYP-5117)

WORK THE OPEN ROUTE (R2 / density-floor AP-minimality) via the PER-RESONANCE-q window
decomposition.  Goal: reduce the sigma-even floor 'AP minimizes mu_{1/7}' to a FINITE
per-q comparison (a sigma-odd handle), or find the precise obstruction.

FRAME (mac-mini-S15 + opus-S134): a circular gap > 1/7 needs the phases to cluster, which
happens near x = p/q with q <= 6 (a gap of 1/q > 1/7 needs q <= 6).  Near x = p/q the config
{frac(e_i x)} collapses onto the residues {e_i mod q}/q.  If E hits r distinct residues mod q,
the largest gap near p/q is ~ (max residue-gap)/q >= 1/q >= 1/6 > 1/7 when r <= 6.  The AP
{1..k} hits ALL q residues (r=q, minimal gap); a non-AP can MISS residues (r<q, WIDER gap).

DECOMPOSITION: mu(E) = sum over the maxgap-superlevel; ATTRIBUTE each x to the resonance q =
the denominator of the nearest low-order rational to the CENTER of the responsible gap.
W_q(E) = measure attributed to q.  TEST: does the AP minimize each W_q?

OUTPUTS:
 (1) per-q window decomposition W_q(E) for the AP and structured non-AP families (k=8,13).
 (2) does W_q(AP) <= W_q(E) for each q?  (=> AP-minimality reduces to per-q).
 (3) residue-deficit vs W_q excess, in PRIMITIVE normalization (gcd of differences = 1).
 (4) the dilation caveat: check mu(cE)=mu(E) but residue structure changes -> which frame
     is dilation-robust.
"""
import random, math
from fractions import Fraction as F

def config_maxgap_and_center(E, x):
    """return (maxgap, center-of-maxgap) for {frac(e_i x)}."""
    ph = sorted((e * x) % 1.0 for e in E); n = len(ph)
    best = -1; ctr = 0
    for i in range(n-1):
        g = ph[i+1] - ph[i]
        if g > best: best = g; ctr = (ph[i] + ph[i+1]) / 2
    gw = ph[0] + 1 - ph[-1]
    if gw > best: best = gw; ctr = ((ph[-1] + ph[0] + 1) / 2) % 1.0
    return best, ctr

def nearest_lowq(xv, Q=6):
    """denominator q<=Q of the nearest p/q to xv (the resonance the gap sits at)."""
    bestq, bestd = 1, 1.0
    for q in range(1, Q+1):
        for p in range(0, q+1):
            d = abs(xv - p/q)
            d = min(d, abs(xv - p/q - 1), abs(xv - p/q + 1))
            if d < bestd: bestd = d; bestq = q
    return bestq

def mu_and_windows(E, res=30000, Q=6):
    """mu(E) and per-q attribution W_q (by the resonance nearest the responsible gap's center)."""
    c = 0; W = {q: 0 for q in range(1, Q+1)}; other = 0
    for r in range(res):
        x = (r + .5) / res
        mg, ctr = config_maxgap_and_center(E, x)
        if mg > 1/7:
            c += 1
            # attribute by the resonance nearest x (the time), which drives the clustering
            q = nearest_lowq(x, Q)
            W[q] += 1
    return c/res, {q: W[q]/res for q in W}

Tk = {8: 0.6185, 13: 0.0565}
def residue_deficit(E, Q=6):
    """dilation-issue-aware: sum over q<=Q of (q - #distinct residues mod q), primitive form."""
    g = 0
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            g = math.gcd(g, E[j]-E[i])
    Ep = [ (e - E[0])//max(g,1) for e in E ]   # primitive
    tot = 0
    for q in range(2, Q+1):
        tot += q - len(set(e % q for e in Ep))
    return tot, Ep

print("=" * 96)
print("PER-q WINDOW DECOMPOSITION of mu_{1/7}: does the AP minimize each W_q?")
print("=" * 96)
for k in (8, 13):
    AP = list(range(1, k+1))
    muAP, WAP = mu_and_windows(AP)
    print(f"\n  k={k}: mu(AP)={muAP:.4f}, per-q windows W_q(AP): " +
          ", ".join(f"q{q}:{WAP[q]:.3f}" for q in range(2,7) if WAP[q]>0.001))
    fams = {
        "spread AP d=2": [1+2*j for j in range(k)],
        "AP-defect (drop mid,+1)": sorted(set(list(range(1,k+1))) - {k//2}) + [k+1],
        "miss-residue mod3 {1,2,4,5,7,8,..}": [e for e in range(1,3*k) if e%3 != 0][:k],
        "miss-residue mod2 (all odd)": [2*j-1 for j in range(1,k+1)],
        "primes": ([2,3,5,7,11,13,17,19,23,29,31,37,41])[:k],
        "geometric": [2**j for j in range(k)],
    }
    for r in range(3):
        fams[f"random-{r}"] = sorted(random.Random(k*10+r).sample(range(1,200), k))
    print(f"    {'family':30s} {'mu':>7} {'mu-muAP':>8} {'W_q per q (2..6)':>28} {'resdef':>7} {'per-q >= AP?':>12}")
    for nm, E in fams.items():
        mu, W = mu_and_windows(E)
        rd, _ = residue_deficit(E)
        perq_ok = all(W[q] >= WAP[q] - 0.006 for q in range(2,7))
        wstr = " ".join(f"{W[q]:.2f}" for q in range(2,7))
        print(f"    {nm:30s} {mu:7.4f} {mu-muAP:+8.4f} {wstr:>28} {rd:>7} {str(perq_ok):>12}")
    print(f"    => AP per-q-minimal (W_q(E) >= W_q(AP) for all q<=6, all families): tested above")

print()
print("=" * 96)
print("DILATION CAVEAT: mu(cE)=mu(E) but residues change -> is per-q window dilation-robust?")
print("=" * 96)
for k in (8,):
    E = list(range(1,k+1))
    for c in (1, 3, 5):
        cE = [c*e for e in E]
        mu, W = mu_and_windows(cE)
        rd, _ = residue_deficit(cE)
        print(f"    {c}*AP: mu={mu:.4f} (dilation-inv: {abs(mu-mu_and_windows(E)[0])<0.01}); "
              f"resdef(primitive)={rd}; W_q: " + " ".join(f"{W[q]:.2f}" for q in range(2,7)))
    print("    => mu dilation-invariant; primitive residue-deficit dilation-invariant (good);")
    print("       but the ATTRIBUTED q shifts under dilation (x->cx relabels resonances) -- the")
    print("       per-q window is NOT individually dilation-invariant, only the SUM mu is.")
print("DONE.")
