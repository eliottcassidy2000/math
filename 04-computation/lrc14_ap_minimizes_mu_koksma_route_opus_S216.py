"""
opus-2026-07-11-S216: validating the convergence point of LRC(14).

FINDING (from 4 concept scouts + THM-685/661/530/657): every route to LRC(14) — the OffLine/Fourier
floor, the t>=3 signed cascade, the Kronecker transfer, the density-floor union bound, the two-scale
mu_infty program — reduces to ONE extremal lemma:

    "AP / consecutive minimizes the safe measure mu(S)"    (THM-530/657/527, verified k>=8, NOT proved).

The t>=3 signed cancellation (THM-684) is SUPERSEDED: the layer series is order-one/alternating/
non-truncatable (THM-685(B)); the transfer |LM(q) - q*mu(S)| <= Sum v (THM-685(A), Lean-formalized)
bypasses it, moving all remaining content to the measure floor mu(S) > 0.

This script validates the CONCRETE closure route for the SPREAD branch of "AP minimizes mu":
  raw safe measure   mu(S) = meas{ alpha in [0,1) : frac(v_i alpha) in [1/14,13/14] for all i }
  pair object        A(a,b) = meas{ alpha : frac(a alpha), frac(b alpha) both in [1/14,13/14] }
  per-pair deviation c(a,b) = A(a,b) - (6/7)^2, PROVED-bounded by Koksma (THM-686(C)):
                     |c(a,b)| <= (24/7)/max(a,b).
Claim under test: (1) AP {1..13} is the mu-minimizer (mu ~ 0, an isolated tight point);
(2) spread/dissociated families have mu ~ (6/7)^13 > 0; (3) the deviation mu(S) - (6/7)^13 is
controlled by the pairwise Koksma sum Sum_{pairs} (24/7)/max(a,b) -> 0 as the family spreads,
so the SPREAD branch decorrelates at rate Sum 1/V. Compact branch = exact census (THM-661). Both = the lemma.
"""
from fractions import Fraction as F
from itertools import combinations

LO = F(1, 14)   # band [1/14, 13/14]
HI = F(13, 14)
W  = HI - LO    # = 6/7

def safe_measure(v):
    """Exact meas{ alpha in [0,1) : frac(v_i alpha) in [LO,HI] for all i }, via a rational breakpoint sweep.
    Breakpoints: alpha where v_i*alpha = k+LO or k+HI (band edges) for some integer k, 0<=alpha<1."""
    bps = {F(0), F(1)}
    for vi in v:
        vi = abs(vi)
        if vi == 0:
            continue
        k = 0
        while True:
            a1 = F(k) / vi + LO / vi
            a2 = F(k) / vi + HI / vi
            if a1 >= 1 and a2 >= 1:
                break
            if 0 <= a1 < 1: bps.add(a1)
            if 0 <= a2 < 1: bps.add(a2)
            k += 1
    pts = sorted(bps)
    total = F(0)
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i+1]
        if b <= a:
            continue
        mid = (a + b) / 2
        ok = True
        for vi in v:
            fr = (vi * mid) % 1     # Fraction % 1 gives fractional part in [0,1)
            if not (LO <= fr <= HI):
                ok = False
                break
        if ok:
            total += (b - a)
    return total

def pair_A(a, b):
    """Exact A(a,b) = meas{ alpha : frac(a alpha), frac(b alpha) in [LO,HI] }."""
    return safe_measure([a, b])

fams = {
    "AP {1..13}            ": list(range(1, 14)),
    "AP dilated 3*{1..13}  ": [3*i for i in range(1, 14)],
    "near-AP (one detune)  ": list(range(1, 13)) + [20],
    "mild spread           ": [1, 3, 5, 8, 11, 15, 19, 24, 29, 35, 41, 48, 55],
    "GEN dissociated       ": [1, 5, 11, 17, 23, 28, 36, 41, 49, 55, 61, 67, 73],
    "wide dissociated      ": [7, 17, 31, 47, 67, 89, 113, 149, 181, 223, 271, 331, 401],
}

iid = (W) ** 13   # (6/7)^13, the decorrelated value
print(f"(6/7)^13 = {float(iid):.6f}   (6/7)^2 = {float(W**2):.6f}\n")
print(f"{'family':>24} {'mu(S) exact':>14} {'mu-(6/7)^13':>12} {'Sum|c_pair|':>12} {'Koksma bnd':>11} {'max v':>7}")
for name, v in fams.items():
    mu = safe_measure(v)
    dev = mu - iid
    # pairwise deviations + Koksma bound
    csum = F(0); ksum = F(0)
    for a, b in combinations(v, 2):
        c = pair_A(a, b) - W**2
        csum += abs(c)
        ksum += F(24, 7) / max(a, b)
    print(f"{name:>24} {float(mu):>14.6f} {float(dev):>+12.6f} {float(csum):>12.5f} {float(ksum):>11.4f} {max(v):>7}")

print("\n--- per-pair Koksma sharpness check: |A(a,b)-(6/7)^2| <= (24/7)/max(a,b)? ---")
viol = 0
for (a, b) in [(1,2),(1,13),(3,7),(5,11),(7,31),(1,100),(50,99)]:
    c = pair_A(a, b) - W**2
    bnd = F(24,7)/max(a,b)
    flag = "" if abs(c) <= bnd else "  <-- VIOLATION"
    viol += (abs(c) > bnd)
    print(f"  ({a:>3},{b:>3}): |c| = {float(abs(c)):.6f}   (24/7)/max = {float(bnd):.6f}{flag}")
print(f"\nKoksma violations: {viol}")
