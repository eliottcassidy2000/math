"""
mac-mini-2026-07-07-S56 (HYP-5207) -- THE CONDITIONAL-TENT c-TABLE (THM-651's named
k=9/10 frontier, kps-S73).

OBJECT: c(d, P) = int_{G_P} f_k(frac(d x)) dx / (meas(G_P) * int f_k), where
  f_k(s) = (s - beta_k)_+ on (0, 1/7], beta_k = (14-k)/(7k)
  (k=9: beta = 5/63, support (5/63, 1/7]; k=10: beta = 4/70 = 2/35).
kps-S73's conditional tent discharges the k-leg if sup_d c(d, P) <= c*_k with
c*_9 = 1.7, c*_10 = 1.29 (their verified arithmetic).

EXACT COMPUTATION: G_P = finite union of rational intervals [l, h].  For winding m,
frac(dx) in (beta, 1/7] <=> x in ((m+beta)/d, (m+1/7)/d]; on the intersection with
[l, h], f = d x - m - beta integrates to an exact rational (quadratic antiderivative).
Sum over m = 0..d-1 and the intervals of G_P.

RUNS: (1) the two k=9-relevant binding P's (|P| = 4: the hard core {10,11,12,13} and
the min-meas P) and k=10's (|P| = 3): exact c(d, P) for d = 1..250; report sup, the
resonant d's (coprime/ladder structure), and the verdict vs c*.
(2) Koksma-tail sanity: c(d) -> 1 as d grows (rate ~ #intervals/d) -- measure the decay.
"""
from fractions import Fraction as F
from math import gcd

def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14*p)
        for j in range(p+1):
            bad.append((F(j,p)-w, F(j,p)+w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort()
    merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else:
            merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def tent_integral_on_GP(P_iv, d, beta, th=F(1,7)):
    """int_{G_P} (frac(dx) - beta)_+ 1[frac(dx) <= th] dx, exact."""
    tot = F(0)
    d = int(d)
    for (l, h) in P_iv:
        # windings m with window ((m+beta)/d, (m+th)/d] intersecting [l, h]
        m_lo = int((l*d - th))  # conservative
        m_hi = int((h*d)) + 1
        for m in range(max(0, m_lo), m_hi+1):
            wl, wh = (F(m)+beta)/d, (F(m)+th)/d
            a, b = max(wl, l), min(wh, h)
            if a >= b: continue
            # integral of (d x - m - beta) dx from a to b
            tot += (d*(b*b - a*a)/2 - (F(m)+beta)*(b - a))
    return tot

def c_table(P, k, dmax=250):
    beta = F(14-k, 7*k)
    th = F(1,7)
    intf = (th - beta)**2 / 2
    iv = GP_intervals(P)
    mGP = sum(h-l for l, h in iv)
    denom = mGP * intf
    rows = []
    for d in range(1, dmax+1):
        num = tent_integral_on_GP(iv, d, beta)
        rows.append((d, num/denom))
    return rows, mGP, len(iv)

print("=== conditional-tent c(d,P): exact table (THM-651 frontier) ===")
CASES = [
    (9,  (10,11,12,13), 1.7),
    (9,  (1,11,12,13),  1.7),     # a min-meas-style |P|=4 comparison
    (10, (11,12,13),    1.29),
    (10, (1,12,13),     1.29),
]
for k, P, cstar in CASES:
    rows, mGP, nIv = c_table(P, k, dmax=250)
    vals = [float(c) for d, c in rows]
    supc = max(vals); supd = rows[vals.index(supc)][0]
    over = [(d, float(c)) for d, c in rows if float(c) > cstar]
    print(f"\n  k={k}, P={P} (meas GP = {float(mGP):.4f}, {nIv} intervals): "
          f"sup_d<=250 c = {supc:.4f} at d = {supd}; c* = {cstar}")
    print(f"    d's exceeding c*: {len(over)}"
          + (f" -> {over[:10]}" if over else "  => CONDITION HOLDS on the small-d table"))
    # resonance structure of the top offenders
    top = sorted(rows, key=lambda r: -float(r[1]))[:8]
    print(f"    top c values: " + " ".join(f"d={d}:{float(c):.3f}" for d, c in top))
    # ladder/coprime signature: gcd of top d's with P's elements and with 14
    for d, c in top[:4]:
        sig = [gcd(d, p) for p in P] + [gcd(d, 14)]
        print(f"      d={d}: gcd with P + 14 = {sig}, c = {float(c):.4f}")
    # decay sanity
    tail = [float(c) for d, c in rows if d > 200]
    print(f"    tail d in (200,250]: c in [{min(tail):.4f}, {max(tail):.4f}] (Koksma -> 1)")
