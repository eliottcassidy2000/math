#!/usr/bin/env python3
"""
klein-2026-07-02-S118 (HYP-4023b) - c=8 CLUSTERED adversary vs bounded-combo citations.

Question: can a clustered c=8 block (speeds w+d_j, d_j small integers) COVER the whole
citation window for all t while DODGING every bounded integer-combo citation?

Model (clustered limit, offsets static on the window): positions psi_j = d_j*t0 mod 1,
arcs [psi_j - 1/14, psi_j + 1/14]. Adversary needs:
  (COVER) the 8 arcs cover the circle ROBUSTLY (max sorted gap of centers <= 1/7 - slack,
          slack > drift budget)
  (DODGE) for every nonzero combo c = sum a_j d_j with |a_j| <= A, |c| <= Cmax:
          dist(c*t0, Z) >= 1/14   (else we cite the combo and kill the config)
d_1 = 0 pinned (block reference); arc 1 sits at psi=0.

Search: small d-vectors (increasing, d_8 <= DMAX), t0 on a fine grid + local refine.
Output: for each d-vector, the best (max over t0) 'survival margin' =
  min( min-gap-slack requirement satisfied?, min over combos of dist(c t0) - 1/14 ).
If ALL d-vectors have survival < 0 => evidence the combo route closes clustered c=8.
Any survivor is printed with its full config = the honest obstruction.
"""
import itertools, math

H = 1/14
def frac(x): return x - math.floor(x)
def cdist(x): f = frac(x); return min(f, 1 - f)

def robust_cover_slack(psis):
    """slack = 1/7 - max gap between consecutive arc CENTERS (sorted, cyclic).
    >0 means arcs of width 1/7 cover with room."""
    ps = sorted(frac(p) for p in psis)
    gaps = [ps[(i+1) % len(ps)] - ps[i] + (1 if i == len(ps)-1 else 0) for i in range(len(ps))]
    return 1/7 - max(gaps)

def combos(ds, A, Cmax):
    """nonzero integer combos of ds with coeff height <= A and 0 < |value| <= Cmax."""
    seen = set()
    n = len(ds)
    for a in itertools.product(range(-A, A+1), repeat=n):
        c = sum(ai*di for ai, di in zip(a, ds))
        if c != 0 and abs(c) <= Cmax and c not in seen and -c not in seen:
            seen.add(c)
    return sorted(seen)

def survival(ds, A=2, Cmax=None, ngrid=20000, refine=6):
    """max over t0 of min(cover_slack, dodge margin). ds includes d_1=0."""
    if Cmax is None: Cmax = 8 * max(ds)  # bounded combos: within ~one block-spread
    cs = combos(ds, A, Cmax)
    P = max(cs) if cs else 1
    best = -1; best_t = None
    # t0 has period 1/gcd of active ds... sweep [0,1) since ds integers
    for k in range(ngrid):
        t0 = (k + 0.5) / ngrid
        psis = [d * t0 for d in ds]
        sl = robust_cover_slack(psis)
        if sl <= best: continue
        dodge = min(cdist(c * t0) - H for c in cs)
        m = min(sl, dodge)
        if m > best: best, best_t = m, t0
    # local refine around best_t
    for _ in range(refine):
        span = 1.0 / ngrid
        cand = [best_t + (i - 50) * span / 50 for i in range(101)]
        for t0 in cand:
            psis = [d * t0 for d in ds]
            m = min(robust_cover_slack(psis), min(cdist(c * t0) - H for c in cs))
            if m > best: best, best_t = m, t0
        ngrid *= 10
    return best, best_t, len(cs)

print("c=8 clustered search: d-vectors (0,d2..d8), coeff height A=2")
print(f"{'d-vector':>34} {'#combos':>8} {'survival':>10} {'verdict':>10}")
survivors = []
DMAX = 12
count = 0
import sys
for d_rest in itertools.combinations(range(1, DMAX + 1), 7):
    ds = (0,) + d_rest
    count += 1
    s, t0, nc = survival(ds, A=2, ngrid=4000, refine=3)
    verdict = "SURVIVES" if s > 0 else "killed"
    if s > 0:
        survivors.append((ds, s, t0))
        print(f"{str(ds):>34} {nc:>8} {s:>10.5f} {verdict:>10}  t0={t0:.6f}")
    if count % 100 == 0:
        print(f"  ... {count} vectors, {len(survivors)} survivors so far", flush=True)

print(f"\nTOTAL: {count} d-vectors, {len(survivors)} survive all A=2 bounded-combo citations")
if survivors:
    print("\nRe-testing survivors with height A=3 combos:")
    for ds, s, t0 in survivors[:40]:
        s3, t03, nc3 = survival(ds, A=3, ngrid=40000, refine=4)
        print(f"  {ds}: A=3 survival {s3:.5f} ({nc3} combos) {'STILL SURVIVES' if s3>0 else 'killed at A=3'}")
