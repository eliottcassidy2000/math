#!/usr/bin/env python3
"""
lrc14_peeling_kps_S13.py -- kind-pasteur-2026-07-05-S13, HYP-4157.

THE PEELING LEMMA for the loose branch of TightLooseDichotomy, and the reduction
to CRITICAL configs.  M(B) = max_t min_{v in B} ||v t|| for a 12-tuple B.

Peeling lemma: if W subset B (|W|<=6) and t0 is an optimal time of B\\W with every
w in W at ||w t0|| >= 2/25, then M(B) >= 2/25 -- because B\\W has M >= 1/(13-|W|) > 2/25
by CITED LRC(13-|W|), so min over B at t0 is >= 2/25.

Non-critical (some redundant runner) => certified loose at depth 1.  The residual =
CRITICAL configs (every runner essential); gap-violators are necessarily critical;
empirically the only critical config with M < 2/25 is the dilated AP.

Verifies: (1) iterated peeling coverage; (2) tight APs not certified; (3) critical
census is gap-free.
"""
from math import gcd
from functools import reduce
from itertools import combinations
from fractions import Fraction as Fr

def distZ(x, q): r = x % q; return min(r, q - r)

def opt_times(vs):
    Q = 2 * max(vs); best = Fr(0); opts = []
    for q in range(2, Q + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            mm = min(distZ(v * a, q) for v in vs); f = Fr(mm, q)
            if f > best: best = f; opts = [(a, q)]
            elif f == best: opts.append((a, q))
    return best, opts

def M_exact(vs): return opt_times(vs)[0]
def wdist(w, a, q): return Fr(distZ(w * a, q), q)
BETA = Fr(2, 25)

def cert_depth(vs, maxdrop=6):
    """smallest peel depth certifying M(vs) >= 2/25, or None."""
    n = len(vs)
    for k in range(1, maxdrop + 1):
        if 13 - k < 7: break                    # floor 1/(13-k) must exceed 2/25
        for drop in combinations(range(n), k):
            sub = [vs[i] for i in range(n) if i not in drop]
            dropped = [vs[i] for i in drop]
            _, opts = opt_times(sub)
            for a, q in opts:
                if all(wdist(w, a, q) >= BETA for w in dropped):
                    return k
    return None

def isAPdil(vs):
    s = sorted(vs); d = s[1] - s[0]
    return all(s[i + 1] - s[i] == d for i in range(11)) and s[0] == d

def is_critical(vs):
    M = M_exact(vs)
    for i in range(12):
        if M_exact(vs[:i] + vs[i + 1:]) <= M:
            return False, M
    return True, M

if __name__ == "__main__":
    gaplo, gaphi = Fr(1, 13), Fr(2, 25)
    print("=== (1) iterated peeling coverage, [1,16] exhaustive ===")
    depth_hist = {}; uncov = []
    for c in combinations(range(1, 17), 12):
        vs = list(c)
        if reduce(gcd, vs) != 1 or isAPdil(vs): continue
        d = cert_depth(vs)
        if d is None: uncov.append(vs)
        else: depth_hist[d] = depth_hist.get(d, 0) + 1
    print(f"  non-AP certified loose by peel-depth: {dict(sorted(depth_hist.items()))}")
    print(f"  uncertified non-AP: {len(uncov)}  {uncov[:3]}")

    print("=== (2) dilated APs must NOT be certified (they are tight) ===")
    for cc in [1, 2, 3, 5]:
        vs = [cc * i for i in range(1, 13)]
        print(f"  {cc}*{{1..12}}: cert_depth={cert_depth(vs)}  M={M_exact(vs)}")

    print("=== (3) critical census is gap-free ([1,16]) ===")
    ncrit = 0; crit_lt = []
    for c in combinations(range(1, 17), 12):
        vs = list(c)
        if reduce(gcd, vs) != 1: continue
        isc, M = is_critical(vs)
        if isc:
            ncrit += 1
            if M < gaphi: crit_lt.append((str(M), vs, isAPdil(vs)))
    print(f"  critical configs: {ncrit}; with M<2/25: {crit_lt}")
    print(f"  critical IN gap (1/13,2/25): {sum(1 for m,_,_ in crit_lt if gaplo<Fr(*map(int,m.split('/')))<gaphi)}")
