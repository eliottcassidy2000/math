#!/usr/bin/env python3
"""Audit G, THM-4609: (1) the Z_7 examples (rank 3 and rank 4), e = 1, T = 16384;
(2) counterexample hunt: random contracting translation-only maps of independent rank >= 3 on Z_5, Z_7, Z_11 --
does any of them coalesce (P(no merge by T) -> 0 with a non-decaying window of debt visits)?
(3) outside the theorem: translation-only rank-3 maps with a repeated non-unit multiplier (not 'independent')."""
import time, math, random, itertools
from mwsim import MW, run_pairs, fmt_run, factor
t0 = time.time()

def indep_rank(ms):
    primes = sorted({p for m in ms for p in factor(m)})
    rows = [[Fr(factor(m).get(p, 0)) for p in primes] for m in ms]
    # rank of differences
    diffs = [[a - b for a, b in zip(row, rows[0])] for row in rows[1:]]
    rk = 0; A = [r[:] for r in diffs]; ncol = len(primes)
    for c in range(ncol):
        pr = next((i for i in range(rk, len(A)) if A[i][c] != 0), None)
        if pr is None: continue
        A[rk], A[pr] = A[pr], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / A[rk][c]; A[i] = [x - f * y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk
from fractions import Fraction as Fr

def report(name, mw, e, T, N, seed):
    cps = [c for c in (64, 256, 1024, 4096, 16384) if c <= T]
    res = run_pairs(mw, e, T, N, seed, cps)
    Lam = sum(math.log(x / mw.d) for x in mw.m) / mw.d
    print(f"{name}: m = {mw.m}, r = {mw.r}, Lambda = {Lam:+.4f}, e = {e}, N = {N}  [{time.time()-t0:.0f}s]")
    for line in fmt_run(res, cps): print("    ", line)
    print(f"     nonzero-debt merges {res['nonzero_debt_merges']}; latest merge {max(res['merge_times']) if res['merge_times'] else None}", flush=True)
    return res

# (1) Z_7 examples
r7 = [0, 6, 5, 4, 3, 2, 1]
report("Z7 rank 3 (THM-4609 ex.)", MW(7, [1, 1, 1, 1, 8, 15, 22], r7), 1, 16384, 500, 71)
report("Z7 rank 4 (THM-4609 ex.)", MW(7, [1, 1, 1, 8, 15, 22, 29], r7), 1, 16384, 500, 72)
# (2) hunt
rnd = random.Random(99)
for d, nmaps in ((5, 8), (7, 8), (11, 4)):
    cand = [k * d + 1 for k in range(1, 30)]
    found = 0; tries = 0
    while found < nmaps and tries < 100000:
        tries += 1
        rho = rnd.choice([3, 3, 4]) if d > 5 else 3
        nonunit = rnd.sample(cand, rho)
        if math.prod(nonunit) >= d ** d: continue
        ms = [1] * (d - rho) + nonunit; rnd.shuffle(ms)
        if indep_rank(ms) != rho: continue
        # independence of the non-unit multipliers (rank must equal number of non-units)
        found += 1
        r = [(-ms[i] * i) % d + d * rnd.randint(-1, 1) for i in range(d)]
        report(f"hunt d={d} #{found} rank {rho}", MW(d, ms, r), rnd.choice([1, 2, 3, d]), 4096, 400, rnd.randrange(10 ** 9))
# (3) outside the theorem: repeated non-unit multiplier
report("Z7 (1,1,1,8,8,15,22) rank 3, NOT independent", MW(7, [1, 1, 1, 8, 8, 15, 22], r7), 1, 4096, 600, 73)
report("Z7 (1,1,8,1,8,15,22) rank 3, NOT independent", MW(7, [1, 1, 8, 1, 8, 15, 22], r7), 1, 4096, 600, 74)
report("Z5 (1,6,6,11,1)? rank 2 control", MW(5, [1, 6, 6, 11, 1], [0, (-6) % 5, (-12) % 5, (-33) % 5, (-4) % 5]), 1, 4096, 600, 75)
