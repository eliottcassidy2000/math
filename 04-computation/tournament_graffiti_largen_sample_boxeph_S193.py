#!/usr/bin/env python3
"""tournament_graffiti_largen_sample_boxeph_S193.py -- boxeph-2026-07-21-S193
Guard the 'lifts to all n' survivors against small-n artifacts: random-sample n=8..12 and check
the candidate inequalities that passed n<=7. (kind-pasteur v2 discipline: n<=7 filter + larger-n
cross-validation.) Uses a fixed LCG (Date/random unavailable), varied per (n,trial)."""
from itertools import combinations
import tournament_graffiti_coupling_boxeph_S193 as G

def lcg(seed):
    x = seed & 0xFFFFFFFF
    while True:
        x = (1103515245*x + 12345) & 0x7FFFFFFF
        yield x

def rand_tournament(n, gen):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if next(gen) & 1: A[i][j] = 1
            else:             A[j][i] = 1
    return A

cands = [
    ("srange <= tr + 1",   lambda d: d["srange"] <= d["tr"] + 1),
    ("srange <= tr",       lambda d: d["srange"] <= d["tr"]),           # known-broken control (should FAIL)
    ("c3 <= H",            lambda d: d["c3"] <= d["H"]),
    ("srange <= 2*(tr-1)", lambda d: d["srange"] <= 2*(d["tr"]-1)),
    ("tr <= smax + 1",     lambda d: d["tr"] <= d["smax"] + 1),
    ("kings <= 2*c3 + 1",  lambda d: d["kings"] <= 2*d["c3"] + 1),
]
gen = lcg(20260721)
TRIALS = 400
print("random-sample validation, %d trials per n" % TRIALS)
for n in range(8, 13):
    fails = {name: None for name, _ in cands}
    for t in range(TRIALS):
        A = rand_tournament(n, gen)
        d = G.invariants(A, n)
        for name, pred in cands:
            if fails[name] is None and not pred(d):
                fails[name] = dict(n=d["n"], c3=d["c3"], H=d["H"], tr=d["tr"], srange=d["srange"], kings=d["kings"])
    print("\n n=%d:" % n)
    for name, _ in cands:
        f = fails[name]
        if f is None:
            print("   OK   %-20s (no counterexample in %d samples)" % (name, TRIALS))
        else:
            print("   FAIL %-20s witness c3=%d H=%d tr=%d srange=%d kings=%d"
                  % (name, f["c3"], f["H"], f["tr"], f["srange"], f["kings"]))
