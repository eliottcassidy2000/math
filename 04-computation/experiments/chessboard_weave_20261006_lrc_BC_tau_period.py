#!/usr/bin/env python3
"""Tasks B and C (chessboard weave, LRC reading), 2026-10-06.  Exact arithmetic.

B  FIRST LONELY TIME  tau(v) = min{ t in (0,1) : min_i ||t v_i|| >= 1/(n+1) }
   (left endpoint of the first component of Safe(v); computed by exact interval
   intersection, cross-checked against a brute grid and a fixed-point chase in
   the core self-test).
C  GOOD PERIOD  q(v) = least q such that some k/q (0<k<q) is in Safe(v)
   (= least Stern-Brocot denominator over the components of Safe(v)).

Universes: all primitive (gcd 1) n-subsets of {1..B}.
   required : n=2 B=200, n=3 B=60, n=4 B=30, n=5 B=18
   extended : n=2 B=400, n=3 B=150, n=4 B=60, n=5 B=36, n=6 B=24, n=7 B=20
Reproduce: python3 chessboard_weave_20261006_lrc_BC_tau_period.py > chessboard_weave_20261006_lrc_BC_tau_period.out
"""
import sys
import time
from fractions import Fraction
from collections import defaultdict, Counter

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from chessboard_weave_20261006_lrc_core import (safe_components, good_period_from,
                                                primitive_sets, check)

TOPK = 10
CAP = 40          # sets stored per top value


def census(n, B):
    N = n + 1
    top = {}                       # tau value -> [count, list of sets (capped)]
    top_no1 = {}                   # same, restricted to v_min >= 2
    qhist = Counter()
    qtop = {}                      # q -> [count, sets]
    rel = Counter()                # denominator of tau vs q
    first_chance = 0               # tau == 1/(N v_min)
    maxq_by_vmax = defaultdict(int)
    num = 0
    for v in primitive_sets(n, B):
        num += 1
        comps, D = safe_components(v)
        check(comps, f"LRC fails at {v}")
        t = Fraction(comps[0][0], D)
        # B checks
        check(t > 0, v)
        if n >= 2:
            check(2 * t < 1, f"tau >= 1/2 at {v}")
        lb = Fraction(1, N * v[0])
        check(t >= lb, f"tau < 1/(N vmin) at {v}")
        if t == lb:
            first_chance += 1
        for store, ok in ((top, True), (top_no1, v[0] >= 2)):
            if not ok:
                continue
            if t in store:
                store[t][0] += 1
                if len(store[t][1]) < CAP:
                    store[t][1].append(v)
            elif len(store) < TOPK or t > min(store):
                store[t] = [1, [v]]
                if len(store) > TOPK:
                    del store[min(store)]
        # C
        q, tq = good_period_from(comps, D)
        qhist[q] += 1
        maxq_by_vmax[v[-1]] = max(maxq_by_vmax[v[-1]], q)
        if q in qtop:
            qtop[q][0] += 1
            if len(qtop[q][1]) < CAP:
                qtop[q][1].append(v)
        elif len(qtop) < 4 or q > min(qtop):
            qtop[q] = [1, [v]]
            if len(qtop) > 4:
                del qtop[min(qtop)]
        dt = t.denominator
        check(dt >= q, f"denominator of tau < q at {v}")
        rel["den(tau) == q" if dt == q else "den(tau) > q"] += 1
    return dict(num=num, top=top, top_no1=top_no1, qhist=qhist, qtop=qtop, rel=rel,
                first_chance=first_chance, maxq_by_vmax=maxq_by_vmax)


def fmt_sets(lst, cnt, k=8):
    s = ", ".join(str(x) for x in lst[:k])
    if cnt > k:
        s += f", ... ({cnt} sets)"
    return s


UNIVERSES = [("required", {2: 200, 3: 60, 4: 30, 5: 18}),
             ("extended", {2: 400, 3: 150, 4: 60, 5: 36, 6: 24, 7: 20})]

T0 = time.time()
results = {}
print("B. FIRST LONELY TIME tau(v)   /   C. GOOD PERIOD q(v)")
print("   universes = all primitive n-subsets of {1..B}; exact arithmetic")
for label, U in UNIVERSES:
    for n, B in U.items():
        t1 = time.time()
        R = census(n, B)
        results[(label, n)] = R
        N = n + 1
        print()
        print("=" * 78)
        print(f"[{label}] n = {n}, max speed <= {B}: {R['num']} primitive sets "
              f"({time.time() - t1:.1f}s); threshold 1/{N}")
        print(f"  checks passed: Safe nonempty (LRC) for all; 0 < tau < 1/2 for all; "
              f"tau >= 1/(N*vmin) for all")
        print(f"  tau == 1/(N*vmin) ('lonely at first chance'): {R['first_chance']} "
              f"of {R['num']}")
        top = R["top"]
        ks = sorted(top, reverse=True)
        print(f"  MAX tau = {ks[0]} = {float(ks[0]):.6f}; extremizers: {top[ks[0]][1]}")
        print(f"  top {TOPK} values of tau:")
        for k in ks:
            c, lst = top[k]
            has1 = sum(1 for x in lst if x[0] == 1)
            hasN = sum(1 for x in lst if any(y % N == 0 for y in x))
            print(f"    {str(k):>9} = {float(k):.6f}  [{c} sets; stored {len(lst)}: "
                  f"contain 1: {has1}, contain a multiple of {N}: {hasN}]  {fmt_sets(lst, c)}")
        tn = R["top_no1"]
        if tn:
            k1 = max(tn)
            print(f"  max tau over sets with vmin >= 2: {k1} = {float(k1):.6f} at "
                  f"{fmt_sets(tn[k1][1], tn[k1][0])}")
        # C
        qh = R["qhist"]
        print(f"  C: q(v) histogram: {dict(sorted(qh.items()))}")
        qk = sorted(R["qtop"], reverse=True)
        print(f"  C: MAX q = {qk[0]} at {fmt_sets(R['qtop'][qk[0]][1], R['qtop'][qk[0]][0], 12)}")
        for k in qk[1:]:
            print(f"     q = {k}: {fmt_sets(R['qtop'][k][1], R['qtop'][k][0], 6)}")
        mb = R["maxq_by_vmax"]
        cum, run = [], 0
        for b in range(1, B + 1):
            run = max(run, mb.get(b, 0))
            if b in (10, 20, 30, 50, 100, 150, 200, 300, 400) or b == B:
                cum.append(f"B={b}:{run}")
        print(f"  C: running max q by max speed: {', '.join(cum)}")
        print(f"  C: denominator of tau vs q(v): {dict(R['rel'])} (den(tau) < q is impossible: tau is in Safe)")
print()
print(f"Total time {time.time() - T0:.1f}s")
