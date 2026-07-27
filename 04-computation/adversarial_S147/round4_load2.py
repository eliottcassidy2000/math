#!/usr/bin/env python3
"""Round 4: exhaust the load-2 doorway.

Take the 5 load-2 survivors; for every 12-subset (one replaced) and
11-subset (two replaced), place tails via vectorized CRT admissibility over
[14,60], then fold checks (all breakpoint moduli >= 55 tied to tails), then
exact M.  Tails range: h in [55, 8500] with v_max >= 65 enforced at the end.
"""
import sys, time
import numpy as np
from math import gcd
from fractions import Fraction
from itertools import combinations
sys.path.insert(0, ".")
from lrc_engine import exact_M, covering_check, d_gap, W_of_q, THR, FLOOR
from round2_crt_tail import pass_at_q

QLO, QHI = 14, 60
HMAX = 8500

def qdata(C):
    """per q in [QLO,QHI]: None if covered/pinned by C alone; else
    (U, ok1 bool-array over residues, pairfix bool) where ok1[r] = residue r
    alone finishes q; also store A sets for pair logic."""
    out = {}
    for q in range(QLO, QHI + 1):
        d = d_gap(q)
        if any(v % q == 0 for v in C):
            continue
        units = [a for a in range(1, q) if gcd(a, q) == 1]
        U = [a for a in units
             if min(min((a * v) % q, q - (a * v) % q) for v in C) > d]
        if not U:
            continue
        A = []
        for r in range(q):
            s = set()
            for a in U:
                x = (a * r) % q
                if min(x, q - x) <= d:
                    s.add(a)
            A.append(s)
        Uset = set(U)
        ok1 = np.array([A[r] >= Uset for r in range(q)], dtype=bool)
        out[q] = (Uset, A, ok1)
    return out

def admissible_single(qd, hmax=HMAX):
    """bool array adm[h] for h in [0,hmax]: every active q finished by h."""
    adm = np.ones(hmax + 1, dtype=bool)
    hs = np.arange(hmax + 1)
    for q, (U, A, ok1) in qd.items():
        adm &= ok1[hs % q]
    return adm

def pair_admissible(qd, h1, hmax=HMAX):
    """bool array over h2: for every active q, U <= A[h1%q] u A[h2%q]."""
    adm = np.ones(hmax + 1, dtype=bool)
    hs = np.arange(hmax + 1)
    for q, (U, A, ok1) in qd.items():
        rem = U - A[h1 % q]
        if not rem:
            continue
        ok2 = np.array([rem <= A[r] for r in range(q)], dtype=bool)
        adm &= ok2[hs % q]
    return adm

def fold_and_M(V, log, label):
    V = sorted(V)
    g = 0
    for x in V:
        g = gcd(g, x)
    if g != 1 or len(set(V)) != 13 or V[-1] < 65 or covering_check(V):
        return None
    folds = set()
    for i, a in enumerate(V):
        for b in V[i+1:]:
            if a + b >= 55:
                folds.add(a + b)
            if b - a >= 55:
                folds.add(b - a)
        folds.add(2 * a)
    for q in sorted(folds):
        W, k = W_of_q(V, q)
        if W > d_gap(q):
            log.append((tuple(V), q, k, W, d_gap(q), label))
            return None
    M, (q, k) = exact_M(V)
    return (M, tuple(V), q, k, label)

def cover_need(C):
    return [q for q in range(2, 14) if all(v % q for v in C)]

def main():
    LOAD2 = [
        [1, 4, 8, 12, 16, 20, 28, 32, 36, 40, 44, 48, 52],
        [4, 7, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 52],
        [4, 8, 12, 16, 20, 24, 27, 28, 32, 36, 40, 44, 52],
        [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 41, 44, 52],
        [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 52, 55],
    ]
    results, log = [], []
    t0 = time.time()
    # ---- Option A: 12-core + one tail
    n_adm = 0
    for S in LOAD2:
        for i in range(13):
            C = [v for j, v in enumerate(S) if j != i]
            need = cover_need(C)
            qd = qdata(C)
            if any(not v[2].any() for v in qd.values()):
                continue
            adm = admissible_single(qd)
            for h in range(55, HMAX + 1):
                if not adm[h] or h in C:
                    continue
                if any(h % x for x in need):
                    continue
                n_adm += 1
                r = fold_and_M(C + [h], log, f"A:{S[:2]}..drop{S[i]}")
                if r:
                    results.append(r)
                    M = r[0]
                    tag = "*** IN GAP ***" if FLOOR < M < THR else ""
                    print(f"A-PASS M={M} V={r[1]} wit {r[3]}/{r[2]} {tag}",
                          flush=True)
    print(f"[optA done {time.time()-t0:.0f}s] admissible tails tried: {n_adm}, "
          f"full-pass: {len(results)}, fold-killed: {len(log)}", flush=True)
    # ---- Option B: 11-core + two tails
    nB = 0
    for S in LOAD2:
        for i, j in combinations(range(13), 2):
            C = [v for k2, v in enumerate(S) if k2 not in (i, j)]
            need = cover_need(C)
            qd = qdata(C)
            # quick pair feasibility: every q must be pair-coverable
            feas = True
            for q, (U, A, ok1) in qd.items():
                if ok1.any():
                    continue
                good = False
                for r1 in range(q):
                    rem = U - A[r1]
                    if any(rem <= A[r2] for r2 in range(q)):
                        good = True
                        break
                if not good:
                    feas = False
                    break
            if not feas:
                continue
            for h1 in range(55, HMAX + 1):
                if h1 in C:
                    continue
                adm2 = pair_admissible(qd, h1)
                lo = max(65, h1 + 1)
                cand = np.nonzero(adm2[lo:])[0] + lo
                for h2 in cand:
                    h2 = int(h2)
                    if h2 in C:
                        continue
                    need2 = [x for x in need if h1 % x and h2 % x]
                    if need2:
                        continue
                    nB += 1
                    r = fold_and_M(C + [h1, h2], log,
                                   f"B:drop({S[i]},{S[j]})")
                    if r:
                        results.append(r)
                        M = r[0]
                        tag = "*** IN GAP ***" if FLOOR < M < THR else ""
                        print(f"B-PASS M={M} V={r[1]} wit {r[3]}/{r[2]} {tag}",
                              flush=True)
        print(f"  [S={S[:3]}.. done {time.time()-t0:.0f}s B-tried={nB} "
              f"pass={len(results)} foldkill={len(log)}]", flush=True)
    # report
    from collections import Counter
    print(f"\nTOTAL full-pass {len(results)}, fold-killed {len(log)}")
    kf = Counter(q for _, q, *_ in log)
    print("fold-kill moduli (top 15):", kf.most_common(15))
    results.sort()
    for M, V, q, k, lab in results[:25]:
        tag = "*** IN GAP ***" if FLOOR < M < THR else ""
        print(f"M={str(M):>9}={float(M):.6f} wit {k}/{q} [{lab}] V={list(V)} {tag}")
    with open("round4_foldlog.txt", "w") as f:
        for V, q, k, W, dg, lab in log:
            f.write(f"{list(V)} q={q} k={k} W={W} need<={dg} [{lab}]\n")

if __name__ == "__main__":
    main()
