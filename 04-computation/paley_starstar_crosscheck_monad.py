#!/usr/bin/env python3
"""
paley_starstar_crosscheck_monad.py
monad-explorer-2026-06-07 (8th session)

CROSS-VALIDATION: my fast INTEGER even-series test (fundamental cycles) vs the
original SVD-based test (paley_catalan_star_star_monad.py).  Resolves the k=6
count discrepancy (my fast enumerator: 2351 even-series patterns; prior-session
OEIS extrapolation A215257: 2345 -- NEVER actually computed at k=6).

Runs BOTH classifiers on every partition (k<=KFULL exhaustive) and on a random
sample at k=KSAMPLE; reports ANY disagreement.  Zero disagreements => my fast
count is correct and the A215257 identification needs re-checking at k=6.
"""
import sys
from math import factorial
from collections import defaultdict
import numpy as np

# ---------- original SVD-based test (verbatim logic) ----------
def edge_flow_lines_svd(edges, nb):
    E = len(edges)
    Bm = np.zeros((nb, E), dtype=np.float64)
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1.0; Bm[u, ei] -= 1.0
    _, s, vh = np.linalg.svd(Bm)
    tol = 1e-9
    rank = int((s > tol).sum()); m = E - rank
    if m == 0:
        return [tuple()] * E, 0
    ns = vh[rank:]
    lines = []
    for e in range(E):
        v = ns[:, e]
        if np.max(np.abs(v)) < 1e-7:
            lines.append(("ZERO",)); continue
        v = v / np.max(np.abs(v))
        for x in v:
            if abs(x) > 1e-7:
                if x < 0: v = -v
                break
        lines.append(tuple(round(float(x), 6) for x in v))
    return lines, m

def is_even_series_svd(edges, nb):
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; stk = [0]
    while stk:
        x = stk.pop()
        for w in adj[x]:
            if w not in seen: seen.add(w); stk.append(w)
    if len(seen) != nb:
        return False
    lines, m = edge_flow_lines_svd(edges, nb)
    if m == 0: return False
    if any(l == ("ZERO",) for l in lines): return False
    groups = defaultdict(int)
    for l in lines: groups[l] += 1
    return all(c % 2 == 0 for c in groups.values())

# ---------- my fast integer test ----------
sys.path.insert(0, "04-computation")
from paley_starstar_triangle_fast_monad import analyze, rgs_iter

def fast_is_even_series(a, L):
    return analyze(a, L) is not None

def build_edges(a, L):
    edges = [(a[i], a[i + 1]) for i in range(L)]
    return edges, max(a) + 1


def main():
    KFULL = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    KSAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    NSAMPLE = int(sys.argv[3]) if len(sys.argv) > 3 else 200000

    for k in range(3, KFULL + 1):
        L = 2 * k
        dis = 0; cnt_fast = 0; cnt_svd = 0; tot = 0
        for a in rgs_iter(L + 1):
            tot += 1
            edges, nb = build_edges(a, L)
            has_loop = any(u == v for (u, v) in edges)
            f = fast_is_even_series(a, L)
            s = (not has_loop) and is_even_series_svd(edges, nb)
            cnt_fast += f; cnt_svd += s
            if f != s:
                dis += 1
                if dis <= 5:
                    print(f"  DISAGREE k={k} rgs={a} fast={f} svd={s}")
        print(f"k={k} EXHAUSTIVE: fast={cnt_fast} svd={cnt_svd} disagreements={dis} (of {tot})")
        sys.stdout.flush()

    # sample at KSAMPLE
    k = KSAMPLE
    L = 2 * k
    import random
    random.seed(12345)
    dis = 0
    nb_pos = L + 1
    seen = 0
    for _ in range(NSAMPLE):
        # random restricted growth string
        a = [0] * nb_pos
        mx = 0
        for i in range(1, nb_pos):
            a[i] = random.randint(0, mx + 1)
            if a[i] > mx: mx = a[i]
        edges, nb = build_edges(a, L)
        has_loop = any(u == v for (u, v) in edges)
        f = fast_is_even_series(a, L)
        s = (not has_loop) and is_even_series_svd(edges, nb)
        seen += 1
        if f != s:
            dis += 1
            if dis <= 10:
                print(f"  DISAGREE k={k} rgs={a} fast={f} svd={s}")
    print(f"k={k} SAMPLE({seen}): disagreements={dis}")
    print("If all disagreements==0, the fast integer test is EQUIVALENT to SVD,")
    print("so the k=6 count 2351 is correct and A215257(->2345) mismatches at k=6.")


if __name__ == "__main__":
    main()
