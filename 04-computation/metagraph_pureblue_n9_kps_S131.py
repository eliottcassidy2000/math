#!/usr/bin/env python3
r"""
metagraph_pureblue_n9_kps_S131.py   (kind-pasteur-2026-07-26-S131, HYP-4997)

Decide pure-blue(9) for the conjecture

    pure-blue(n) = floor((n+1)/2) - [n even]     (HYP-4997 / kps-S66 PB1)

verified at n=3..8 (2,1,3,2,4,3); prediction pure-blue(9) = 5.

Method (kps-S66 blue-sub-cube trick, upgraded canonicalizer):
pure-blue classes are SC classes whose ENTIRE tiling fiber is
grid-symmetric, i.e. blue-mult == tc = H/|Aut| (LEM-003).  Only the
blue sub-cube Fix(sigma) (2^e tilings, e = floor((n-1)^2/4)) can
contain blue tilings, so enumerate it, group by isomorphism class,
and compare blue multiplicity with the class tiling count.

n=9: e=16 -> 65,536 grid-symmetric tilings.  Exact isomorphism via
iterated WL refinement fingerprint buckets + backtracking iso tests
inside buckets; |Aut| by colour-pruned backtracking; H by Ham-path DP.

Controls: n=5,6,7 reproduce 3,2,4 (gn_lines_parity_census, kps-S66,
codex atlases) and n=8 must reproduce 3 (THM-790 S7 / opus-S305
metagraph_flow_n8_check).  purity_n8.out's PB(8)=0 is the known
refuted artifact (tile-ordering bug) and is NOT a control.

Universe: all grid-symmetric tilings at each n in {5,...,9}; exact
integer arithmetic throughout; no sampling.
"""
from itertools import permutations
import sys

def tiles_of(n):
    T = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2:
                T.append((x, y))
    return T

def grid_orbits(n, T):
    pos = {t: i for i, t in enumerate(T)}
    seen = set(); orbits = []
    for i, t in enumerate(T):
        if i in seen:
            continue
        x, y = t
        j = pos[(n - y + 1, n - x + 1)]
        orb = {i, j}; seen |= orb; orbits.append(sorted(orb))
    return orbits

def tour_from_bits(n, T, bits):
    """adjacency matrix, 0-indexed vertices 0..n-1 for speed; base path
    (k+1)->k in 1-indexed becomes k -> k-1."""
    A = [[0] * n for _ in range(n)]
    for k in range(1, n):
        A[k][k - 1] = 1
    for (x, y), b in zip(T, bits):
        if b == 0:
            A[x - 1][y - 1] = 1
        else:
            A[y - 1][x - 1] = 1
    return A

def wl_colors(n, A, rounds=3):
    col = [sum(A[i]) for i in range(n)]           # scores
    for _ in range(rounds):
        nxt = []
        for i in range(n):
            outs = tuple(sorted(col[j] for j in range(n) if A[i][j]))
            ins = tuple(sorted(col[j] for j in range(n) if j != i and not A[i][j]))
            nxt.append(hash((col[i], outs, ins)))
        # compress
        uniq = {c: k for k, c in enumerate(sorted(set(nxt)))}
        col = [uniq[c] for c in nxt]
    return col

def fingerprint(n, A):
    col = wl_colors(n, A)
    per = sorted(zip(col, [sum(A[i]) for i in range(n)]))
    return tuple(per)

def iso_map(n, A, B):
    """return True iff A ~ B (tournament isomorphism), colour-pruned."""
    ca, cb = wl_colors(n, A), wl_colors(n, B)
    if sorted(ca) != sorted(cb):
        return False
    # candidate targets per vertex by colour
    cand = [[j for j in range(n) if cb[j] == ca[i]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i]))
    mapping = [-1] * n
    used = [False] * n

    def bt(k):
        if k == n:
            return True
        i = order[k]
        for j in cand[i]:
            if used[j]:
                continue
            ok = True
            for k2 in range(k):
                i2 = order[k2]
                if A[i][i2] != B[j][mapping[i2]] or A[i2][i] != B[mapping[i2]][j]:
                    ok = False
                    break
            if ok:
                mapping[i] = j; used[j] = True
                if bt(k + 1):
                    return True
                used[j] = False; mapping[i] = -1
        return False

    return bt(0)

def aut_size(n, A):
    col = wl_colors(n, A)
    cand = [[j for j in range(n) if col[j] == col[i]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i]))
    mapping = [-1] * n
    used = [False] * n
    count = 0

    def bt(k):
        nonlocal count
        if k == n:
            count += 1
            return
        i = order[k]
        for j in cand[i]:
            if used[j]:
                continue
            ok = True
            for k2 in range(k):
                i2 = order[k2]
                if A[i][i2] != A[j][mapping[i2]] or A[i2][i] != A[mapping[i2]][j]:
                    ok = False
                    break
            if ok:
                mapping[i] = j; used[j] = True
                bt(k + 1)
                used[j] = False; mapping[i] = -1

    bt(0)
    return count

def ham(n, A):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            Al = A[last]
            for nx in range(n):
                if Al[nx] and not (mask >> nx) & 1:
                    dp[mask | (1 << nx)][nx] += c
    return sum(dp[full][i] for i in range(n))

def run(n):
    T = tiles_of(n); m = len(T); orbits = grid_orbits(n, T)
    e = len(orbits)
    assert e == ((n - 1) ** 2) // 4, (n, e)
    buckets = {}          # fingerprint -> list of (rep_matrix, class_index)
    classes = []          # class_index -> [rep, blue_mult]
    for code in range(1 << e):
        bits = [0] * m
        for oi, orb in enumerate(orbits):
            v = (code >> oi) & 1
            for idx in orb:
                bits[idx] = v
        A = tour_from_bits(n, T, bits)
        fp = fingerprint(n, A)
        hit = None
        for rep, ci in buckets.setdefault(fp, []):
            if iso_map(n, A, rep):
                hit = ci
                break
        if hit is None:
            classes.append([A, 0])
            buckets[fp].append((A, len(classes) - 1))
            hit = len(classes) - 1
        classes[hit][1] += 1
    # decide pure-blue per touched class
    pure = []
    total_blue = 0
    for A, bm in classes:
        total_blue += bm
        H = ham(n, A); au = aut_size(n, A)
        assert H % au == 0, (H, au)
        tc = H // au
        if bm == tc:
            pure.append((H, au, tc))
    assert total_blue == 1 << e
    pred = ((n + 1) // 2) - (1 if n % 2 == 0 else 0)
    print(f"n={n}: blue cube 2^{e}; classes touched={len(classes)}; "
          f"PURE-BLUE={len(pure)} predicted={pred} match={len(pure) == pred}")
    for H, au, tc in sorted(pure):
        print(f"    pure-blue class: H={H} |Aut|={au} tc={tc}")
    return len(pure), pred

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [5, 6, 7, 8, 9]
    ok = True
    for n in ns:
        got, pred = run(n)
        if n <= 8:
            expected = {3: 2, 4: 1, 5: 3, 6: 2, 7: 4, 8: 3}[n]
            if got != expected:
                print(f"CONTROL FAILURE at n={n}: got {got}, canon {expected}")
                ok = False
    print("ALL CONTROLS PASSED" if ok else "CONTROL FAILURES PRESENT")
