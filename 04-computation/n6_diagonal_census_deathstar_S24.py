#!/usr/bin/env python3
"""death-star-2026-07-16-S24 (HYP-7051): THE n=6 DIAGONAL CENSUS.
d(W) = diagonal complexity (min degree with delta functions realizable by F2 polynomials
on W, i.e. evaluation-matrix rank = |W|) for all 56 iso classes of 6-tournaments in the
tiling cube F2^10; internal wiggly degrees; Lemma 2.3 (2607.14068) bound tightness at the
larger regime; correlations with |W|, H, c3, |Aut|, SC/NS; counting-bound exceptionality."""
from itertools import permutations, combinations
from math import comb
import sys, time

N = 6
TILES = [(x, y) for y in range(1, N - 1) for x in range(N, y + 1, -1) if x - y >= 2]
M = len(TILES)  # 10
PERMS = list(permutations(range(1, N + 1)))

def tourn(bits_int):
    T = [[0] * (N + 1) for _ in range(N + 1)]
    for a in range(N, 1, -1):
        T[a][a - 1] = 1
    for idx, (a, b) in enumerate(TILES):
        if (bits_int >> idx) & 1:
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def key_of(T, perm):
    s = 0
    for i in range(1, N + 1):
        pi = perm[i - 1]
        for j in range(1, N + 1):
            if i != j:
                s = (s << 1) | T[pi][perm[j - 1]]
    return s

def canon_and_aut(T):
    best = None; aut = 0
    ident = key_of(T, tuple(range(1, N + 1)))
    for perm in PERMS:
        s = key_of(T, perm)
        if best is None or s < best:
            best = s
        if s == ident:
            aut += 1
    return best, aut

def canon_only(T):
    best = None
    for perm in PERMS:
        s = key_of(T, perm)
        if best is None or s < best:
            best = s
    return best

def H_count(T):
    # Hamiltonian paths via DP over subsets
    full = (1 << N) - 1
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1 << N):
        for v in range(N):
            if dp[mask][v]:
                for w in range(N):
                    if not (mask >> w) & 1 and T[v + 1][w + 1]:
                        dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(N))

def c3_of(T):
    outs = [sum(T[i][j] for j in range(1, N + 1)) for i in range(1, N + 1)]
    return comb(N, 3) - sum(comb(o, 2) for o in outs)

def rev(T):
    R = [[0] * (N + 1) for _ in range(N + 1)]
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            R[i][j] = T[j][i]
    return R

def diag_degree(W, m=M):
    for d in range(0, m + 1):
        rows = []
        for dd in range(0, d + 1):
            for sup in combinations(range(m), dd):
                row = 0
                for i, w in enumerate(W):
                    val = 1
                    for b in sup:
                        val &= (w >> b) & 1
                    row |= val << i
                rows.append(row)
        basis = []
        for r in rows:
            x = r
            for b in basis:
                hb = b.bit_length() - 1
                if (x >> hb) & 1:
                    x ^= b
            if x:
                basis.append(x); basis.sort(key=lambda z: -z.bit_length())
        if len(basis) == len(W):
            return d
    return m

def counting_lb(sz, m=M):
    tot = 0
    for d in range(m + 1):
        tot += comb(m, d)
        if tot >= sz:
            return d
    return m

if __name__ == "__main__":
    t0 = time.time()
    classes = {}
    for bits in range(1 << M):
        T = tourn(bits)
        key = canon_only(T)
        classes.setdefault(key, []).append(bits)
        if bits % 128 == 0:
            sys.stderr.write(f"\r{bits}/1024 [{time.time()-t0:.0f}s]")
    sys.stderr.write("\n")
    print(f"#classes = {len(classes)} (A000568(6) = 56 expected)  [{time.time()-t0:.0f}s]")
    rows = []
    for key, W in classes.items():
        T = tourn(W[0])
        H = H_count(T); c3 = c3_of(T)
        _, aut = canon_and_aut(T)
        sc = (canon_only(rev(T)) == key)
        Wset = set(W)
        ie = sum(1 for w in W for b in range(M) if (w ^ (1 << b)) in Wset) // 2
        avg = 2 * ie / len(W)
        dW = diag_degree(W)
        clb = counting_lb(len(W))
        rows.append(dict(sz=len(W), H=H, c3=c3, aut=aut, sc=sc, ie=ie,
                         avg=avg, dW=dW, clb=clb))
    rows.sort(key=lambda r: (r['sz'], r['H']))
    print(f"{'|W|':>4} {'H':>4} {'c3':>3} {'Aut':>4} {'SC':>3} {'intE':>5} {'avgdeg':>7} "
          f"{'d(W)':>5} {'cntLB':>5} {'exc':>4} {'2d-vs-avg':>9}")
    for r in rows:
        exc = r['dW'] - r['clb']
        print(f"{r['sz']:>4} {r['H']:>4} {r['c3']:>3} {r['aut']:>4} {str(r['sc'])[:1]:>3} "
              f"{r['ie']:>5} {r['avg']:>7.3f} {r['dW']:>5} {r['clb']:>5} {exc:>4} "
              f"{2*r['dW'] - r['avg']:>9.3f}")
    # summaries
    n_exc = sum(1 for r in rows if r['dW'] > r['clb'])
    n_tight = sum(1 for r in rows if abs(2 * r['dW'] - r['avg']) < 1e-9 and r['sz'] > 1)
    zero_ie = sum(1 for r in rows if r['ie'] == 0)
    print(f"\nSUMMARY: {len(rows)} classes; counting-LB exceeded by {n_exc}; "
          f"Lemma-2.3 tight (nontrivial): {n_tight}; zero-internal-edge classes: {zero_ie}")
    # equal-size separations by d(W)
    from collections import defaultdict
    bysz = defaultdict(set)
    for r in rows:
        bysz[r['sz']].add(r['dW'])
    seps = {sz: ds for sz, ds in bysz.items() if len(ds) > 1}
    print(f"sizes where d(W) separates equal-size classes: {seps}")
    # SC vs NS d(W) distribution
    for scv in [True, False]:
        ds = [r['dW'] for r in rows if r['sc'] == scv]
        if ds:
            print(f"SC={scv}: n={len(ds)} d(W) values {sorted(set(ds))} mean={sum(ds)/len(ds):.2f}")
    print(f"[total {time.time()-t0:.1f}s]")
