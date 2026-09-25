#!/usr/bin/env python3
"""Independent audit of duck_decoder section 4: for an even tournament T, sum over
F2 the perfect matchings {p0p1, p2p3, ...} of all Hamiltonian paths p.  Claimed:
order 6 ranks 2/4/6 with counts 1680/17520/13568; rank-2 kernels are K_{1,5}
(960) or K_{3,3} (720); orders 2 and 4 always full rank; hostile: edges larger ->
smaller except 0->4 gives 9 paths, rank 4, null vectors e0+e4, e2+e5.
Also: rank-2 symmetric zero-diagonal F2 matrices with odd degrees <=> K_{a,b}, a,b odd.
Own implementation: enumerate permutations, for each permutation the 2^(C(n,2)-(n-1))
tournaments containing it as a directed path.
"""
from itertools import permutations, combinations


def rank_f2(rows):
    rows = [r for r in rows]
    rank = 0
    ncols = max((r.bit_length() for r in rows), default=0)
    for col in range(ncols):
        piv = None
        for i in range(rank, len(rows)):
            if (rows[i] >> col) & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> col) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def census(n):
    pairs = list(combinations(range(n), 2))
    idx = {p: i for i, p in enumerate(pairs)}
    E = len(pairs)
    # tournament encoded by bitmask: bit i = 1 means pairs[i]=(a,b) oriented a->b (a<b)
    W = [0] * (1 << E)  # kernel as bitmask over pairs
    H = [0] * (1 << E)
    for p in permutations(range(n)):
        fixed_mask = 0
        fixed_val = 0
        for a, b in zip(p, p[1:]):
            if a < b:
                i = idx[(a, b)]
                fixed_val |= 1 << i
            else:
                i = idx[(b, a)]
            fixed_mask |= 1 << i
        match = 0
        for k in range(0, n, 2):
            a, b = p[k], p[k + 1]
            match |= 1 << idx[(min(a, b), max(a, b))]
        free = [i for i in range(E) if not (fixed_mask >> i) & 1]
        # iterate all completions
        for sub in range(1 << len(free)):
            t = fixed_val
            for j, i in enumerate(free):
                if (sub >> j) & 1:
                    t |= 1 << i
            W[t] ^= match
            H[t] += 1
    return pairs, W, H


def matrix_rows(n, pairs, w):
    rows = [0] * n
    for i, (a, b) in enumerate(pairs):
        if (w >> i) & 1:
            rows[a] |= 1 << b
            rows[b] |= 1 << a
    return rows


def graph_type(n, rows):
    # complete bipartite check
    deg = [bin(r).count("1") for r in rows]
    # find parts: vertices with identical rows form independent classes
    classes = {}
    for v in range(n):
        classes.setdefault(rows[v], []).append(v)
    sizes = sorted(len(c) for c in classes.values())
    return tuple(sizes), tuple(sorted(deg))


def main():
    for n in (2, 4):
        pairs, W, H = census(n)
        ranks = {}
        for t in range(1 << len(pairs)):
            r = rank_f2(matrix_rows(n, pairs, W[t]))
            ranks[r] = ranks.get(r, 0) + 1
            assert H[t] % 2 == 1
        print("order", n, "rank distribution", ranks)
    n = 6
    pairs, W, H = census(n)
    ranks = {}
    types = {}
    odddeg = True
    for t in range(1 << len(pairs)):
        rows = matrix_rows(n, pairs, W[t])
        if any(bin(r).count("1") % 2 == 0 for r in rows):
            odddeg = False
        r = rank_f2(rows)
        ranks[r] = ranks.get(r, 0) + 1
        if r == 2:
            gt = graph_type(n, rows)
            types[gt] = types.get(gt, 0) + 1
    print("order 6 rank distribution", dict(sorted(ranks.items())), "total", sum(ranks.values()))
    print("order 6 all kernel degrees odd:", odddeg)
    print("order 6 rank-2 kernel types (identical-row class sizes, degree sequence):", types)
    # hostile: every edge larger -> smaller except 0->4
    t = 0
    for i, (a, b) in enumerate(pairs):
        # bit=1 means a->b with a<b ; larger->smaller means b->a => bit 0
        if (a, b) == (0, 4):
            t |= 1 << i
    rows = matrix_rows(n, pairs, W[t])
    print("hostile: H =", H[t], "rank =", rank_f2(rows))
    # null vectors
    nulls = []
    for v in range(1, 1 << n):
        ok = True
        for r in rows:
            if bin(r & v).count("1") % 2:
                ok = False
                break
        if ok:
            nulls.append([i for i in range(n) if (v >> i) & 1])
    print("hostile null space vectors:", nulls)


if __name__ == "__main__":
    main()
