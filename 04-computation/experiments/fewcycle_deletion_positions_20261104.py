"""Direct check of THM-4874's case table, vertex by vertex, on explicitly built one-arc joins T = Q[z <- (D, u)].

Prefixes: every labelled rooted uniquely Hamiltonian (Q, z) of order 3-5 with z a sole exit of Q (the head l is z's
only out-neighbour), plus the one-vertex prefix.
Cores: the eight blocks (W, A) of type <= 5 (THM-4869 B5, as re-derived by fewcycle_block_deletions_20261104.py),
chains of length p = 0..PMAX (p = 0 only when A != W), with every choice of the free arcs between u and w_2..w_p.
For every vertex v of every T, hc(T - v) is computed by subset DP and compared with the prediction by position:
  v in R - l: k * s_l(R - v);  v = u: 0;  v = l: 0 if |R| >= 2, else k (p >= 1) or h0(W, A) (p = 0);
  v = w_i: k;  v in W: delta(v) = sum_{x in A - v} end_x(W - v).
Also checks hc(T) = k and that T has no vertex of out-degree 1.
kind-pasteur-2026-11-04-S10."""
from itertools import permutations, product


def hc(rows):
    n = len(rows)
    if n < 3:
        return 0
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    dp[1][0] = 1
    for S in range(1, 1 << n, 2):
        for v in range(n):
            c = dp[S][v]
            if not c:
                continue
            nb = rows[v] & ~S & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                dp[S | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(1, n) if rows[v] & 1)


def paths_from(rows, s):
    """number of Hamiltonian paths starting at s (1 if the tournament is the single vertex s)"""
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    dp[1 << s][s] = 1
    for S in range(1, 1 << n):
        if not (S >> s) & 1:
            continue
        for v in range(n):
            c = dp[S][v]
            if not c:
                continue
            nb = rows[v] & ~S & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                dp[S | (1 << w)][w] += c
    return sum(dp[full])


def ends(rows):
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c:
                continue
            nb = rows[v] & ~S & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                dp[S | (1 << w)][w] += c
    return dp[full]


def induced(rows, keep):
    idx = {x: i for i, x in enumerate(keep)}
    return [sum(1 << idx[w] for w in keep if (rows[x] >> w) & 1) for x in keep]


def from_arcs(n, arcs):
    rows = [0] * n
    for a, b in arcs:
        rows[a] |= 1 << b
    return rows


def labelled_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield rows


# the eight blocks, in the labels printed by fewcycle_block_deletions_20261104.out (row strings: bit j of vertex i)
def rows_from_strings(strs):
    return [sum(1 << j for j, ch in enumerate(s) if ch == "1") for s in strs]


C3 = rows_from_strings(["001", "100", "010"])
W4 = rows_from_strings(["0001", "1000", "1100", "0110"])
AT5 = rows_from_strings(["00001", "10000", "11000", "11100", "01110"])
W5 = rows_from_strings(["00011", "10000", "11000", "01100", "01110"])
BLOCKS = [
    ("type 3: (C3, all)", C3, [0, 1, 2]),
    ("type 3: (W4, {c,d})", W4, [0, 1]),
    ("type 4: (W4, {b,c,d})", W4, [0, 1, 2]),
    ("type 4: (W4, {a,c,d})", W4, [0, 1, 3]),
    ("type 5: (W4, all)", W4, [0, 1, 2, 3]),
    ("type 5: (AT5, sole exits)", AT5, [0, 1]),
    ("type 5: (W5, A1)", W5, [0, 1]),
    ("type 5: (W5, A2)", W5, [1, 4]),
]


def build(Qrows, z, l, W, A, p, uarcs):
    """T on vertices: R = Q - z (relabelled 0..r-1, l first), then chain w_1..w_p, then W, then u."""
    if Qrows is None:
        R = ["l"]
        Rrows = [0]
        lr = 0
    else:
        q = len(Qrows)
        keep = [l] + [x for x in range(q) if x not in (z, l)]
        Rrows = induced(Qrows, keep)
        R = keep
        lr = 0
    r = len(Rrows)
    w = len(W)
    n = r + p + w + 1
    chain = list(range(r, r + p))
    Wv = list(range(r + p, r + p + w))
    u = n - 1
    rows = [0] * n
    for i in range(r):
        rows[i] = Rrows[i]
    for i in range(r):                    # R beats all of D, except u -> l
        for y in chain + Wv + [u]:
            if not (i == lr and y == u):
                rows[i] |= 1 << y
    rows[u] |= 1 << lr
    for a in range(p):                     # chain transitive, chain beats W
        for b in range(a + 1, p):
            rows[chain[a]] |= 1 << chain[b]
        for y in Wv:
            rows[chain[a]] |= 1 << y
    for i in range(w):                     # the block
        for j in range(w):
            if (W[i] >> j) & 1:
                rows[Wv[i]] |= 1 << Wv[j]
    for i in range(w):                     # u: beaten by A, beats W - A
        if i in A:
            rows[Wv[i]] |= 1 << u
        else:
            rows[u] |= 1 << Wv[i]
    if p >= 1:
        rows[u] |= 1 << chain[0]
        for t in range(1, p):              # free arcs u - w_(t+1)
            if uarcs[t - 1]:
                rows[u] |= 1 << chain[t]
            else:
                rows[chain[t]] |= 1 << u
    return rows, list(range(r)), lr, chain, Wv, u


def main(pmax=2, qmax=5):
    prefixes = [(None, None, None)]
    seen = set()
    for q in range(3, qmax + 1):
        for Qrows in labelled_tournaments(q):
            if hc(Qrows) != 1:
                continue
            for z in range(q):
                if bin(Qrows[z]).count("1") != 1:
                    continue
                key = min((tuple(sum(1 << pm[w] for w in range(q) if (Qrows[x] >> w) & 1) for x in sorted(range(q), key=lambda x: pm[x])), pm[z])
                          for pm in permutations(range(q)))
                if key in seen:
                    continue
                seen.add(key)
                l = Qrows[z].bit_length() - 1
                prefixes.append((Qrows, z, l))
    print(f"rooted uniquely Hamiltonian prefixes (Q, z) up to isomorphism, orders 3-{qmax}, plus the one-vertex prefix: {len(prefixes)}")
    total = bad = 0
    for name, W, A in BLOCKS:
        k = sum(e for i, e in enumerate(ends(W)) if i in A)
        delta = []
        for v in range(len(W)):
            keep = [x for x in range(len(W)) if x != v]
            e2 = ends(induced(W, keep))
            delta.append(sum(e2[i] for i, x in enumerate(keep) if x in A))
        h0 = 0
        # h0 = Hamiltonian paths of W from W - A to A
        n_w = len(W)
        full = (1 << n_w) - 1
        for s in range(n_w):
            if s in A:
                continue
            dp = [[0] * n_w for _ in range(1 << n_w)]
            dp[1 << s][s] = 1
            for S in range(1, 1 << n_w):
                for v in range(n_w):
                    c = dp[S][v]
                    if not c:
                        continue
                    nb = W[v] & ~S & full
                    while nb:
                        x = (nb & -nb).bit_length() - 1
                        nb &= nb - 1
                        dp[S | (1 << x)][x] += c
            h0 += sum(dp[full][t] for t in A)
        nT = 0
        for p in range(0, pmax + 1):
            if p == 0 and sorted(A) == list(range(len(W))):
                continue
            for uarcs in product([0, 1], repeat=max(p - 1, 0)):
                for Qrows, z, l in prefixes:
                    rows, R, lr, chain, Wv, u = build(Qrows, z, l, W, A, p, uarcs)
                    n = len(rows)
                    assert all(bin(rows[i]).count("1") + sum(1 for j in range(n) if (rows[j] >> i) & 1) == n - 1 for i in range(n))
                    if hc(rows) != k or min(bin(r).count("1") for r in rows) < 2:
                        bad += 1
                        print("  BAD T:", name, p, uarcs, "hc", hc(rows))
                        continue
                    for v in range(n):
                        keep = [x for x in range(n) if x != v]
                        got = hc(induced(rows, keep))
                        if v in R and v != lr:
                            Rk = [x for x in R if x != v]
                            pred = k * paths_from(induced(rows, Rk), Rk.index(lr))
                        elif v == u:
                            pred = 0
                        elif v == lr:
                            pred = 0 if len(R) >= 2 else (k if p >= 1 else h0)
                        elif v in chain:
                            pred = k
                        else:
                            pred = delta[Wv.index(v)]
                        total += 1
                        if got != pred:
                            bad += 1
                            print("  MISMATCH", name, "p", p, "v", v, "got", got, "pred", pred)
                    nT += 1
        print(f"{name}: k = {k}, delta = {delta}, h0 = {h0}; tournaments built: {nT}")
    print(f"vertex deletions checked: {total}; mismatches or bad tournaments: {bad}")


if __name__ == "__main__":
    main()
