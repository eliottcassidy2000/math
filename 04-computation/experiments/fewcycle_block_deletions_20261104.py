"""Blocks of opus-S75's core theory (THM-4869 B4-B5) and their uniquely Hamiltonian deletion data.

A block (W, A) of type k: W strong, A a subset of W with |A| >= 2 and sum_{x in A} end_x(W) = k, every vertex of W - A has
out-degree >= 2 in W, and no non-empty Z inside W has (#arcs Z -> W - Z) + |Z & A| = 1.

For every block of type k <= 5 and order <= 6 (all labelled tournaments, deduplicated up to isomorphism of (W, A)) this
prints:
  delta(v) = sum_{x in A - v} end_x(W - v) for v in W  (= hc(T - v) for every T with this block, THM-4874 case (v));
  h0(W, A) = #{Hamiltonian paths of W from W - A to A}  (= hc(T - l) when the chain is empty and T - D = {l});
and the predicted numbers of v with hc(T - v) = 1 and = 2.
Re-derives THM-4869 B5's list (k = 2: none; k = 3: 2 blocks; k = 4: 2; k = 5: 4) as a sanity check.
kind-pasteur-2026-11-04-S10."""
from itertools import permutations
from collections import defaultdict


def end_counts(n, rows):
    """end[x] = number of Hamiltonian paths of the tournament (rows = out-neighbour bitmasks) ending at x."""
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


def path_counts_from_to(n, rows):
    """P[s][t] = number of Hamiltonian paths from s to t."""
    full = (1 << n) - 1
    P = [[0] * n for _ in range(n)]
    for s in range(n):
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
        for t in range(n):
            P[s][t] = dp[full][t] if t != s else (1 if n == 1 else 0)
    return P


def strong(n, rows):
    full = (1 << n) - 1
    for s in range(n):
        seen = 1 << s
        stack = [s]
        while stack:
            v = stack.pop()
            nb = rows[v] & ~seen & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                seen |= 1 << w
                stack.append(w)
        if seen != full:
            return False
    return True


def delete(n, rows, v):
    keep = [x for x in range(n) if x != v]
    idx = {x: i for i, x in enumerate(keep)}
    return keep, [sum(1 << idx[w] for w in keep if (rows[x] >> w) & 1) for x in keep]


def is_block(n, rows, A):
    full = (1 << n) - 1
    for v in range(n):
        if not (A >> v) & 1 and bin(rows[v]).count("1") < 2:
            return False
    for Z in range(1, full + 1):
        leave = sum(bin(rows[v] & ~Z & full).count("1") for v in range(n) if (Z >> v) & 1)
        if leave + bin(Z & A).count("1") == 1:
            return False
    return True


def canon(n, rows, A):
    best = None
    for p in permutations(range(n)):
        # p[i] = new label of old vertex i
        inv = [0] * n
        for i, x in enumerate(p):
            inv[x] = i
        code = tuple(sum(1 << p[w] for w in range(n) if (rows[inv[j]] >> w) & 1) for j in range(n))
        a = sum(1 << p[x] for x in range(n) if (A >> x) & 1)
        key = (code, a)
        if best is None or key < best:
            best = key
    return best


def name(n, rows):
    sc = sorted(bin(r).count("1") for r in rows)
    if n == 3:
        return "C3"
    if n == 4:
        return "W4"
    return f"order {n}, scores {sc}"


def main(kmax=5, nmax=6):
    blocks = {}
    for n in range(3, nmax + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for bits in range(1 << len(pairs)):
            rows = [0] * n
            for b, (i, j) in enumerate(pairs):
                if (bits >> b) & 1:
                    rows[i] |= 1 << j
                else:
                    rows[j] |= 1 << i
            if not strong(n, rows):
                continue
            e = end_counts(n, rows)
            if sorted(e)[0] + sorted(e)[1] > kmax:
                continue
            for A in range(1, 1 << n):
                if bin(A).count("1") < 2:
                    continue
                k = sum(e[x] for x in range(n) if (A >> x) & 1)
                if k > kmax or not is_block(n, rows, A):
                    continue
                key = canon(n, rows, A)
                if key not in blocks:
                    blocks[key] = (n, rows, A, k)
    by_type = defaultdict(list)
    for key, (n, rows, A, k) in blocks.items():
        by_type[k].append((n, rows, A))
    for k in range(2, kmax + 1):
        print(f"type {k}: {len(by_type[k])} blocks")
        for n, rows, A in sorted(by_type[k], key=lambda t: (t[0], t[2])):
            e = end_counts(n, rows)
            deltas = []
            for v in range(n):
                keep, r2 = delete(n, rows, v)
                e2 = end_counts(n - 1, r2)
                deltas.append(sum(e2[i] for i, x in enumerate(keep) if (A >> x) & 1))
            P = path_counts_from_to(n, rows)
            h0 = sum(P[s][t] for s in range(n) for t in range(n) if not (A >> s) & 1 and (A >> t) & 1)
            hcW = sum(P[s][t] for s in range(n) for t in range(n) if s != t and (rows[t] >> s) & 1) // n if n >= 3 else 0
            sole_exits = [v for v in range(n) if bin(rows[v]).count("1") == 1]
            print(f"  W = {name(n, rows)}, rows {[bin(r)[2:].zfill(n)[::-1] for r in rows]}, A = {[v for v in range(n) if (A >> v) & 1]}"
                  f" (A = W: {A == (1 << n) - 1}), end counts {e}, hc(W) = {hcW}, sole exits {sole_exits}")
            print(f"    delta(v) for v = 0..{n - 1}: {deltas}  ->  #(hc(T-v) = 1) = {deltas.count(1)}, #(hc(T-v) = 2) = {deltas.count(2)};"
                  f"  h0 = {h0}" + ("  (chain empty allowed)" if A != (1 << n) - 1 else "  (chain must be non-empty)"))


if __name__ == "__main__":
    main()
