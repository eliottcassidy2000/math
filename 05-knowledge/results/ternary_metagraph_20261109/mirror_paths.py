"""Mirror-path decomposition of H mod 4 on blue tilings at odd order n = 2c + 1 (0-based vertices 0..n-1, path i -> i+1).

A blue (sigma-fixed) tiling t gives a tournament T with the anti-automorphism iota(i) = n - 1 - i, fixing only c.
Odd cycles are paired by C -> iota(C) reversed; the fixed ones are the symmetric cycles
    c -> v_1 -> ... -> v_h = u -> iota(u) -> iota(v_(h-1)) -> ... -> iota(v_1) -> c,
with {v_1..v_h} disjoint from its mirror image ("mirror-free" path).  Hence
    alpha_1(T) = sum over u != c with u -> iota(u) of N(u)   (mod 2),
where N(u) is the number of mirror-free directed paths from c to u.  N(u) never uses a mirror arc {i, n-1-i}, so with
x_i = [tile (i, n-1-i) reversed] and y the remaining blue coordinates,
    f(t) = alpha_1 mod 2 = A(y) + sum_i x_i B_i(y),   A(y) = sum_{i<c} N(i),   B_i(y) = N(i) + N(n-1-i)   (mod 2).
This script checks the decomposition against H directly and tests whether B(y) is never 0 and never equal to B(y + 1).
usage: python mirror_paths.py n
"""
import sys
import itertools
import importlib.util
import os

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("p", os.path.join(here, "blue_flip_probe.py"))
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def adjacency(n, tiles, bits):
    out = [0] * n
    for i in range(n - 1):
        out[i] |= 1 << (i + 1)
    for k, (i, j) in enumerate(tiles):
        if bits >> k & 1:
            out[j] |= 1 << i
        else:
            out[i] |= 1 << j
    return out


def mirror_free_counts(n, out):
    c = n // 2
    N = [0] * n
    # DFS over mirror-free paths from c
    def dfs(v, used):
        for u in range(n):
            if out[v] >> u & 1 and not (used >> u & 1) and not (used >> (n - 1 - u) & 1) and u != c:
                N[u] += 1
                dfs(u, used | (1 << u))
    dfs(c, 1 << c)
    return N


def main():
    n = int(sys.argv[1]); assert n % 2 == 1
    c = n // 2
    tiles, orbits = p.setup(n)
    m = len(tiles)
    allones = (1 << m) - 1
    idx = {t: k for k, t in enumerate(tiles)}
    mirror = [idx[(i, n - 1 - i)] for i in range(c)]
    pairs = [o for o in orbits if len(o) == 2]
    masks = [sum(1 << k for k in o) for o in pairs]
    ok = 0; bad = 0; B_zero = 0; B_equal = 0; ny = 0
    for code in range(1 << len(pairs)):
        y = 0
        for q in range(len(pairs)):
            if code >> q & 1:
                y |= masks[q]
        # B(y) from x = 0 tiling
        Nt = mirror_free_counts(n, adjacency(n, tiles, y))
        Ny = mirror_free_counts(n, adjacency(n, tiles, y ^ allones ^ sum(1 << k for k in mirror)))
        B = tuple((Nt[i] + Nt[n - 1 - i]) % 2 for i in range(c))
        Bp = tuple((Ny[i] + Ny[n - 1 - i]) % 2 for i in range(c))
        ny += 1
        if not any(B):
            B_zero += 1
        if B == Bp:
            B_equal += 1
        # verify the decomposition against H for every x
        A = sum(Nt[i] for i in range(c)) % 2
        for xs in range(1 << c):
            t = y
            for i in range(c):
                if xs >> i & 1:
                    t |= 1 << mirror[i]
            f = (A + sum((xs >> i & 1) * B[i] for i in range(c))) % 2
            h = p.H(n, tiles, t)
            if ((h - 1) // 2) % 2 == f:
                ok += 1
            else:
                bad += 1
    print(f"n = {n}: decomposition f = A(y) + x.B(y) checked on {ok + bad} blue tilings: {ok} agree, {bad} disagree")
    print(f"  y-classes {ny}: B(y) = 0 for {B_zero}; B(y) = B(y + 1) for {B_equal}")


if __name__ == "__main__":
    main()
