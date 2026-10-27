"""vt15_noncayley.py -- mac-mini-2026-10-27-potentials

Vertex-transitive tournaments of order 15 that are not Cayley tournaments of Z15 must be invariant under one of the two
transitive odd-order groups without a regular subgroup (THM-4812, VT lemma at order 15):
  G_A = A_3 x| Z_5, A_3 = {k in F_3^5 : sum k = 0}, acting on (i, j) in Z_5 x Z_3 by (i, j) -> (i + 1, j) and (i, j) -> (i, j + k_i);
  G_B = A_5 x| Z_3, A_5 = {k in F_5^3 : sum k = 0}, acting on (i, j) in Z_3 x Z_5 likewise.
This script enumerates every tournament invariant under G_A or G_B (orientation chosen per pair of mutually inverse
orbitals), and reports H, |Aut| and tc = H/|Aut| (exact Held-Karp; automorphisms by backtracking with degree pruning).
"""
import itertools


def ham(adj, n):
    full = (1 << n) - 1
    out = [sum(1 << w for w in range(n) if adj[v][w]) for v in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for m in range(1, full + 1):
        row = dp[m]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            o = out[v] & ~m
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                dp[m | (1 << w)][w] += c
    return sum(dp[full])


def aut_count(adj, n):
    outs = [frozenset(w for w in range(n) if adj[v][w]) for v in range(n)]
    # refine by score of out-neighbourhood multiset
    sig = [tuple(sorted(len(outs[w]) for w in outs[v])) for v in range(n)]
    cnt = 0
    perm = [-1] * n
    used = [False] * n

    def bt(i):
        nonlocal cnt
        if i == n:
            cnt += 1
            return
        for x in range(n):
            if used[x] or sig[x] != sig[i]:
                continue
            if all(adj[i][j] == adj[x][perm[j]] for j in range(i)):
                perm[i] = x
                used[x] = True
                bt(i + 1)
                used[x] = False
        perm[i] = -1
    bt(0)
    return cnt


def group_elements(blocks, size, p_block):
    # points (i, j), i in Z_blocks, j in Z_size; generators: shift i -> i+1; translations by sum-zero vectors k in F_size^blocks
    pts = [(i, j) for i in range(blocks) for j in range(size)]
    idx = {pt: t for t, pt in enumerate(pts)}
    gens = []
    gens.append([idx[((i + 1) % blocks, j)] for (i, j) in pts])
    for c in range(blocks - 1):
        k = [0] * blocks
        k[c] = 1
        k[c + 1] = size - 1
        gens.append([idx[(i, (j + k[i]) % size)] for (i, j) in pts])
    # closure
    n = len(pts)
    ident = tuple(range(n))
    G = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                gh = tuple(h[g[x]] for x in range(n))
                if gh not in G:
                    G.add(gh)
                    new.append(gh)
        frontier = new
    return n, list(G)


def invariant_tournaments(n, G):
    # orbitals on ordered pairs
    seen = {}
    orbitals = []
    for u in range(n):
        for v in range(n):
            if u == v or (u, v) in seen:
                continue
            orb = {(g[u], g[v]) for g in G}
            for pr in orb:
                seen[pr] = len(orbitals)
            orbitals.append(orb)
    pairs = []
    done = set()
    for t, orb in enumerate(orbitals):
        if t in done:
            continue
        u, v = next(iter(orb))
        tinv = seen[(v, u)]
        assert tinv != t, "self-paired orbital: no invariant tournament"
        done.add(t)
        done.add(tinv)
        pairs.append((t, tinv))
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for (t, tinv), b in zip(pairs, bits):
            for (u, v) in orbitals[t if b == 0 else tinv]:
                adj[u][v] = 1
        yield adj


def main():
    for name, blocks, size in (("G_A = A_3 x| Z_5 (five blocks of size 3)", 5, 3), ("G_B = A_5 x| Z_3 (three blocks of size 5)", 3, 5)):
        n, G = group_elements(blocks, size, None)
        print(f"{name}: |G| = {len(G)}")
        mintc = None
        for adj in invariant_tournaments(n, G):
            s = [sum(r) for r in adj]
            assert all(x == (n - 1) // 2 for x in s)
            H = ham(adj, n)
            A = aut_count(adj, n)
            assert H % A == 0
            tc = H // A
            mintc = tc if mintc is None else min(mintc, tc)
            print(f"   invariant tournament: H = {H}, |Aut| = {A}, tc = {tc}")
        print(f"   minimum tc = {mintc}")


if __name__ == "__main__":
    main()
