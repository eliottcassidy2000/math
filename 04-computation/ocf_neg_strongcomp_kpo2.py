#!/usr/bin/env python3
"""
ocf_neg_strongcomp_kpo2.py -- THREAD C mini-probe 3 (HYP-2380)
================================================================
LEMMA (strong-component factorization of the OCF polynomial; trivial proof:
every directed cycle of a tournament lies inside one strong component, and
cycles in different components are automatically vertex-disjoint):

    I(Omega(T), x)  =  prod_{strong components C of T}  I(Omega(T|_C), x)

Consequences:  sign(I(-2)) = prod of component signs;  H = prod H(C) (known);
positive I(-2) at small n comes from an even number of negative components.

Checks:
 1. Verify the factorization on 800 random n=7 + 200 random n=8 tournaments
    (exact alpha-vector convolution vs direct census).
 2. The tight family S4 (+) S4 (n=8, block-transitive):
    I = (1+2x)^2 = [1,4,4] (discriminant EXACTLY 0, double root -1/2),
    I(-2) = 9, I(-1) = 1, H = 25.  Shows n=8 positives exist (sample missed
    them) and that I(-1) = +1 occurs beyond the transitive class.
 3. S4 (+) C3 (+) pt (n=8): I = (1+2x)(1+x) = [1,3,2], I(-2) = 3.
Pure integers; fresh cycle code (MISTAKE-001/023/054 guards).
"""

import itertools, random

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, plist):
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(plist):
        if (mask >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def directed_odd_cycles(adj, verts=None):
    n = len(adj)
    if verts is None:
        verts = list(range(n))
    out = []
    for k in range(3, len(verts) + 1, 2):
        for S in itertools.combinations(verts, k):
            v0 = S[0]
            rest = S[1:]
            for perm in itertools.permutations(rest):
                if not adj[v0][perm[0]]:
                    continue
                ok = True
                for a in range(k - 2):
                    if not adj[perm[a]][perm[a + 1]]:
                        ok = False
                        break
                if ok and adj[perm[-1]][v0]:
                    m = (1 << v0)
                    for v in perm:
                        m |= (1 << v)
                    out.append((m, k))
    return out

def alpha_from_cycles(cycles, kmax):
    masks = [m for m, _ in cycles]
    alpha = [1, len(masks)]
    if kmax >= 2:
        cnt = {}
        L = len(masks)
        def dfs(start, used, depth):
            if depth >= 2:
                cnt[depth] = cnt.get(depth, 0) + 1
            if depth == kmax:
                return
            for t in range(start, L):
                if not (masks[t] & used):
                    dfs(t + 1, used | masks[t], depth + 1)
        for s in range(L):
            dfs(s + 1, masks[s], 1)
        for k in range(2, kmax + 1):
            alpha.append(cnt.get(k, 0))
    while len(alpha) > 1 and alpha[-1] == 0:
        alpha.pop()
    return alpha

def strong_components(adj):
    """Strong components of a tournament via reachability (n small)."""
    n = len(adj)
    reach = [set([i]) for i in range(n)]
    for i in range(n):
        stack = [i]
        seen = {i}
        while stack:
            u = stack.pop()
            for v in range(n):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    stack.append(v)
        reach[i] = seen
    comps = []
    assigned = [False] * n
    for i in range(n):
        if assigned[i]:
            continue
        comp = [j for j in range(n) if j in reach[i] and i in reach[j]]
        for j in comp:
            assigned[j] = True
        comps.append(sorted(comp))
    return comps

def poly_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def I_eval(alpha, x):
    v = 0
    for k in reversed(range(len(alpha))):
        v = v * x + alpha[k]
    return v

def block_transitive(blocks):
    sizes = [len(b) for b in blocks]
    n = sum(sizes)
    adj = [[0] * n for _ in range(n)]
    off = []
    o = 0
    for s in sizes:
        off.append(o)
        o += s
    for bi, B in enumerate(blocks):
        for i in range(len(B)):
            for j in range(len(B)):
                adj[off[bi] + i][off[bi] + j] = B[i][j]
    for bi in range(len(blocks)):
        for bj in range(bi + 1, len(blocks)):
            for i in range(sizes[bi]):
                for j in range(sizes[bj]):
                    adj[off[bi] + i][off[bj] + j] = 1
    return adj

def held_karp_H(adj):
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full])

C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
S4 = [[0, 1, 1, 0],
      [0, 0, 1, 1],
      [0, 0, 0, 1],
      [1, 0, 0, 0]]
SINGLE = [[0]]

def main():
    rng = random.Random(777002)

    print("=" * 74)
    print("CHECK 1: strong-component factorization of I(Omega(T), x)")
    print("=" * 74)
    bad = 0
    tested = 0
    for n, NS in ((7, 800), (8, 200)):
        plist = pairs_of(n)
        m = len(plist)
        nontriv = 0
        for _ in range(NS):
            mask = rng.randrange(1 << m)
            adj = adj_from_mask(n, mask, plist)
            alpha = alpha_from_cycles(directed_odd_cycles(adj), n // 3)
            comps = strong_components(adj)
            prod = [1]
            for comp in comps:
                ac = alpha_from_cycles(directed_odd_cycles(adj, comp),
                                       max(1, len(comp) // 3))
                prod = poly_mul(prod, ac)
            tested += 1
            if len(comps) > 1:
                nontriv += 1
            if prod != alpha:
                bad += 1
                print(f"  *** FACTORIZATION FAIL n={n} mask={mask}: "
                      f"direct={alpha} product={prod} comps={comps}")
        print(f"n={n}: {NS} random tournaments, "
              f"{nontriv} with >1 strong component -- factorization "
              f"{'ALL OK' if bad == 0 else 'FAILURES!'}")
    print(f"TOTAL: {tested} tested, {bad} failures")

    print()
    print("=" * 74)
    print("CHECK 2/3: n=8 positive-I(-2) constructions (sample missed them)")
    print("=" * 74)
    for name, blocks, exp_alpha in (
            ("S4 (+) S4        [n=8]", [S4, S4], [1, 4, 4]),
            ("S4 (+) C3 (+) pt [n=8]", [S4, C3, SINGLE], [1, 3, 2]),
            ("C3 (+) C3 (+) 2pt[n=8]", [C3, C3, SINGLE, SINGLE], [1, 2, 1])):
        adj = block_transitive(blocks)
        n = len(adj)
        alpha = alpha_from_cycles(directed_odd_cycles(adj), n // 3)
        H = held_karp_H(adj)
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0
        disc = a1 * a1 - 4 * a2
        ok = "OK" if I_eval(alpha, 2) == H else "*** OCF FAIL ***"
        print(f"{name} alpha={alpha} (expect {exp_alpha}) H={H} "
              f"I(2)={I_eval(alpha,2)}[{ok}] I(-2)={I_eval(alpha,-2)} "
              f"I(-1)={I_eval(alpha,-1)} disc={disc}")
    print()
    print("S4(+)S4: I=(1+2x)^2, double root -1/2, disc=0, I(-2)=9>3, I(-1)=+1:")
    print("  positive I(-2) values at n=8 exist and exceed the n=7 max (3);")
    print("  I(-1)=+1 occurs beyond the transitive class (alpha=[1,4,4]).")
    print("done.")

if __name__ == "__main__":
    main()
