#!/usr/bin/env python3
"""Adversarial audit of lane square_sum_hamiltonicity (wave 6, 2026-09-22).

Independent recomputation, written without reusing the lane's Ham class:
  A1  components (union-find) and degree<=1 vertices, n = 1..40; self-loop values
  A2  Hamiltonian path / cycle existence 1..40 with a DFS whose only pruning is
      "every unvisited vertex keeps >= 1 available neighbour" (weaker than the
      lane's three-rule prune, hence an independent oracle); exact undirected
      path counts for 15 <= n <= 33 and cycle counts for 32 <= n <= 36
  A3  cut sets |S| <= 3 for n = 18..24, brute force
  A4  the n = 24 certificate: the forced fragment printed by the lane in DFS
      order is NOT a path (4-24 is not an edge); the path order is checked
  A5  n = 15 census at vertex 4; n = 23 endpoint distribution
  A6  Gerbicz's 25-fold blow-up re-implemented from his PARI pseudo-code
      (F(a,b,c,ty)), independent of the lane's GLUE/blow_up; the partition
      claim [13, 25n+12] checked for m = 1, 2
  A7  NEW (PROVED): Q_n is connected for every n >= 14, because for n >= 5 a
      square lies in [n+1, 2n-1]; the inequality is checked to 200000 and
      proved for n >= 10 by (ceil(sqrt(n+1)))^2 <= (sqrt(n+1)+1)^2 <= 2n-1
  A8  boundary cases of the pasted table (Q_1 is connected)
No timing output; identical under python3 -O (no bare asserts).
"""
import sys
from itertools import combinations
from math import isqrt

sys.setrecursionlimit(10000)


def fail(msg):
    raise RuntimeError(msg)


def sq(m):
    r = isqrt(m)
    return r * r == m


def nbrs(n):
    return {v: [w for w in range(1, n + 1) if w != v and sq(v + w)] for v in range(1, n + 1)}


def uf_components(n, adj):
    parent = list(range(n + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for v in adj:
        for w in adj[v]:
            a, b = find(v), find(w)
            if a != b:
                parent[a] = b
    comps = {}
    for v in range(1, n + 1):
        comps.setdefault(find(v), []).append(v)
    return sorted(sorted(c) for c in comps.values())


def dfs_paths(n, adj, start, count_all, need_cycle=False):
    """Directed Hamiltonian paths from start; weak prune only."""
    nb = {v: sum(1 << w for w in adj[v]) for v in adj}
    full = sum(1 << v for v in range(1, n + 1))
    res = [0, None]
    path = [start]

    def rec(cur, vis):
        if vis == full:
            if need_cycle and not (nb[cur] >> start) & 1:
                return False
            res[0] += 1
            if res[1] is None:
                res[1] = list(path)
            return not count_all
        unv = full & ~vis
        m = unv
        while m:
            b = m & -m
            m ^= b
            v = b.bit_length() - 1
            if nb[v] & (unv | (1 << cur)) == 0:
                return False
        cand = [(bin(nb[v] & unv).count("1"), v) for v in adj[cur] if not (vis >> v) & 1]
        cand.sort()
        for _, v in cand:
            path.append(v)
            if rec(v, vis | (1 << v)):
                return True
            path.pop()
        return False

    rec(start, 1 << start)
    return res


def valid_chain(n, seq):
    return sorted(seq) == list(range(1, n + 1)) and all(sq(seq[i] + seq[i + 1]) for i in range(len(seq) - 1))


# Gerbicz post #19 (Mersenneforum 2018-01-17), transcribed from the archived thread:
#   F(a,b,c,ty): ty=0 append the integer c; ty=1 append T(c); ty=-1 append R(T(c))
GERBICZ_F_CALLS = [(1, 0), (-1, 1), (1, 1), (-7, -1), (6, 1), (-6, 1), (0, -1), (11, 0), (-5, -1),
                   (5, 0), (4, 0), (12, 0), (-12, 1), (12, 1), (7, -1), (-8, 1), (2, -1), (-3, 1),
                   (9, 0), (7, 0), (4, 1), (-4, 1), (10, 0), (6, 0), (5, 1), (-11, -1), (2, 0),
                   (-2, 1), (8, 0), (3, 1), (-9, -1), (9, -1), (-10, 1), (10, 1), (11, 1), (8, -1),
                   (3, 0)]
GERBICZ_V35 = [1, 8, 28, 21, 4, 32, 17, 19, 6, 30, 34, 15, 10, 26, 23, 13, 12, 24, 25, 11, 5, 20,
               29, 35, 14, 2, 7, 18, 31, 33, 16, 9, 27, 22, 3]


def gerbicz_fun_odd(a):
    n = len(a)
    if n % 2 == 0 or a[0] != 1 or a[-1] != 3:
        fail("fun_odd input")
    r = []
    for c, ty in GERBICZ_F_CALLS:
        if ty == 0:
            r.append(c)
        else:
            v = [25 * a[i] + c * (-1) ** i for i in range(n)]   # (-1)^(i+1) with 1-based i
            if ty == -1:
                v = v[::-1]
            r.extend(v)
    return r


def main():
    P = print
    P("# collatz_mod6_20260922_w6_square_sum_hamiltonicity_audit")
    NMAX = 40
    G = {n: nbrs(n) for n in range(1, NMAX + 1)}

    # ---- A1 ----
    P("== A1  components / leaves (union-find, independent)")
    ncomp = {}
    low = {}
    for n in range(1, NMAX + 1):
        cs = uf_components(n, G[n])
        ncomp[n] = len(cs)
        low[n] = [v for v in range(1, n + 1) if len(G[n][v]) <= 1]
        if n in (12, 13, 14):
            P(f"  n={n}: components {cs}")
    P(f"  #components n=1..40: {[ncomp[n] for n in range(1, NMAX + 1)]}")
    if not (all(ncomp[n] == 3 for n in range(4, 13)) and ncomp[13] == 2
            and all(ncomp[n] == 1 for n in range(14, NMAX + 1)) and ncomp[1] == 1
            and ncomp[2] == 2 and ncomp[3] == 2):
        fail("A1 components")
    P(f"  deg<=1: n=18 {low[18]}, n=19 {low[19]}, n=20..30 all [18]: "
      f"{all(low[n] == [18] for n in range(20, 31))}, n=31 {low[31]}, n=32 {low[32]}")
    if low[18] != [16, 17, 18] or low[19] != [16, 18] or low[31] or low[32]:
        fail("A1 leaves")
    P(f"  neighbours of 18 in Q_30: {G[30][18]}; in Q_31: {G[31][18]}")
    P(f"  2x square, x<=40: {[x for x in range(1, 41) if sq(2 * x)]} (2+2=4: pasted Lean Adj is reflexive at value 2)")
    P(f"  pasted ladder 3,7,11,17 differences: {[7 - 3, 11 - 7, 17 - 11]}")
    P(f"  Q_1 components: {ncomp[1]} (so 'N<=13 disconnected' fails at N=1; true for 2<=N<=13: "
      f"{all(ncomp[n] >= 2 for n in range(2, 14))})")
    P()

    # ---- A2 ----
    P("== A2  Hamiltonian path / cycle existence 1..40 (weak-prune DFS, independent)")
    path_yes, cyc_yes = [], []
    for n in range(1, NMAX + 1):
        if n == 1:
            path_yes.append(1)
            continue
        found = None
        for s in sorted(range(1, n + 1), key=lambda v: (len(G[n][v]), v)):
            r = dfs_paths(n, G[n], s, False)
            if r[0]:
                found = r[1]
                break
        if found is not None:
            if not valid_chain(n, found):
                fail(f"A2 witness {n}")
            path_yes.append(n)
        if n >= 3:
            r = dfs_paths(n, G[n], 1, False, need_cycle=True)
            if r[0]:
                if not (valid_chain(n, r[1]) and sq(r[1][0] + r[1][-1])):
                    fail(f"A2 cycle witness {n}")
                cyc_yes.append(n)
    P(f"  path: {path_yes}")
    P(f"  no path (2..40): {[n for n in range(2, NMAX + 1) if n not in path_yes]}")
    P(f"  cycle: {cyc_yes}")
    exp = [1, 15, 16, 17, 23] + list(range(25, 41))
    if path_yes != exp or cyc_yes != list(range(32, 41)):
        fail("A2 sets")
    P("  agrees with lane S2 and with {1} u A090461, A078107, support of A071984: OK")
    P("  session lead paths valid: "
      f"{valid_chain(15, [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9])}, "
      f"{valid_chain(23, [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22])}")
    A071983 = [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 392, 669]   # n=15..33
    A071984 = [1, 1, 11, 57, 31]                                                  # n=32..36
    A090460 = [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 361, 637]
    P("  exact undirected path counts n=15..33 and cycle counts n=32..36:")
    pc = {}
    for n in range(15, 34):
        tot = sum(dfs_paths(n, G[n], s, True)[0] for s in range(1, n + 1))
        if tot % 2:
            fail("odd")
        pc[n] = tot // 2
    cc = {}
    for n in range(32, 37):
        c = dfs_paths(n, G[n], 1, True, need_cycle=True)[0]
        if c % 2:
            fail("odd cycles")
        cc[n] = c // 2
    P(f"  paths : {[pc[n] for n in range(15, 34)]}")
    P(f"  cycles: {[cc[n] for n in range(32, 37)]}")
    if [pc[n] for n in range(15, 34)] != A071983 or [cc[n] for n in range(32, 37)] != A071984:
        fail("A2 counts")
    a90 = [pc[n] - (n - 1) * (cc.get(n, 0)) for n in range(15, 34)]
    if a90 != A090460:
        fail("A2 A090460")
    P("  match A071983 (15..33), A071984 (32..36), A090460 = A071983 - (n-1)*A071984 (15..33): OK")
    P("  n=32: 392 - 31*1 = 361; n=33: 669 - 32*1 = 637")
    P()

    # ---- A3 ----
    P("== A3  cut sets |S|<=3 with c(Q_n - S) > |S| + 1 (brute force)")
    for n in [18, 19, 20, 21, 22, 23, 24]:
        best = None
        for k in range(1, 4):
            for S in combinations(range(1, n + 1), k):
                sub = {v: [w for w in G[n][v] if w not in S] for v in range(1, n + 1) if v not in S}
                c = len(uf_components(n, {v: sub[v] for v in sub}) ) - k  # uf over 1..n counts S as singletons
                if c > k + 1 and (best is None or c > best[1]):
                    best = (S, c)
            if best:
                break
        P(f"  n={n}: {best if best else 'none'}")
        if n in (18, 19) and best[0] != (7,):
            fail("A3")
        if n in (23, 24) and best is not None:
            fail("A3 23/24")
    P()

    # ---- A4 ----
    P("== A4  n=24 certificate details")
    adj = G[24]
    frag_lane = [1, 8, 17, 19, 6, 10, 15, 21, 4, 24, 12]
    frag_path = [4, 21, 15, 10, 6, 19, 17, 8, 1, 24, 12]
    lane_ok = all(sq(frag_lane[i] + frag_lane[i + 1]) for i in range(10))
    path_ok = all(sq(frag_path[i] + frag_path[i + 1]) for i in range(10))
    P(f"  lane's printed fragment order {frag_lane} is a path: {lane_ok} (4+24={4 + 24} not a square)")
    P(f"  correct path order {frag_path}: {path_ok}; ends 4,12; 4+12=16 is an edge: {12 in adj[4]}")
    if lane_ok or not path_ok:
        fail("A4 fragment")
    deg2 = [v for v in adj if len(adj[v]) == 2]
    forced = set()
    for v in adj:
        if len(adj[v]) <= 2:
            for w in adj[v]:
                forced.add((min(v, w), max(v, w)))
    P(f"  degree-2 vertices of Q_24: {deg2} ({len(deg2)} of them) + leaf 18 -> {len(forced)} forced edges in round 1")
    if len(forced) != 21:
        fail("A4 forced count")
    deleted = {(1, 3), (1, 15), (4, 5), (3, 6), (2, 7), (2, 14), (4, 12)}
    red = {v: [w for w in adj[v] if (min(v, w), max(v, w)) not in deleted] for v in adj}
    P(f"  after the seven deletions: N(2)={red[2]}, N(4)={red[4]}, N(18)={red[18]} -> three leaves")
    if not (red[2] == [23] and red[4] == [21] and red[18] == [7]):
        fail("A4 leaves")
    for n in [19, 21, 22]:
        P(f"  n={n}: N(2)={G[n][2]} N(9)={G[n][9]} N(18)={G[n][18]} -> 7 carries 2-7, 7-9, 7-18")
        if not (G[n][2] == [7, 14] and G[n][9] == [7, 16] and G[n][18] == [7]):
            fail("A4 deg3")
    P(f"  n=20: N(4)={G[20][4]} N(11)={G[20][11]} N(20)={G[20][20]} -> 5 carries 4-5, 5-11, 5-20")
    P()

    # ---- A5 ----
    P("== A5  n=15 and n=23")
    P(f"  N_15(4) = {G[15][4]}; degree-2 vertices {[v for v in G[15] if len(G[15][v]) == 2]}; leaves {low[15]}")
    P(f"  chain squares used by the unique n=15 chain: "
      f"{sorted({a + b for a, b in zip([8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9][:-1], [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9][1:])})}")
    ends = {s: dfs_paths(23, G[23], s, True)[0] for s in range(1, 24)}
    ends = {s: c for s, c in ends.items() if c}
    P(f"  n=23 directed path counts by start: {ends}")
    if ends != {2: 1, 9: 1, 18: 3, 22: 1}:
        fail("A5 ends")
    P()

    # ---- A6 ----
    P("== A6  Gerbicz blow-up re-implemented from post #19 pseudo-code")
    a = GERBICZ_V35
    if not valid_chain(35, a):
        fail("A6 v35")
    P(f"  v35 valid chain 1..3: True; 35 = (71*25^0-1)/2 = {(71 - 1) // 2}")
    for m in (1, 2):
        b = gerbicz_fun_odd(a)
        N = len(b)
        blocks = [x for x in b if x >= 13]
        ok = (valid_chain(N, b) and b[0] == 1 and b[-1] == 3 and N == (71 * 25 ** m - 1) // 2
              and sorted(blocks) == list(range(13, N + 1)) and sorted(x for x in b if x <= 12) == list(range(1, 13)))
        P(f"  m={m}: N={N}, chain 1..3 and cycle: {ok}; block values partition [13,{N}]: {sorted(blocks) == list(range(13, N + 1))}")
        if not ok:
            fail("A6")
        a = b
    junctions = len(GERBICZ_F_CALLS) - 1
    P(f"  glue word length {len(GERBICZ_F_CALLS)} items, {junctions} junctions; singletons {sorted(c for c, t in GERBICZ_F_CALLS if t == 0)}; "
      f"c-values {sorted(c for c, t in GERBICZ_F_CALLS if t != 0)}")
    if junctions != 36:
        fail("A6 junctions")
    P()

    # ---- A7 ----
    P("== A7  NEW PROVED: Q_n connected for every n >= 14")
    bad = [n for n in range(5, 200001) if isqrt(2 * n - 1) ** 2 <= n]   # no square in [n+1, 2n-1]
    P(f"  n in [5, 200000] with no square in [n+1, 2n-1]: {bad}")
    if bad:
        fail("A7")
    small = {n: [k * k for k in range(1, 2 * n) if n + 1 <= k * k <= 2 * n - 1] for n in range(5, 10)}
    P(f"  n=5..9 squares in [n+1,2n-1]: {small}")
    P("  n>=10: (sqrt(n+1)+1)^2 = n+2+2sqrt(n+1) <= 2n-1 iff 4(n+1) <= (n-3)^2 iff n^2-10n+5 >= 0, true for n >= 10 "
      f"(check: {[n * n - 10 * n + 5 for n in (9, 10)]}); so vertex n has a neighbour k^2-n in [1, n-1]")
    P("  induction: Q_14 connected (A1) and every n >= 15 attaches to an earlier vertex => Q_n connected for all n >= 14.")
    P("  the pasted 'N>=25 conjectured connected' is therefore PROVED for all N>=14, not conjectural.")
    P()
    P("== A8  labels")
    P("  S2a (three leaves) and S2b (cut criterion c(G-S) <= |S|+1) are classical textbook facts; the attribution of the")
    P("  cycle form of S2b to Chvatal 1973 (toughness) is UNCITED-RECOLLECTION; neither is new to the repository.")
    P()
    P("== done")


if __name__ == "__main__":
    main()
