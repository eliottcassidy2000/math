#!/usr/bin/env python3
"""
paley_starstar_triangle_fast_monad.py
monad-explorer-2026-06-07 (deep-research, 8th session)

FAST integer enumerator for the cycle-rank (genus-like) triangle of (**):

   S_k = sum_{sigma: EVEN-SERIES pattern of path [0..2k]} mu(0,sigma) = (-1)^k C_k
   t(k,m) = sum_{even-series sigma, cycle-rank m} prod_v (|B_v|-1)!      (POSITIVE triangle)
   S_k = sum_m (-1)^m t(k,m)                                            (sign identity, ADD-5)

Replaces the per-partition numpy SVD (THM-438 ADDENDUM-2/3 scripts) with:
  - iterative RGS (restricted-growth-string) set-partition generation,
  - cheap prefilters (no self-loop, no-leaf/min-degree>=2, connected via union-find),
  - INTEGER fundamental-cycle line detection for the even-series test
    (each edge -> its coefficient vector over the m fundamental cycles; group by
     sign-normalized primitive direction; even-series <=> every group even).

This is pure-Python integer (deterministic, no float tol), validated against the
known triangle rows (k<=5), and pushes to k=6 (Bell(13)=27.6M) / k=7 if time allows.

Outputs, per k:
  - S_k and match to (-1)^k C_k
  - row t(k,1..k)
  - refined: number of even-series patterns by #flow-lines e (for the e-grading),
    and signed-weighted Psi(k,e) = sum (-1)^m W over patterns with e lines.

Usage: python3 paley_starstar_triangle_fast_monad.py [KMAX]
"""
import sys
from math import comb, factorial
from collections import defaultdict


def catalan(k):
    return comb(2 * k, k) // (k + 1)


# ----- iterative RGS generation of set partitions of {0,...,n-1} -----
def rgs_iter(n):
    """Yield restricted growth strings a[0..n-1] (a[0]=0, a[i]<=max(a[:i])+1)."""
    a = [0] * n
    b = [0] * n          # b[i] = max(a[0..i-1])+1 cap helper: b[i]=a[i] running max+1
    # standard Knuth Algorithm H (set partitions in RGS form)
    if n == 1:
        yield a
        return
    m = [0] * n          # m[i] = 1 + max(a[0..i])
    yield a
    while True:
        # find rightmost j we can increment
        j = n - 1
        while j > 0 and a[j] == m[j - 1] + 1:
            j -= 1
        if j == 0:
            return
        a[j] += 1
        # recompute m[j]
        mj = m[j - 1] if a[j] <= m[j - 1] else a[j]
        m[j] = mj
        for i in range(j + 1, n):
            a[i] = 0
            m[i] = mj
        yield a


def line_signature(edges, V):
    """Return list of sign-normalized integer primitive line-vectors for each edge,
       computed via fundamental cycles of a spanning tree.  edges: list of (u,v)
       DIRECTED (u=tail head=v). Returns (sigs, m) where m=cycle rank E-V+1.
       Tree edges & non-tree edges all get an m-vector; proportional <=> same line."""
    E = len(edges)
    # build adjacency with edge indices for spanning tree (BFS), undirected
    adj = [[] for _ in range(V)]
    for ei, (u, v) in enumerate(edges):
        adj[u].append((v, ei))
        adj[v].append((u, ei))
    parent = [-1] * V
    parent_edge = [-1] * V
    visited = [False] * V
    # BFS from 0 (graph assumed connected by caller)
    order = []
    stack = [0]
    visited[0] = True
    tree_edge = [False] * E
    while stack:
        x = stack.pop()
        order.append(x)
        for (y, ei) in adj[x]:
            if not visited[y]:
                visited[y] = True
                parent[y] = x
                parent_edge[y] = ei
                tree_edge[ei] = True
                stack.append(y)
    # depth for tree-path computations
    depth = [0] * V
    for x in order:
        if parent[x] != -1:
            depth[x] = depth[parent[x]] + 1
    non_tree = [ei for ei in range(E) if not tree_edge[ei]]
    m = len(non_tree)            # = E - (V-1) = E-V+1
    cyc_index = {ei: j for j, ei in enumerate(non_tree)}
    # coefficient vector per edge over the m fundamental cycles
    vecs = [[0] * m for _ in range(E)]
    # each non-tree edge f=(u,v): fundamental cycle = f + tree-path(v->u).
    # orient cycle in direction of edge f (tail->head along f).
    for f in non_tree:
        j = cyc_index[f]
        u, v = edges[f]
        vecs[f][j] += 1          # the chord itself, +1 in direction tail->head
        # tree path from v up to LCA and u up to LCA; we add edges with orientation
        a, b = u, v
        # We want cycle: tail u -> head v via chord, then v -> u via tree.
        # Add tree edges on path v..u oriented so that the cycle is consistent
        # (head v back to tail u). Walk both up to common ancestor.
        pa, pb = a, b
        path_a = []
        path_b = []
        while depth[pa] > depth[pb]:
            path_a.append(pa); pa = parent[pa]
        while depth[pb] > depth[pa]:
            path_b.append(pb); pb = parent[pb]
        while pa != pb:
            path_a.append(pa); pa = parent[pa]
            path_b.append(pb); pb = parent[pb]
        lca = pa
        # tree edges from u up to lca: edge (parent_edge[node]) connecting node-parent
        # cycle direction: u -> ... -> lca -> ... -> v  is the "tree part" of chord-cycle
        # but our chord is u->v, so tree part is v->u. We'll just assign consistent signs;
        # only PROPORTIONALITY (line) matters, exact sign is normalized away.
        for node in path_a:           # u side, going up: orient node->parent  => +1
            te = parent_edge[node]
            tu, tv = edges[te]
            # te connects node and parent[node]; figure its + direction (tail->head)
            s = +1
            # we add coefficient for cycle j; sign by whether te direction matches up-walk
            if (tu, tv) == (node, parent[node]):
                s = +1
            else:
                s = -1
            vecs[te][j] += s
        for node in path_b:           # v side: orient parent->node (so overall consistent)
            te = parent_edge[node]
            tu, tv = edges[te]
            s = -1
            if (tu, tv) == (node, parent[node]):
                s = -1
            else:
                s = +1
            vecs[te][j] += s
    # sign-normalize each vec to a primitive canonical form (gcd, leading sign +)
    from math import gcd
    sigs = []
    for vec in vecs:
        g = 0
        for c in vec:
            g = gcd(g, abs(c))
        if g == 0:
            sigs.append(None)        # bridge: zero vector -> not even-series-saturating
            continue
        red = tuple(c // g for c in vec)
        # leading-sign normalize
        for c in red:
            if c != 0:
                if c < 0:
                    red = tuple(-x for x in red)
                break
        sigs.append(red)
    return sigs, m


def analyze(a, L):
    """a = RGS of {0..L}; return (m, weight) if even-series pattern else None."""
    n = L + 1
    V = max(a) + 1
    # build directed edges i->i+1 on block labels
    edges = [(a[i], a[i + 1]) for i in range(L)]
    # no self-loop
    for (u, v) in edges:
        if u == v:
            return None
    # degrees (no-leaf): each block incident to >=2 edge-endpoints
    deg = [0] * V
    for (u, v) in edges:
        deg[u] += 1
        deg[v] += 1
    if any(d < 2 for d in deg):
        return None
    # connected via union-find
    par = list(range(V))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (u, v) in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            par[ru] = rv
    root = find(0)
    if any(find(x) != root for x in range(V)):
        return None
    sigs, m = line_signature(edges, V)
    if m == 0:
        return None
    if any(s is None for s in sigs):
        return None                  # bridge present
    grp = defaultdict(int)
    for s in sigs:
        grp[s] += 1
    if any(c % 2 for c in grp.values()):
        return None
    e = len(grp)                     # number of flow-lines (topological edges)
    # weight = prod_v (|B_v|-1)! ; |B_v| = #positions in block v
    size = [0] * V
    for x in a:
        size[x] += 1
    W = 1
    for s in size:
        W *= factorial(s - 1)
    return m, W, e


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    print("=" * 72)
    print("FAST cycle-rank triangle for (**)  S_k = sum_{even-series} mu = (-1)^k C_k")
    print("=" * 72)
    A215257 = {1: 1, 2: 3, 3: 13, 4: 67, 5: 383, 6: 2345, 7: 16319}
    for k in range(1, KMAX + 1):
        L = 2 * k
        n = L + 1
        S = 0
        cnt = 0
        t_row = defaultdict(int)      # t(k,m) positive
        psi_e = defaultdict(int)      # signed Psi(k,e)
        cnt_e = defaultdict(int)      # #patterns by e (unsigned)
        for a in rgs_iter(n):
            res = analyze(a, L)
            if res is None:
                continue
            m, W, e = res
            sign = -1 if (m & 1) else 1
            S += sign * W
            cnt += 1
            t_row[m] += W
            psi_e[e] += sign * W
            cnt_e[e] += 1
        tgt = (-1) ** k * catalan(k)
        match = (S == tgt)
        cmatch = (cnt == A215257.get(k))
        row = " ".join(str(t_row[m]) for m in range(1, k + 1))
        print(f"k={k}: S_k={S:+d}  target={tgt:+d}  match={match}   "
              f"#even-series={cnt} (A215257(k+1)? {cmatch})")
        print(f"      t(k,m) m=1..{k}: {row}")
        print(f"      check sum_m (-1)^m t(k,m) = "
              f"{sum((-1)**m * t_row[m] for m in range(1, k+1)):+d}")
        psi_str = ", ".join(f"e={e}:{psi_e[e]:+d}" for e in sorted(psi_e))
        print(f"      Psi(k,e) signed-by-#lines: {psi_str}")
        cnt_str = ", ".join(f"e={e}:{cnt_e[e]}" for e in sorted(cnt_e))
        print(f"      #patterns by #lines e:     {cnt_str}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
