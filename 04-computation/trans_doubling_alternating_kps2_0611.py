#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HYP-2360 stress test: the alternating configuration A_l.
kind-pasteur-2026-06-11-S2 (Erdős–Moser #1216 thread, THM-455/HYP-2360).

D(T) arc rules (THM-447 skew doubling, derived from H_2n = [[H,H],[-H^T,H^T]],
H = I + S; verified against THM-455's lower-bound proof and the n=5 exception):
  (u,0)->(v,0) iff u->v          (copy 0 = T)
  (u,1)->(v,1) iff v->u          (copy 1 = T^op)
  (u,0)->(v,1) iff u->v or u=v   (cross follows T; twin 0->1)
  (u,1)->(v,0) iff u->v, u != v  (cross follows T)

A transitive k-chain in D(T) with plains U (eps=0) and primes W (eps=1) forces in T:
  plains ascend, primes (read backwards) ascend, cross: plain-before-prime u->w (or =),
  prime-before-plain w->u.  The FULLY ALTERNATING pattern w1' u1 w2' u2 ... w_{l+1}'
  forces exactly the tournament A_l on {u1..ul, w1..w_{l+1}}:
    u_i -> u_j (i<j);  w_j -> w_i (i<j... i.e. later prime beats earlier);
    w_i -> u_j iff i <= j;  u_i -> w_j iff j >= i+1.
  By construction D(A_l) contains a transitive (2l+1)-chain.
  QUESTION: trans(A_l) = ?  If l+1, then delta(A_l) >= l, refuting HYP-2360 (delta<=2)
  for l >= 3 — beyond the n=6 census horizon (A_3 lives on n=7).

Exact: trans via bitset DP over subsets (n <= 25-ish fine for n <= 15 here);
trans(D(A_l)) exact for small l, and the explicit chain VERIFIED as transitive.
"""

import itertools, time


def trans_exact(n, adj):
    """Maximum transitive subtournament via DP: longest chain in the DAG of
    'extension' — classic: for each subset... use simpler O(2^n * n): f[S] =
    max over v in S that BEATS ALL of S-v ... that's for transitive with source...
    Standard trick: a transitive subtournament = a chain v1->...->vk with all
    forward arcs. DP over subsets is exponential; instead use: longest path in
    the 'domination order' is wrong. Use recursion with memo on (last, allowed)?
    For n <= 15 just do DP: g[v] over orderings... Simplest correct: max chain
    where chain condition = all pairs; equivalent to longest sequence v1..vk
    with v_i -> v_j for all i<j. This equals the longest path in the partial
    structure: greedy doesn't work; do DFS with pruning over descending sets:
    chain extension: next vertex must be beaten by ALL current members — track
    the set of common out-neighbors as a bitmask. Complete search with memo on
    the candidate mask (the chain's future is determined by the mask alone)."""
    out = [0] * n
    for u in range(n):
        for v in range(n):
            if u != v and adj[u][v]:
                out[u] |= 1 << v
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def best(mask):
        # longest chain we can build choosing vertices from mask, where each
        # chosen vertex must beat all later ones: pick the source v in mask,
        # recurse on mask & out[v]
        b = 0
        m = mask
        while m:
            v = (m & -m).bit_length() - 1
            m &= m - 1
            b = max(b, 1 + best(mask & out[v]))
        return b
    return best((1 << n) - 1)


def build_A(l):
    """A_l on n = 2l+1 vertices: u_1..u_l = 0..l-1, w_1..w_{l+1} = l..2l."""
    n = 2 * l + 1
    adj = [[False] * n for _ in range(n)]
    U = list(range(l))
    W = list(range(l, 2 * l + 1))
    for i in range(l):
        for j in range(i + 1, l):
            adj[U[i]][U[j]] = True
    for i in range(l + 1):
        for j in range(i + 1, l + 1):
            adj[W[j]][W[i]] = True  # later prime beats earlier
    for i in range(l + 1):          # w_{i+1} (1-indexed i+1)
        for j in range(l):          # u_{j+1}
            # w_a -> u_b iff a <= b (1-indexed); here a = i+1, b = j+1
            if i <= j:
                adj[W[i]][U[j]] = True
            else:
                adj[U[j]][W[i]] = True
    return n, adj, U, W


def double(n, adj):
    """D(T) on 2n vertices: (v,0) = v, (v,1) = n+v."""
    N = 2 * n
    dadj = [[False] * N for _ in range(N)]
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            if adj[u][v]:
                dadj[u][v] = True              # copy0 = T
                dadj[n + v][n + u] = True      # copy1 = T^op
                dadj[u][n + v] = True          # cross 0->1 follows T
                dadj[n + u][v] = True          # cross 1->0 follows T
            # else handled by the (v,u) iteration
    for v in range(n):
        dadj[v][n + v] = True                  # twin
    # sanity: tournament?
    for a in range(N):
        for b in range(a + 1, N):
            assert dadj[a][b] != dadj[b][a], (a, b)
    return N, dadj


def verify_chain(dadj, chain):
    for i in range(len(chain)):
        for j in range(i + 1, len(chain)):
            if not dadj[chain[i]][chain[j]]:
                return False
    return True


def trans_brute(n, adj, k):
    """Independent check: does a transitive subtournament of size k exist?
    Pure combinations + pairwise check (different method from the DP)."""
    import itertools as it
    for sub in it.combinations(range(n), k):
        # try to topologically order sub by wins within sub
        wins = {v: sum(1 for w in sub if w != v and adj[v][w]) for v in sub}
        order = sorted(sub, key=lambda v: -wins[v])
        if all(adj[order[i]][order[j]] for i in range(k) for j in range(i + 1, k)):
            return True
    return False


def main():
    t0 = time.time()
    print("=== A_l: alternating configuration vs HYP-2360 (delta <= 2?) ===", flush=True)
    for l in range(2, 7):
        n, adj, U, W = build_A(l)
        t = trans_exact(n, adj)
        N, dadj = double(n, adj)
        # the designed chain: w1', u1, w2', u2, ..., w_{l+1}'
        chain = []
        for i in range(l):
            chain.append(n + W[i])  # w_{i+1} primed
            chain.append(U[i])      # u_{i+1} plain
        chain.append(n + W[l])
        ok = verify_chain(dadj, chain)
        tD = trans_exact(N, dadj) if N <= 22 else None
        delta = (tD - t) if tD is not None else None
        # independent brute-force cross-checks (different method, MISTAKE-067)
        bf_t = trans_brute(n, adj, t) and not trans_brute(n, adj, t + 1)
        print(f"   l={l} (n={n}): trans(A_l)={t} [brute confirms: {bf_t}], "
              f"designed {2*l+1}-chain in D verifies={ok}, trans(D)={tD}, delta={delta}"
              + ("   *** REFUTES HYP-2360 ***" if delta is not None and delta > 2
                 else (f"   *** chain alone gives delta >= {2*l+1-t} ***" if ok and 2*l+1-t > 2 else "")),
              flush=True)
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
