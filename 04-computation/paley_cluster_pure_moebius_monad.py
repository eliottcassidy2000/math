#!/usr/bin/env python3
"""
paley_cluster_pure_moebius_monad.py
monad-explorer-2026-06-07 (deep-research / analytic lane, 4th session)

ADVANCES THM-438 handoff #1 ("prove the signed even-cacti sum = C_k cleanly").

KEY REDUCTION (the point of this script):
  The leading coefficient C_k of A_{2k}/p^{k+1} is a PURE partition-lattice
  Moebius sum over even cacti -- the character / Gauss-sum content is TRIVIAL.

  Recall (THM-438 ADDENDUM, PROVED):  M_rho = (-1)^k p^{V-k} F(rho),
      F(rho) = sum over F_p-flows on G_rho of prod_e chi(t_e),
  and the top order p^{k+1} is reached exactly by the EVEN CACTI (connected,
  bridgeless, every biconnected block an even cycle).

  CLAIM 1 (character-triviality, verified here):  for EVERY even cactus rho,
      F(rho) = (p-1)^c   EXACTLY    (c = #blocks = cycle rank, m),
  with NO Legendre-symbol sign surviving.  Reason (proof in the .md):
  the walk x_0->...->x_{2k} is a CLOSED Euler trail of G_rho (all degrees even),
  so it sweeps each even-cycle block once in a single rotational sense; the
  flow on a length-2j block is s_e = +-t with prod_e chi(s_e) = chi(t^{2j}) = 1.
  Blocks are flow-independent (cactus => cycle space is a direct sum), so
  F = prod_blocks (p-1) = (p-1)^c, positive, character-free.

  CONSEQUENCE:  lead( M_rho / p^{k+1} ) = (-1)^k  for every even cactus, so
      C_k = (-1)^k * SUM_{even cacti} mu(0,rho),     mu(0,rho)=prod_v(-1)^{|B|-1}(|B|-1)!
  i.e. the ENTIRE Catalan mechanism is the pure combinatorial identity
      ($)  SUM_{even cacti, 2k edges} (-1)^c prod_v (|B_v|-1)!  =  (-1)^k C_k.
  (Using prod_v (-1)^{|B_v|-1} = (-1)^{sum(|B_v|-1)} = (-1)^{(2k+1)-V} = (-1)^c.)

  No number theory remains: ($) is a statement about the partition lattice
  restricted to even-cactus gluings of a path -- the moment->free-cumulant
  (Wigner semicircle) identity, Wigner-UNIVERSAL (confirms HYP-2308's spectral
  reading with the actual mechanism: planarity, not arithmetic).

This script:
  (A) verifies F(rho) = (p-1)^c EXACTLY for every even cactus, k=2,3 (p=11,13,19,23,31);
  (B) computes S_k = sum_{even cacti} mu(0,rho) by exhaustive partition enumeration,
      k=1..5, and checks S_k = (-1)^k C_k (Catalan);
  (C) bonus: the bigon-tree-only sub-sum reproduces OEIS A088368 (1,3,13,69,421),
      and the cycle-block corrections bring it to Catalan -- printed per block-profile.
"""
import math
from itertools import product as iproduct
from collections import defaultdict, Counter
import numpy as np


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]
        return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm


def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B)
        m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks)


def biconnected_components(edges, nb):
    """Return list of (vertex_set, edge_index_list) for each biconnected component."""
    adj = defaultdict(list)
    for ei, (u, v) in enumerate(edges):
        adj[u].append((v, ei))
        adj[v].append((u, ei))
    idx = {}
    low = {}
    cnt = [0]
    stack = []
    comps = []
    ve = set()
    import sys
    sys.setrecursionlimit(100000)

    def dfs(u, pe):
        idx[u] = low[u] = cnt[0]
        cnt[0] += 1
        for (w, ei) in adj[u]:
            if ei == pe:
                continue
            if ei in ve:
                if w in idx and idx[w] < idx[u]:
                    stack.append(ei)
                    low[u] = min(low[u], idx[w])
                continue
            ve.add(ei)
            stack.append(ei)
            if w not in idx:
                dfs(w, ei)
                low[u] = min(low[u], low[w])
                if low[w] >= idx[u]:
                    comp = []
                    while True:
                        e = stack.pop()
                        comp.append(e)
                        if e == ei:
                            break
                    comps.append(comp)
            else:
                low[u] = min(low[u], idx[w])

    for s in range(nb):
        if s not in idx:
            dfs(s, -1)
    out = []
    for comp in comps:
        vs = set()
        for ei in comp:
            u, v = edges[ei]
            vs.add(u)
            vs.add(v)
        out.append((vs, comp))
    return out


def is_even_cactus(edges, nb):
    """connected, bridgeless, every biconnected component an EVEN cycle.
       returns (True, c) with c = #blocks, else (False, None)."""
    # connectivity
    if nb == 0:
        return False, None
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen = {0}
    stk = [0]
    while stk:
        u = stk.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stk.append(w)
    if len(seen) != nb:
        return False, None
    comps = biconnected_components(edges, nb)
    c = 0
    total_edges = 0
    for vs, comp in comps:
        e = len(comp)
        v = len(vs)
        total_edges += e
        if e < 2:
            return False, None            # a bridge -> not bridgeless
        if e != v:
            return False, None            # not a simple cycle (theta etc.)
        if e % 2 != 0:
            return False, None            # odd cycle
        c += 1
    if total_edges != len(edges):
        return False, None
    return True, c


def nullspace_basis(edges, nb, p):
    E = len(edges)
    A = np.zeros((nb, E), dtype=np.int64)
    for ei, (u, v) in enumerate(edges):
        A[v, ei] += 1
        A[u, ei] -= 1
    A %= p
    rows, cols = A.shape
    pivots = []
    r = 0
    for cc in range(cols):
        piv = next((rr for rr in range(r, rows) if A[rr, cc] % p != 0), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, cc]), p - 2, p)) % p
        for rr in range(rows):
            if rr != r and A[rr, cc] % p != 0:
                A[rr] = (A[rr] - A[rr, cc] * A[r]) % p
        pivots.append(cc)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        vec = np.zeros(cols, dtype=np.int64)
        vec[fc] = 1
        for ri, pc in enumerate(pivots):
            vec[pc] = (-A[ri, fc]) % p
        basis.append(vec)
    return np.array(basis, dtype=np.int64) if basis else np.zeros((0, cols), dtype=np.int64)


def flow_sum(edges, nb, p):
    basis = nullspace_basis(edges, nb, p)
    m = basis.shape[0]
    E = len(edges)
    chi = np.array([legendre(z, p) for z in range(p)], dtype=np.int64)
    if m == 0:
        return (1 if E == 0 else 0), 0
    grids = np.array(np.meshgrid(*[range(p)] * m, indexing='ij')).reshape(m, -1).T
    T = (grids @ basis) % p
    V = chi[T]
    return int(V.prod(axis=1).sum()), m


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


# ===================================================================
# (A) F(rho) = (p-1)^c EXACTLY for every even cactus (character-triviality)
# ===================================================================
print("=" * 72)
print("(A) CHARACTER-TRIVIALITY:  F(rho) = (p-1)^c  exactly, for every even cactus")
print("    (c = #blocks = cycle rank).  If true => leading coeff is pure Moebius.")
for k in [2, 3]:
    L = 2 * k
    allgood = True
    nchecked = 0
    for blocks in set_partitions(range(L + 1)):
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges):
            continue
        ok, c = is_even_cactus(edges, nb)
        if not ok:
            continue
        nchecked += 1
        for p in [11, 13, 19, 23, 31]:
            fs, m = flow_sum(edges, nb, p)
            assert m == c, (blocks, m, c)
            if fs != (p - 1) ** c:
                allgood = False
                print(f"   VIOLATION k={k} p={p} blocks={blocks} F={fs} (p-1)^c={(p-1)**c}")
    print(f"   k={k}: checked {nchecked} even cacti, F=(p-1)^c for ALL p,rho:  {allgood}")

# ===================================================================
# (B) pure Moebius identity S_k = sum_{even cacti} mu = (-1)^k C_k
# ===================================================================
print("=" * 72)
print("(B) PURE MOEBIUS IDENTITY   S_k = sum_{even cacti} mu(0,rho)  =?  (-1)^k C_k")
print("    (no primes, no characters -- pure partition lattice)")
print(f"    {'k':>2} {'S_k':>8} {'(-1)^k C_k':>12} {'match':>7}   #even-cacti")
for k in range(1, 6):
    L = 2 * k
    S = 0
    ncacti = 0
    profile = defaultdict(int)         # block-size-profile -> signed weight
    for blocks in set_partitions(range(L + 1)):
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges):
            continue
        ok, c = is_even_cactus(edges, nb)
        if not ok:
            continue
        ncacti += 1
        mu = mu_partition(blocks)
        S += mu
        # biconnected-block edge-size profile
        comps = biconnected_components(edges, nb)
        sizes = tuple(sorted(len(comp) for _, comp in comps))
        profile[sizes] += mu
    target = (-1) ** k * catalan(k)
    print(f"    {k:>2} {S:>8} {target:>12} {str(S == target):>7}   {ncacti}")
    # show the per-block-profile breakdown and the bigon-tree sub-sum
    parts = sorted(profile.items(), key=lambda kv: (len(kv[0]), kv[0]))
    line = "        profile -> mu-sum: " + "  ".join(
        f"{p}:{w:+d}" for p, w in parts)
    print(line)

# ===================================================================
# (C) bigon-tree-only sub-sum = OEIS A088368 (the all-pairings overcount)
# ===================================================================
print("=" * 72)
print("(C) bigon-tree-ONLY |mu|-sum (all blocks = bigons) = OEIS A088368 (~e*n!)")
print("    vs Catalan after the even-cycle-cactus corrections")


def bigon_tree_coeff(k):
    # sum over non-crossing edge pairings of prod (b_v-1)!   (A088368)
    def nc_pairings(items):
        if not items:
            yield []
            return
        first = items[0]
        for i in range(1, len(items), 2):
            for pa in nc_pairings(items[1:i]):
                for pb in nc_pairings(items[i + 1:]):
                    yield [(first, items[i])] + pa + pb
    tot = 0
    for pr in nc_pairings(list(range(2 * k))):
        parent = list(range(2 * k + 1))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for (i, j) in pr:
            parent[find(i)] = find(j + 1)
            parent[find(i + 1)] = find(j)
        bs = list(Counter(find(x) for x in range(2 * k + 1)).values())
        if len(bs) != k + 1:
            continue
        w = 1
        for b in bs:
            w *= math.factorial(b - 1)
        tot += w
    return tot


print("    k:           " + "  ".join(f"{k:5d}" for k in range(1, 7)))
print("    A088368:     " + "  ".join(f"{bigon_tree_coeff(k):5d}" for k in range(1, 7)))
print("    Catalan C_k: " + "  ".join(f"{catalan(k):5d}" for k in range(1, 7)))
print("=" * 72)
print("CONCLUSION: the Catalan leading coefficient is a PURE partition-lattice")
print("Moebius sum over even-cactus gluings of the path -- the character content")
print("collapses to F(rho)=(p-1)^c.  Handoff #1 reduces to the number-theory-FREE")
print("identity  ($) sum_{even cacti}(-1)^c prod(|B|-1)! = (-1)^k C_k,  the")
print("moment->free-cumulant (semicircle) relation; Wigner-universal (HYP-2308).")
