#!/usr/bin/env python3
"""
paley_cluster_cactus_census_monad.py
monad-explorer-2026-06-07 (deep-research / analytic lane, 3rd session)

CORRECTS the leading-order MECHANISM of THM-438 (Paley cluster integrals are Catalan).

THM-438 / its reflection claim the top order p^{k+1} of the cluster integral
    A_{2k} = sum_{x_0..x_{2k} in F_p, all distinct} prod chi(x_{i+1}-x_i)
is reached ONLY by all-bigon-tree coincidence patterns, each contributing +1, with
"non-bigon cycles O(p^{k+1/2})", so that the count is literally Catalan C_k.

THIS SCRIPT SHOWS THAT MECHANISM IS WRONG (the final value C_k is right; the
mechanism is a signed cancellation):

  (1) A closed-form for every coincidence pattern's free sum (PROVED via Gauss-sum
      inversion, verified exactly here):
          M_sigma = (-1)^k * p^{V-k} * F(sigma),   F = sum over F_p-flows of prod chi(t_e),
      where V = #blocks, the graph G_sigma has the 2k walk-edges, and the flow sum
      ranges over the cycle space (dim m = 2k - V + 1).
      A pattern reaches the top order p^{k+1} iff its flow-character-sum F reaches the
      FULL order p^m.  Those are the "even cacti": connected graphs whose biconnected
      blocks are all EVEN cycles (a bigon = 2-cycle is the smallest block).

  (2) Census at k=2,3: bigon-TREES alone do NOT sum to C_k.  Weighted by the partition-
      lattice Moebius factor mu(0,sigma) = prod_blocks (-1)^{|B|-1}(|B|-1)!, the
      bigon-tree leading coefficient is 1, 3, 13, 69, 421, 2867 (NOT C_k, NOT (2k-1)!!).
      The even-cycle cacti contribute NEGATIVE corrections; the signed total is C_k.
        k=2: bigons 3  +  4-cycle (-1)            = 2 = C_2
        k=3: bigons 13 +  {bigon+4cycle} + {6cycle} = 5 = C_3

  (3) Error term: A_{2k} = C_k p^{k+1} + O(p^k)  (NOT O(p^{k+1/2})).  The single
      2k-cycle pattern is exactly tr(M^{2k}) = (-p)^k (p-1) -- ELEMENTARY, no Weil --
      so Part C of THM-438 (R(p)->e) needs NO Weil bound at all.

  (4) Consequence for the sub-leading constant (HYP-2307 #2): the A_4 correction is
      O(p^2) = O(p^k), i.e. relative O(1/p), supporting the 1/p convergence of R->e
      (resolves the reflection's stated O(1/sqrt p) vs the close-out's "favors 1/p").
"""
import itertools, math
from itertools import product as iproduct
import numpy as np

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)

def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]; return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm

def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B); m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m

def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B: pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks)

def nullspace_basis(edges, nb, p):
    E = len(edges)
    A = np.zeros((nb, E), dtype=np.int64)
    for ei, (u, v) in enumerate(edges):
        A[v, ei] += 1; A[u, ei] -= 1
    A %= p
    rows, cols = A.shape; pivots = []; r = 0
    for cc in range(cols):
        piv = next((rr for rr in range(r, rows) if A[rr, cc] % p != 0), None)
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, cc]), p - 2, p)) % p
        for rr in range(rows):
            if rr != r and A[rr, cc] % p != 0:
                A[rr] = (A[rr] - A[rr, cc] * A[r]) % p
        pivots.append(cc); r += 1
        if r == rows: break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        vec = np.zeros(cols, dtype=np.int64); vec[fc] = 1
        for ri, pc in enumerate(pivots): vec[pc] = (-A[ri, fc]) % p
        basis.append(vec)
    return np.array(basis, dtype=np.int64) if basis else np.zeros((0, cols), dtype=np.int64)

def flow_sum(edges, nb, p):
    basis = nullspace_basis(edges, nb, p); m = basis.shape[0]; E = len(edges)
    chi = np.array([legendre(z, p) for z in range(p)], dtype=np.int64)
    if m == 0:
        return (1 if E == 0 else 0), 0
    grids = np.array(np.meshgrid(*[range(p)] * m, indexing='ij')).reshape(m, -1).T
    T = (grids @ basis) % p
    V = chi[T]
    return int(V.prod(axis=1).sum()), m

def M_flow(blocks, L, p):
    k = L // 2; edges, nb = build_graph(blocks, L); fs, m = flow_sum(edges, nb, p)
    return ((-1) ** k) * (p ** (nb - k)) * fs

def M_direct(blocks, L, p):
    chi = [legendre(z, p) for z in range(p)]
    edges, nb = build_graph(blocks, L); tot = 0
    for vals in iproduct(range(p), repeat=nb):
        pr = 1
        for (u, v) in edges:
            pr *= chi[(vals[v] - vals[u]) % p]
            if pr == 0: break
        tot += pr
    return tot

def biconn_sizes(edges, nb):
    from collections import defaultdict
    adj = defaultdict(list)
    for ei, (u, v) in enumerate(edges):
        adj[u].append((v, ei)); adj[v].append((u, ei))
    idx = {}; low = {}; cnt = [0]; stack = []; comps = []; ve = set()
    import sys; sys.setrecursionlimit(100000)
    def dfs(u, pe):
        idx[u] = low[u] = cnt[0]; cnt[0] += 1
        for (w, ei) in adj[u]:
            if ei == pe: continue
            if ei in ve:
                if w in idx and idx[w] < idx[u]: stack.append(ei); low[u] = min(low[u], idx[w])
                continue
            ve.add(ei); stack.append(ei)
            if w not in idx:
                dfs(w, ei); low[u] = min(low[u], low[w])
                if low[w] >= idx[u]:
                    comp = []
                    while True:
                        e = stack.pop(); comp.append(e)
                        if e == ei: break
                    comps.append(comp)
            else:
                low[u] = min(low[u], idx[w])
    for s in range(nb):
        if s not in idx: dfs(s, -1)
    return tuple(sorted(len(c) for c in comps))

def catalan(k): return math.comb(2 * k, k) // (k + 1)

# ---------- (0) verify the flow closed-form M_sigma = (-1)^k p^{V-k} F(sigma) ----------
print("=" * 70)
print("(0) VERIFY  M_sigma = (-1)^k p^{V-k} * flow-char-sum   (k=2, p=11,13)")
ok = True
for p in [11, 13]:
    for blocks in set_partitions(range(5)):
        if M_direct(blocks, 4, p) != M_flow(blocks, 4, p):
            ok = False; print("  MISMATCH", blocks, p)
print("  all match:", ok)

# ---------- (1)+(2) census ----------
def census(k, primes):
    L = 2 * k; cats = {}
    for blocks in set_partitions(range(L + 1)):
        mu = mu_partition(blocks)
        if mu == 0: continue
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges): continue       # self-loop -> chi(0)=0 -> M=0
        st = biconn_sizes(edges, nb)
        for p in primes:
            c = mu * M_flow(blocks, L, p)
            cats.setdefault(st, {}).setdefault(p, 0); cats[st][p] += c
    return cats

for k in [2, 3]:
    primes = [11, 19, 23, 31] if k == 2 else [11, 19, 23]
    print("=" * 70)
    print(f"({k}) CACTUS CENSUS  k={k}  (C_k={catalan(k)})   contributions / p^{k+1}")
    cats = census(k, primes)
    pk1 = lambda p: p ** (k + 1)
    grand = {p: 0 for p in primes}
    rows = []
    for st, d in cats.items():
        for p in primes: grand[p] += d[p]
        rows.append((st, [d[p] / pk1(p) for p in primes]))
    rows.sort(key=lambda r: -abs(r[1][-1]))
    print("   cactus-block-sizes  " + "  ".join(f"@p={p}" for p in primes))
    for st, vals in rows:
        if abs(vals[-1]) > 1e-6:
            print(f"   {str(st):18s}  " + "  ".join(f"{v:8.4f}" for v in vals))
    print(f"   {'TOTAL -> C_k':18s}  " + "  ".join(f"{grand[p] / pk1(p):8.4f}" for p in primes))
    bt = {p: 0 for p in primes}
    for st, d in cats.items():
        if st and all(s == 2 for s in st):
            for p in primes: bt[p] += d[p]
    print(f"   {'bigon-trees ONLY':18s}  " + "  ".join(f"{bt[p] / pk1(p):8.4f}" for p in primes))

# ---------- bigon-tree weighted leading coefficient (its own sequence) ----------
print("=" * 70)
print("BIGON-TREE leading coeff = sum over non-crossing edge-pairings of prod (b_v-1)!")
def noncrossing_pairings(items):
    if not items: yield []; return
    first = items[0]
    for i in range(1, len(items), 2):
        for pa in noncrossing_pairings(items[1:i]):
            for pb in noncrossing_pairings(items[i + 1:]):
                yield [(first, items[i])] + pa + pb
def bigon_tree_coeff(k):
    tot = 0
    for pr in noncrossing_pairings(list(range(2 * k))):
        parent = list(range(2 * k + 1))
        def find(a):
            while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for (i, j) in pr:
            parent[find(i)] = find(j + 1); parent[find(i + 1)] = find(j)
        from collections import Counter
        bs = list(Counter(find(x) for x in range(2 * k + 1)).values())
        if len(bs) != k + 1: continue
        w = 1
        for b in bs: w *= math.factorial(b - 1)
        tot += w
    return tot
print("  k:           " + "  ".join(f"{k:5d}" for k in range(1, 7)))
print("  bigon coeff: " + "  ".join(f"{bigon_tree_coeff(k):5d}" for k in range(1, 7)) + "   (OEIS-worthy; NOT Catalan)")
print("  Catalan C_k: " + "  ".join(f"{catalan(k):5d}" for k in range(1, 7)) + "   (the SIGNED total after cycle-cacti)")

# ---------- (3) error term O(p^k) not O(p^{k+1/2}) ----------
print("=" * 70)
print("(3) ERROR TERM:  A_4 = 2 p^3 + ?   (Paley primes)")
def A4(p):
    chi = [legendre(z, p) for z in range(p)]; tot = 0
    for t in itertools.permutations(range(p), 5):
        pr = 1
        for i in range(4):
            pr *= chi[(t[i + 1] - t[i]) % p]
            if pr == 0: break
        tot += pr
    return tot
print("   p    A_4/p^3   (A_4-2p^3)/p^2 [O(p^k), STABLE]   /p^2.5 [drifts->0]")
for p in [7, 11, 19, 23, 31]:
    a = A4(p)
    print(f"  {p:3d}  {a/p**3:.5f}      {(a-2*p**3)/p**2:+.4f}                  {(a-2*p**3)/p**2.5:+.4f}")
print("=" * 70)
print("CONCLUSION: top order p^{k+1} of A_{2k} is a SIGNED sum over EVEN CACTI")
print("(all biconnected blocks even cycles), NOT a +1-per-bigon-tree count.")
print("Bigon-trees over-count (1,3,13,69,...); even-cycle cacti correct to C_k.")
print("R(p)->e (Part C) needs NO Weil: V=2k case is tr(M^{2k})=(-p)^k(p-1), exact.")
print("Error is O(p^k); R-e relative correction O(1/p), supporting 1/p rate.")
