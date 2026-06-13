#!/usr/bin/env python3
"""
paley_cluster_topterm_monad.py
monad-explorer-2026-06-07 (deep-research, 4th session)

SETTLE the leading coefficient c0 = lim A_{2k}/p^{k+1}  EXACTLY, by summing the
top-order density of EVERY connected coincidence pattern (not just even cacti):

    c0 = (-1)^k * sum_{rho connected} mu(0,rho) * g(rho),
    g(rho) := lim_p F(rho)/p^m,   m = cycle rank = 2k+1-V.

A pattern contributes iff g(rho) != 0.  For even cacti g=1 (PROVED: F=(p-1)^c).
The QUESTION: do NON-cactus even patterns (even theta graphs, etc.) have g!=0?
If yes, the prior "only even cacti" accounting (which gave 6,25,132) is incomplete,
and the extra patterns correct it back to Catalan 5,14,42 -- OR the true c0 is NOT
Catalan at all.  We measure g(rho) at two primes to detect nonzero limits and to
get the exact rational, then assemble c0 and compare to C_k.
"""
import math
from collections import defaultdict
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


def connected(edges, nb):
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
    return len(seen) == nb


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


def nice_rational(x, maxden=24):
    from fractions import Fraction
    return Fraction(x).limit_denominator(maxden)


from fractions import Fraction


def richardson_g(f1, m, p1, f2, p2):
    """g = lim F/p^m, eliminating the O(1/p) term via two-point Richardson.
       For perfect-square (top-order) patterns F/p^m = g(1+O(1/p)); for genuine
       character sums F/p^m -> 0. Returns float estimate."""
    r1 = f1 / p1 ** m
    r2 = f2 / p2 ** m
    return (p2 * r2 - p1 * r1) / (p2 - p1)


def round_nice(x):
    cands = [Fraction(a, b) for b in (1, 2, 3, 4, 6) for a in range(-2 * b, 2 * b + 1)]
    best = min(cands, key=lambda q: abs(float(q) - x))
    return best if abs(float(best) - x) < 0.03 else None


for k in [2, 3, 4]:
    L = 2 * k
    print("=" * 72)
    print(f"k={k}: assembling c0 = lim A_{2*k}/p^{k+1} from per-pattern top densities g(rho)")
    p1, p2 = {2: (37, 53), 3: (31, 43), 4: (19, 23)}[k]
    c0 = Fraction(0)
    contrib_by_type = defaultdict(lambda: Fraction(0))
    cnt_by_type = defaultdict(int)
    ncontrib = 0
    nbad = 0
    for blocks in set_partitions(range(L + 1)):
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges):
            continue
        if not connected(edges, nb):
            continue
        mu = mu_partition(blocks)
        if mu == 0:
            continue
        # contributing (perfect-square) patterns have edges in even series classes
        # => cycle rank m = E - V + 1 <= k.  m>k patterns have g=0 (skip; infeasible
        # to brute-force anyway, e.g. 2 vertices/6 parallel edges has m=5).
        m_cyc = L - nb + 1
        if m_cyc > k:
            continue
        f1, m1 = flow_sum(edges, nb, p1)
        f2, m2 = flow_sum(edges, nb, p2)
        m = m1
        gest = richardson_g(f1, m, p1, f2, m and p2)
        if abs(gest) < 0.04:
            continue                  # g -> 0
        gr = round_nice(gest)
        if gr is None:
            nbad += 1
            print(f"   !! unresolved g={gest:.4f} blocks={blocks} V={nb} m={m}")
            continue
        V = nb
        # is it an even cactus? (cycle rank == #biconnected cycle comps, all even & e==v)
        c0 += Fraction((-1) ** k) * mu * gr
        ncontrib += 1
        key = (V, m, float(gr))
        contrib_by_type[key] += Fraction((-1) ** k) * mu * gr
        cnt_by_type[key] += 1
    print(f"   g estimated at primes {p1},{p2} (Richardson 1/p); "
          f"#patterns with g!=0: {ncontrib}; unresolved: {nbad}")
    print(f"   c0 = {c0} = {float(c0):.6f}      Catalan C_{k} = {catalan(k)}   "
          f"MATCH={c0 == catalan(k)}")
    print(f"   breakdown by (V, cyclerank m, g) [#patterns]:")
    for key in sorted(contrib_by_type):
        V, m, g = key
        print(f"      V={V} m={m} g={g:+.3f}  [{cnt_by_type[key]:3d} patt]  "
              f"signed contrib = {float(contrib_by_type[key]):+.3f}")
print("=" * 72)
print("If c0=C_k and the breakdown shows g!=0 for NON-cactus even patterns (e.g. m=2")
print("with V<k+1 'theta'), then the prior 'only even cacti' top-order census was")
print("INCOMPLETE -- the extra even (non-cactus) patterns correct 6,25,132 -> Catalan.")
