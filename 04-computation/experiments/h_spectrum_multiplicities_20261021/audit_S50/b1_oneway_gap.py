"""Audit S50, second pass - b1: coverage of the written proof of THM-4745 (5) (one-way lemma, parts in {2,3,4}).

The written proof: Step 1 excludes junctions; Step 2 excludes positions interior to a block of >= 4 vertices
("a 3-cycle counts as an end"); Step 3 (finite) says the local gluing T_loc of the blocks containing v is exact only
when one-way, forward with v the end of A_loc and the start of B_loc, or backward (a C_3 block: any vertex); Step 4
then assumes v = x_N = y_0 (forward) or v = x_0 = y_M (backward) GLOBALLY.

(G1) Enumerate every chain pair (parts in {2,3,4}, |A|, |B| <= 9) and every position (i, j), and list the cases that
     survive Steps 1-3 but are not of the form Step 4 assumes ("uncovered").  Types:
       MID_A / MID_B : v is the middle vertex of a C_3 block of a chain with >= 2 blocks;
       END_A / END_B : v is the 'wrong' end of a terminal C_3 block for the local orientation
                       (e.g. v = x_0, first block C_3, local arcs forward).
(G2) For every uncovered case, search for a certificate: a Hamiltonian path that uses only arcs of A, of B and the
     LOCAL cross arcs (oriented as Step 3 forces) and whose restriction to A or to B is not a path.  Such a path is a
     Hamiltonian path of every completion, so no completion is exact (THM-4745 (1)).
(G3) For uncovered cases with <= 14 free (non-local) cross arcs and |T| <= 12, check every completion by exact DP:
     none is exact.
(G4) Step 4's identities for every chain with parts in {2,3,4}, N <= 14, k >= 2:
       Init(A - x_N) = {x_0..x_(l_k)}, N^-(x_N) = A' - {x_(l_k)};  Term(B - y_0) = {y_(m_2)..y_M},
       N^+(y_0) = B' - {y_(m_2)}; and the arcs not forced by (B1), (B2) are exactly the local ones.
(G5) The Init/Term lemma of (T3), by brute force on random tournaments with n <= 8."""
import random, time, itertools
from collections import Counter, defaultdict
import numpy as np
import numba as nb
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def bounds(parts):
    ls = [0]
    for p in parts:
        ls.append(ls[-1] + p)
    return ls


def pos_info(parts, i):
    """('J',) junction; ('I', j) interior of a block with >= 4 vertices; ('E', j, s, e) a position in block j
    (start s, end e) that is an end of the block or any non-junction vertex of a C_3 block."""
    ls = bounds(parts)
    if i in ls[1:-1]:
        return ('J',)
    j = max(t for t in range(len(parts)) if ls[t] <= i)
    if ls[j] < i < ls[j + 1] and parts[j] >= 3:
        return ('I', j)
    return ('E', j, ls[j], ls[j + 1])


@nb.njit(parallel=True, cache=True)
def scan(base, free_u, free_w, n, out):
    nf = free_u.shape[0]
    for mask in nb.prange(1 << nf):
        adj = base.copy()
        for f in range(nf):
            if (mask >> f) & 1:
                adj[free_u[f]] |= 1 << free_w[f]
            else:
                adj[free_w[f]] |= 1 << free_u[f]
        out[mask] = L.ham_dp(n, adj)


def find_certificate(adj, n, Aset, Bset):
    """DFS for a Hamiltonian path of the partial digraph adj whose restriction to A or B is not a path."""
    full = (1 << n) - 1
    res = [None]

    def ok_restr(path, S):
        seq = [x for x in path if x in S]
        return all((adj[seq[t]] >> seq[t + 1]) & 1 for t in range(len(seq) - 1))

    def dfs(path, used):
        if res[0] is not None:
            return
        if used == full:
            if not ok_restr(path, Aset) or not ok_restr(path, Bset):
                res[0] = list(path)
            return
        nxt = adj[path[-1]] & ~used
        while nxt:
            low = nxt & -nxt
            w = low.bit_length() - 1
            path.append(w)
            dfs(path, used | low)
            path.pop()
            nxt ^= low
            if res[0] is not None:
                return

    for s in range(n):
        dfs([s], 1 << s)
        if res[0] is not None:
            break
    return res[0]


comps = [c for N in range(2, 9) for c in L.compositions(N, 2, 4)]
NMAX_CERT = 13
cases = []          # uncovered: (pa, pb, i, j, orient, types)
n_pos = 0; n_cov = Counter()
for pa in comps:
    for pb in comps:
        NA, NB = sum(pa), sum(pb)
        for i in range(NA + 1):
            for j in range(NB + 1):
                n_pos += 1
                ia, ib = pos_info(pa, i), pos_info(pb, j)
                if ia[0] == 'J' or ib[0] == 'J':
                    n_cov['step1'] += 1; continue
                if ia[0] == 'I' or ib[0] == 'I':
                    n_cov['step2'] += 1; continue
                _, ja, sa, ea = ia
                _, jb, sb, eb = ib
                c3a, c3b = pa[ja] == 2, pb[jb] == 2
                wholeA, wholeB = pa == (2,), pb == (2,)
                for o in ('fwd', 'bwd'):
                    if o == 'fwd':
                        allowed = (c3a or i == ea) and (c3b or j == sb)
                        covered = (wholeA or i == NA) and (wholeB or j == 0)
                    else:
                        allowed = (c3a or i == sa) and (c3b or j == eb)
                        covered = (wholeA or i == 0) and (wholeB or j == NB)
                    if not allowed:
                        n_cov['step3 (local gluing not exact)'] += 1
                        continue
                    if covered:
                        n_cov['step4'] += 1
                        continue
                    types = []
                    if not wholeA:
                        if i not in (0, NA):
                            types.append('MID_A')
                        elif (o == 'fwd' and i != NA) or (o == 'bwd' and i != 0):
                            types.append('END_A')
                    if not wholeB:
                        if j not in (0, NB):
                            types.append('MID_B')
                        elif (o == 'fwd' and j != 0) or (o == 'bwd' and j != NB):
                            types.append('END_B')
                    assert types
                    cases.append((pa, pb, i, j, o, tuple(types)))
tc = Counter(t for c in cases for t in c[5])
print(f'     {len(comps)} chains, {n_pos} (A, B, i, j) positions; covered: {dict(n_cov)}')
print(f'     UNCOVERED (survive Steps 1-3, not of the form assumed in Step 4): {len(cases)} (case, orientation) pairs; '
      f'types {dict(tc)}')
print('     smallest examples:', [(c[0], c[1], f'v=x_{c[2]}', f'v=y_{c[3]}', c[4], c[5]) for c in sorted(cases, key=lambda c: sum(c[0]) + sum(c[1]))[:6]])
check(len(cases) > 0, '(G1) the written proof of THM-4745 (5) does not reach these configurations (C_3 blocks at non-end '
      'positions, or the wrong end of a terminal C_3 block); Step 3 allows them locally and Step 4 assumes v = x_N = y_0')


# ------------------------------------------------------------------ (G2), (G3)
def glue(pa, pb, i, j, orient):
    A, B = L.chain(pa), L.chain(pb)
    a, b = len(A), len(B)
    n = a + b - 1
    mA = {}; k = 0
    for x in range(a):
        if x == i:
            mA[x] = a - 1
        else:
            mA[x] = k; k += 1
    mB = {}; k = a
    for y in range(b):
        if y == j:
            mB[y] = a - 1
        else:
            mB[y] = k; k += 1
    base = [0] * n
    for x in range(a):
        for x2 in range(a):
            if (A[x] >> x2) & 1:
                base[mA[x]] |= 1 << mA[x2]
    for y in range(b):
        for y2 in range(b):
            if (B[y] >> y2) & 1:
                base[mB[y]] |= 1 << mB[y2]
    ia, ib = pos_info(pa, i), pos_info(pb, j)
    locA = [mA[x] for x in range(ia[2], ia[3] + 1) if x != i]
    locB = [mB[y] for y in range(ib[2], ib[3] + 1) if y != j]
    for x in locA:
        for y in locB:
            if orient == 'fwd':
                base[x] |= 1 << y
            else:
                base[y] |= 1 << x
    free = [(x, y) for x in range(a - 1) for y in range(a, n) if not (x in locA and y in locB)]
    Aset = set(range(a)); Bset = set([a - 1] + list(range(a, n)))
    return base, free, n, Aset, Bset, L.H(A) * L.H(B)


ncert = 0; nocert = []; nexh = 0; exact_found = []
for (pa, pb, i, j, o, types) in cases:
    base, free, n, Aset, Bset, target = glue(pa, pb, i, j, o)
    if n <= NMAX_CERT:
        cert = find_certificate(base, n, Aset, Bset)
        if cert is None:
            nocert.append((pa, pb, i, j, o, types))
        else:
            ncert += 1
    if len(free) <= 14 and n <= 12:
        out = np.zeros(1 << len(free), dtype=np.int64)
        scan(np.array(base, dtype=np.int64), np.array([f[0] for f in free], dtype=np.int64),
             np.array([f[1] for f in free], dtype=np.int64), n, out)
        nexh += 1
        if np.any(out == target):
            exact_found.append((pa, pb, i, j, o))
        assert np.all(out >= target)
n15 = sum(1 for c in cases if sum(c[0]) + sum(c[1]) + 1 <= NMAX_CERT)
check(not nocert, f'(G2) every uncovered case with |T| <= {NMAX_CERT} ({n15} cases) has a certificate: a Hamiltonian path using only '
      f'arcs of A, B and the local cross arcs whose restriction is not a path ({ncert} found)' + (f'; NONE for {nocert[:4]}' if nocert else ''))
check(not exact_found, f'(G3) exhaustive completions: {nexh} uncovered cases (<= 14 free arcs, |T| <= 12) have no exact completion'
      + (f'; EXACT {exact_found[:4]}' if exact_found else ''))


# ------------------------------------------------------------------ (G4)
def components(adj, verts):
    """strong components of the subtournament on verts, in order (initial first)."""
    verts = list(verts)
    idx = {v: r for r, v in enumerate(verts)}
    reach = {u: {u} for u in verts}
    for u in verts:
        st = [u]
        while st:
            x = st.pop()
            for w in verts:
                if (adj[x] >> w) & 1 and w not in reach[u]:
                    reach[u].add(w); st.append(w)
    comps_, seen = [], set()
    for u in verts:
        if u in seen:
            continue
        c = {w for w in verts if w in reach[u] and u in reach[w]}
        seen |= c
        comps_.append(c)
    # order: initial component reaches all others
    comps_.sort(key=lambda c: -len(set().union(*[reach[x] for x in c])))
    return comps_


ok4 = True; n4 = 0
for N in range(4, 15):
    for p in L.compositions(N, 2, 4):
        if len(p) < 2:
            continue
        n4 += 1
        adj = L.chain(p); ls = bounds(p); lk = ls[-2]; m2 = ls[1]
        # A side: v = x_N
        cs = components(adj, range(N))
        init = cs[0]
        nin = {x for x in range(N) if (adj[x] >> N) & 1}
        ok4 &= init == set(range(lk + 1)) and nin == set(range(N)) - {lk}
        # B side: v = y_0
        cs = components(adj, range(1, N + 1))
        term = cs[-1]
        nout = {y for y in range(1, N + 1) if (adj[0] >> y) & 1}
        ok4 &= term == set(range(m2, N + 1)) and nout == set(range(1, N + 1)) - {m2}
# partition of the cross arcs
okp = True
for pa in [c for N in range(2, 9) for c in L.compositions(N, 2, 4)]:
    for pb in [c for N in range(2, 9) for c in L.compositions(N, 2, 4)]:
        NA, NB = sum(pa), sum(pb)
        lsA, lsB = bounds(pa), bounds(pb)
        Aprime = set(range(NA)); Bprime = set(range(1, NB + 1))
        forcedA = set(range(lsA[-2])) if len(pa) >= 2 else set()            # {x_0..x_(l_k - 1)} x B'
        forcedB = set(range(lsB[1] + 1, NB + 1)) if len(pb) >= 2 else set()  # A' x {y_(m_2+1)..y_M}
        rest = {(x, y) for x in Aprime for y in Bprime if x not in forcedA and y not in forcedB}
        local = {(x, y) for x in range(lsA[-2], NA) for y in range(1, lsB[1] + 1)}
        okp &= rest == local
check(ok4 and okp, f'(G4) Step 4 identities hold for all {n4} chains with parts in {{2,3,4}}, N <= 14, k >= 2, and the arcs left '
      'after (B1), (B2) are exactly the local ones (all pairs with N <= 8)')

# ------------------------------------------------------------------ (G5)
rng = random.Random(5)
ok5 = True
for _ in range(400):
    n = rng.randint(2, 8)
    adj = [0] * n
    for x in range(n):
        for y in range(x + 1, n):
            if rng.random() < 0.5:
                adj[x] |= 1 << y
            else:
                adj[y] |= 1 << x
    cs = components(adj, range(n))
    ends = set(); starts = set()
    for perm in itertools.permutations(range(n)):
        if all((adj[perm[t]] >> perm[t + 1]) & 1 for t in range(n - 1)):
            starts.add(perm[0]); ends.add(perm[-1])
    ok5 &= ends == cs[-1] and starts == cs[0]
check(ok5, '(G5) in 400 random tournaments (n <= 8): the possible last (first) vertices of Hamiltonian paths are exactly the terminal (initial) strong component')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
