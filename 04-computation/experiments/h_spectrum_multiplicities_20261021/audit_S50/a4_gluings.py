"""Audit S50 - a4: THM-4745 checked with independent code.
  (G1) exact-gluing census: for chain pairs (A, B), every shared-vertex position (i in A, j in B) and every orientation
       of the (|A|-1)(|B|-1) cross arcs: H(T) >= H(A)H(B) always; list the exact gluings; check they are one-way and
       end-to-start (any vertex of a C_3 factor).  Pairs: all 28 ordered pairs of chains with parts in {2,3,4,5},
       <= 6 vertices and (|A|-1)(|B|-1) <= 16, plus C_3 / R_4 with IC(3,3) both ways (the session's 32), plus extra
       pairs outside the session's range (at most 18 cross bits): C_3 with IC(2,2,2) both ways, IC(4,2), IC(2,3,2), IC(6),
       IC(2,4); R_4 with IC(2,2,2).
  (G2) failure modes for C_3 u C_3, C_3 u R_4, R_4 u C_3, R_4 u R_4 by explicit path enumeration.
  (G3) the bad-split constructions (a), (b), (c) of the proof of THM-4745 (2), built literally as written, for every
       chain with N = 3..10 and every spine position t < N; and absence of any bad split at t = N (exhaustive, N <= 7)."""
import time, itertools
from collections import Counter
import numpy as np
import numba as nb
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


@nb.njit(parallel=True, cache=True)
def glue_scan(base, a, b, out):
    """base: adjacency of T without cross arcs (A' = 0..a-2, v = a-1, B' = a..a+b-2). For every mask over the
    (a-1)(b-1) cross pairs (bit x*(b-1)+y set => x -> a+y) store H(T)."""
    n = a + b - 1
    nb_ = (a - 1) * (b - 1)
    tot = 1 << nb_
    for mask in nb.prange(tot):
        adj = base.copy()
        for x in range(a - 1):
            for y in range(b - 1):
                if (mask >> (x * (b - 1) + y)) & 1:
                    adj[x] |= 1 << (a + y)
                else:
                    adj[a + y] |= 1 << x
        out[mask] = L.ham_dp(n, adj)


def glue_base(A, B, i, j):
    a, b = len(A), len(B)
    n = a + b - 1
    ma = {}; k = 0
    for x in range(a):
        if x == i:
            ma[x] = a - 1
        else:
            ma[x] = k; k += 1
    mb = {}; k = a
    for y in range(b):
        if y == j:
            mb[y] = a - 1
        else:
            mb[y] = k; k += 1
    base = [0] * n
    for x in range(a):
        for x2 in range(a):
            if (A[x] >> x2) & 1:
                base[ma[x]] |= 1 << ma[x2]
    for y in range(b):
        for y2 in range(b):
            if (B[y] >> y2) & 1:
                base[mb[y]] |= 1 << mb[y2]
    return np.array(base, dtype=np.int64)


# ------------------------------------------------------------------ (G1)
chains = [c for N in range(2, 6) for c in L.compositions(N, 2, 5)]
pairs = [(cA, cB) for cA in chains for cB in chains if sum(cA) * sum(cB) <= 16]
session = pairs + [((2,), (3, 3)), ((3, 3), (2,)), ((3,), (3, 3)), ((3, 3), (3,))]
extra = [((2,), (2, 2, 2)), ((2, 2, 2), (2,)), ((2,), (4, 2)), ((2,), (2, 3, 2)), ((3,), (2, 2, 2)),
         ((2,), (6,)), ((6,), (2,)), ((2,), (2, 4))]
allpairs = session + extra
n_exact = 0; below = 0; bad = []; seen_both = True; per_pair = []
for (cA, cB) in allpairs:
    A, B = L.chain(cA), L.chain(cB)
    a, b = len(A), len(B)
    target = L.H(A) * L.H(B)
    full = (1 << ((a - 1) * (b - 1))) - 1
    ex_here = 0; ab = ba = False
    for i in range(a):
        for j in range(b):
            out = np.zeros(full + 1, dtype=np.int64)
            glue_scan(glue_base(A, B, i, j), a, b, out)
            below += int(np.sum(out < target))
            for m in np.nonzero(out == target)[0]:
                m = int(m); ex_here += 1
                if m == full:
                    ok = (i == a - 1 or cA == (2,)) and (j == 0 or cB == (2,))
                    ab |= ok
                elif m == 0:
                    ok = (j == b - 1 or cB == (2,)) and (i == 0 or cA == (2,))
                    ba |= ok
                else:
                    ok = False
                if not ok:
                    bad.append((cA, cB, i, j, m))
    expect = 2 * (3 if cA == (2,) else 1) * (3 if cB == (2,) else 1)
    if ex_here != expect:
        bad.append(('count', cA, cB, ex_here, expect))
    seen_both &= ab and ba
    n_exact += ex_here if (cA, cB) in session else 0
    per_pair.append(ex_here)
print(f'     {len(session)} session pairs + {len(extra)} extra pairs; exact gluings in the session range: {n_exact} [{time.time() - T0:.0f}s]')
check(below == 0, "(G1) H(T) >= H(A)H(B) for every gluing (all positions, all cross orientations)")
check(not bad and seen_both, '(G1) every exact gluing is one-way and end-to-start (C_3: any vertex); exact count per pair = 2 x (3 per C_3 factor); '
      f'session total {n_exact} (session: 136)' + (f'  BAD {bad[:6]}' if bad else ''))


# ------------------------------------------------------------------ (G2)
def paths_of(adj, verts):
    verts = list(verts)
    out = []
    S = set(verts)

    def dfs(p, rem):
        if not rem:
            out.append(tuple(p)); return
        u = p[-1]
        for w in list(rem):
            if (adj[u] >> w) & 1:
                p.append(w); rem.remove(w)
                dfs(p, rem)
                rem.add(w); p.pop()
    for s in verts:
        dfs([s], S - {s})
    return out


def is_path(adj, seq):
    return all((adj[seq[t]] >> seq[t + 1]) & 1 for t in range(len(seq) - 1))


rows = []; ok_c33 = True; ok_r4 = True; least = {}
for cA, cB in [((2,), (2,)), ((2,), (3,)), ((3,), (2,)), ((3,), (3,))]:
    A, B = L.chain(cA), L.chain(cB)
    a, b = len(A), len(B)
    hAB = L.H(A) * L.H(B)
    for i in range(a):
        for j in range(b):
            base = [int(x) for x in glue_base(A, B, i, j)]
            v = a - 1
            Av = list(range(a)); Bv = [v] + list(range(a, a + b - 1))
            end_pos = ((i == a - 1 or cA == (2,)) and (j == 0 or cB == (2,))) or ((j == b - 1 or cB == (2,)) and (i == 0 or cA == (2,)))
            full = (1 << ((a - 1) * (b - 1))) - 1
            only_U = 0; mixed_min = None; allR = True
            for mask in range(1, full):
                adj = list(base)
                for x in range(a - 1):
                    for y in range(b - 1):
                        if (mask >> (x * (b - 1) + y)) & 1:
                            adj[x] |= 1 << (a + y)
                        else:
                            adj[a + y] |= 1 << x
                P = paths_of(adj, range(a + b - 1))
                rv = 0; pc = Counter()
                for p in P:
                    pa = tuple(x for x in p if x in Av); pb = tuple(x for x in p if x in Bv)
                    if is_path(adj, pa) and is_path(adj, pb):
                        pc[(pa, pb)] += 1
                    else:
                        rv += 1
                extra_m = sum(c - 1 for c in pc.values())
                assert len(P) - hAB == rv + extra_m          # the counting identity of THM-4745 (1)
                assert len(pc) == hAB                         # every pair has >= 1 merge
                if rv == 0 and extra_m > 0:
                    only_U += 1
                if rv == 0:
                    allR = False
                s = len(P) - hAB
                mixed_min = s if mixed_min is None else min(mixed_min, s)
            if (cA, cB) == ((2,), (2,)):
                ok_c33 &= only_U == 4
            elif end_pos:
                ok_r4 &= allR
                least[(cA, cB)] = min(least.get((cA, cB), 10 ** 9), mixed_min)
check(ok_c33, '(G2) C_3 u C_3: at every position exactly 4 of the 14 mixed patterns violate only (U) (no path breaks (R))')
check(ok_r4, '(G2) with an R_4 factor, at the end-to-start positions every mixed pattern has a Hamiltonian path breaking (R)')
check(least.get(((2,), (3,))) == 4 and least.get(((3,), (2,))) == 4 and least.get(((3,), (3,))) == 6,
      f'(G2) least surplus of a mixed pattern at end-to-start positions: {least} (session: 4 for C_3 with R_4, 6 for R_4 with R_4)')
check(True, '(G2) for every pattern: H(T) - H(A)H(B) = #(R-violating paths) + #(surplus merges), and every pair has a merge (asserted)')


# ------------------------------------------------------------------ (G3)
def ham_path_ending(adj, verts, end):
    """a Hamiltonian path of adj[verts] ending at `end` (DFS), or None"""
    verts = set(verts)
    def dfs(p, rem):
        if not rem:
            return list(p) if p[-1] == end else None
        u = p[-1]
        for w in sorted(rem):
            if (adj[u] >> w) & 1 and (w != end or len(rem) == 1):
                p.append(w); rem.remove(w)
                r = dfs(p, rem)
                if r:
                    return r
                rem.add(w); p.pop()
        return None
    for s in sorted(verts):
        if s == end and len(verts) > 1:
            continue
        r = dfs([s], verts - {s})
        if r:
            return r
    return None


def bad_split_ok(adj, v, Ap, Xa, Ya):
    if not Xa or set(Xa) & set(Ya) or set(Xa) | set(Ya) != set(Ap) or len(Xa) + len(Ya) != len(Ap):
        return False
    if not is_path(adj, Xa) or not is_path(adj, Ya):
        return False
    if not (adj[v] >> Xa[-1]) & 1:
        return False
    if Ya and not (adj[v] >> Ya[0]) & 1:
        return False
    return True


ok_cases = True; ncase = Counter(); fails = []
for N in range(3, 11):
    for p in L.compositions(N):
        adj = L.chain(p)
        k = len(p)
        ls = [0]
        for q in p:
            ls.append(ls[-1] + q)
        starts = ls[:-1]
        for t in range(N):
            v = t
            Ap = [u for u in range(N + 1) if u != v]
            blk = max(jj for jj in range(k) if ls[jj] <= t)       # block index containing t as start or interior (0-based)
            interior = t not in ls
            if interior:
                case = 'a'
                assert L.is_strong([0] * 0) if False else True
                sub = [adj[u] for u in Ap]
                # A - v strong (as claimed): check
                idx = {u: r for r, u in enumerate(Ap)}
                sub = [sum(1 << idx[w] for w in Ap if (adj[u] >> w) & 1) for u in Ap]
                strong = L.is_strong(sub)
                outn = [u for u in Ap if (adj[v] >> u) & 1]
                Xa = ham_path_ending(adj, Ap, outn[0]) if strong else None
                Ya = []
                good = strong and Xa is not None and bad_split_ok(adj, v, Ap, Xa, Ya)
            elif t in starts and starts.index(t) < k - 1:
                case = 'b'
                j = starts.index(t)
                tail = list(range(ls[j + 1], N + 1))
                rest = [u for u in Ap if u not in tail]
                into = all((adj[x] >> y) & 1 for x in rest for y in tail)
                Xa = ham_path_ending(adj, Ap, N)
                good = into and (adj[v] >> N) & 1 and Xa is not None and bad_split_ok(adj, v, Ap, Xa, [])
            else:
                case = 'c'
                lk = ls[k - 1]
                assert t == lk
                if p[-1] >= 3:
                    Xa = list(range(0, lk)) + [lk + 1]
                    Ya = list(range(lk + 2, N + 1))
                else:
                    assert k >= 2
                    lk1 = ls[k - 2]
                    c = lk + 1
                    Xa = list(range(0, lk1)) + [c]
                    Ya = list(range(lk1, lk)) + [N]
                good = bad_split_ok(adj, v, Ap, Xa, Ya)
            ncase[case] += 1
            if not good:
                ok_cases = False; fails.append((p, t, case))
check(ok_cases, f'(G3) the proof\'s bad split exists at every position t < N, built as written, for every chain with N = 3..10: '
      f'cases a/b/c used {dict(ncase)}' + (f'  FAILS {fails[:5]}' if fails else ''))

# no bad split at t = N (exhaustive), N <= 7
ok_none = True
for N in range(3, 8):
    for p in L.compositions(N):
        adj = L.chain(p)
        v = N
        Ap = list(range(N))
        for r in range(1, N + 1):
            for Xs in itertools.combinations(Ap, r):
                Ys = [u for u in Ap if u not in Xs]
                for Xa in itertools.permutations(Xs):
                    if not is_path(adj, Xa) or not (adj[v] >> Xa[-1]) & 1:
                        continue
                    for Ya in (itertools.permutations(Ys) if Ys else [()]):
                        if is_path(adj, Ya) and (not Ya or (adj[v] >> Ya[0]) & 1):
                            ok_none = False
check(ok_none, '(G3) at t = N (the last spine vertex) no bad split exists (exhaustive over all ordered pairs of paths, N <= 7)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
