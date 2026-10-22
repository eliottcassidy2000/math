"""Test of the pre-registered prediction E(13) = 1 (HYP-9356, committed in fad6a40581) by exact ears
(session opus-2026-10-21-S50).  Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md.

Completeness (as in strong_ears.py): by Moon's strong-deletion theorem and the insertion lemma (THM-1370), every strong
m-tournament with H <= X is a nonconstant one-vertex ear of a strong (m-1)-tournament with H <= X.  Starting from the
census (n <= 9) this builds every strong class with H <= 625 at orders 10, 11, 12, and every strong 13-class with
H = 625 = f(13) (Busch, CITED), checking that no strong 13-tournament has H < 625.  Canonical forms: the S36
refined-invariant code for n <= 11, an exact individualisation-refinement code (Python) for n = 12, 13.
Reports E(13) and whether the extremal class is the chain IC(3,3,3,3), and the strong polynomials P_d for odd
375 <= d < 625 at orders <= 12 (complete for d < 625 since f(13) = 625).
Prints ALL CHECKS PASSED.  Runtime about 40 minutes."""
import sys, os, time, itertools
from math import factorial
from collections import defaultdict, Counter
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 's36_audit_20261014'))
import h_multiplicities as HM                                            # noqa: E402
from audit_tournaments_classes import pair_index, code_to_adj           # noqa: E402
from audit_tournaments_classes_n9 import canon_refined, vertex_invariants  # noqa: E402

FAIL = []
T0 = time.time()
LIMIT = 625


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def busch(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]


fact = np.array([factorial(k) for k in range(10)], dtype=np.int64)
permoff = np.zeros(10, dtype=np.int64); chunks = []; off = 0
for k in range(10):
    permoff[k] = off
    P_ = np.array(list(itertools.permutations(range(k))), dtype=np.int64).reshape(-1) if k > 0 else np.zeros(0, dtype=np.int64)
    chunks.append(P_); off += len(P_)
permtab = np.concatenate(chunks)


@nb.njit(cache=True)
def ham_buf(n, adj, dp):
    N = 1 << n
    for S in range(N):
        for v in range(n):
            dp[S, v] = 0
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            for u in range(n):
                if (S >> u) & 1:
                    continue
                if (adj[v] >> u) & 1:
                    dp[S | (1 << u), u] += c
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


@nb.njit(parallel=True, cache=True)
def ears(n, parents, lim):
    """H of every nonconstant ear (or -1 when the partial count already exceeds lim is NOT used: exact DP)."""
    p = parents.shape[0]
    full = (1 << (n - 1)) - 1
    Hs = np.zeros((p, 1 << (n - 1)), dtype=np.int64)
    for t in nb.prange(p):
        dp = np.zeros((1 << n, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64)
        for S in range(1, full):
            for v in range(n - 1):
                a = parents[t, v]
                if not (S >> v) & 1:
                    a |= 1 << (n - 1)
                adj[v] = a
            adj[n - 1] = S
            Hs[t, S] = ham_buf(n, adj, dp)
    return Hs


@nb.njit(parallel=True, cache=True)
def max_block(n, cand):
    m = cand.shape[0]
    out = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        inv = np.empty(n, dtype=np.int64)
        vertex_invariants(n, cand[t], inv)
        s = np.sort(inv)
        best = 1; run = 1
        for i in range(1, n):
            run = run + 1 if s[i] == s[i - 1] else 1
            best = max(best, run)
        out[t] = best
    return out


def canon_py(adj):
    n = len(adj)
    outn = [[(adj[u] >> v) & 1 for v in range(n)] for u in range(n)]

    def refine(col):
        while True:
            sig = [(col[u], tuple(sorted(col[v] for v in range(n) if outn[u][v])),
                    tuple(sorted(col[v] for v in range(n) if v != u and not outn[u][v]))) for u in range(n)]
            keys = sorted(set(sig))
            new = [keys.index(s) for s in sig]
            if len(set(new)) == len(set(col)):
                return new
            col = new

    best = [None]

    def search(col):
        col = refine(col)
        if len(set(col)) == n:
            order = sorted(range(n), key=lambda v: col[v])
            c = 0; e = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if outn[order[i]][order[j]]:
                        c |= 1 << e
                    e += 1
            if best[0] is None or c < best[0]:
                best[0] = c
            return
        cells = defaultdict(list)
        for v in range(n):
            cells[col[v]].append(v)
        target = min((len(vs), k) for k, vs in cells.items() if len(vs) > 1)[1]
        for v in cells[target]:
            search([2 * c + (1 if c > target or (c == target and u != v) else 0) for u, c in enumerate(col)])

    search([bin(a).count('1') for a in adj])
    return best[0]


def code_to_adj_big(n, code):
    adj = [0] * n
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> e) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            e += 1
    return adj


def level(n, parents, lim, exact=None):
    """Children (ears) with H <= lim (or == exact): dict canonical code -> H.  Also returns the least ear value."""
    res = {}
    least = None
    B = 256
    for b0 in range(0, len(parents), B):
        par = np.array(parents[b0:b0 + B], dtype=np.int64)
        Hs = ears(n, par, lim)
        inner = Hs[:, 1:-1]
        mn = int(inner.min())
        least = mn if least is None else min(least, mn)
        sel = np.argwhere((Hs > 0) & ((Hs == exact) if exact else (Hs <= lim)))
        if len(sel) == 0:
            continue
        cand = np.zeros((len(sel), n), dtype=np.int64)
        hv = np.zeros(len(sel), dtype=np.int64)
        for r, (t, S) in enumerate(sel):
            for v in range(n - 1):
                a = par[t, v]
                if not (S >> v) & 1:
                    a |= 1 << (n - 1)
                cand[r, v] = a
            cand[r, n - 1] = S
            hv[r] = Hs[t, S]
        if n <= 11:
            assert max_block(n, cand).max() <= 9
            codes, _ = canon_refined(n, cand, pair_index(n), permtab, permoff, fact)
            for c, h in zip(codes, hv):
                res.setdefault(int(c), int(h))
        else:
            for row, h in zip(cand, hv):
                res.setdefault(canon_py([int(x) for x in row]), int(h))
        print(f'       order {n}: parents {b0 + len(par)}/{len(parents)}, classes so far {len(res)} [{time.time() - T0:.0f}s]', flush=True)
    return res, least


def chain(c):
    n = 1 + sum(c)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    l = 0
    for p in c:
        r = l + p
        adj[l] &= ~(1 << r); adj[r] |= 1 << l
        l = r
    return adj


data = HM.classes_upto(9)
A9 = data[9]['A']
H9, st9 = HM.ham_all(9, A9)
par = [A9[i].copy() for i in range(len(H9)) if st9[i] and H9[i] <= LIMIT]
inv = {m: Counter() for m in range(3, 13)}
for m in range(3, 10):
    A = data[m]['A']; H, st = HM.ham_all(m, A)
    for h, s in zip(H, st):
        if s and h < LIMIT:
            inv[m][int(h)] += 1
print(f'census loaded; strong 9-classes with H <= {LIMIT}: {len(par)} [{time.time() - T0:.0f}s]', flush=True)
for n in (10, 11, 12):
    res, least = level(n, par, LIMIT)
    check(least == busch(n), f'order {n}: least ear value {least} = f({n}) = {busch(n)}')
    print(f'     order {n}: {len(res)} strong classes with H <= {LIMIT} [{time.time() - T0:.0f}s]', flush=True)
    for c, h in res.items():
        if h < LIMIT:
            inv[n][h] += 1
    par = [code_to_adj(n, c) for c in res] if n <= 11 else [code_to_adj_big(n, c) for c in res]
res13, least13 = level(13, par, LIMIT, exact=LIMIT)
check(least13 == 625, f'no strong 13-tournament has H < 625 (least ear value {least13})')
E13 = len(res13)
ch = canon_py(chain((3, 3, 3, 3)))
print(f'     E(13) = {E13}; the extremal class is IC(3,3,3,3): {set(res13) == {ch}}')
check(E13 >= 1 and ch in res13, 'IC(3,3,3,3) is extremal at order 13 (THM-4744 (1))')
print(f'     PRE-REGISTERED PREDICTION E(13) = 1: {"PASSED" if E13 == 1 and set(res13) == {ch} else "FAILED"}')
zeros = [d for d in range(1, LIMIT, 2) if not any(inv[m][d] for m in inv)]
check(zeros == [1, 7, 21], f'odd d < 625 with no strong tournament of any order: {zeros}')
print('     P_d(1) for odd d = 375..623:', [sum(inv[m][d] for m in inv) for d in range(375, LIMIT, 2)])

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
