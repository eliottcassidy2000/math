"""Strong tournaments with few Hamiltonian paths beyond order 9, by exact ears (session opus-2026-10-21-S50).
Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md; THM-4744, HYP-9356.

Usage: python strong_ears.py [LIMIT]      LIMIT = 225 (default; orders 10, 11) or 375 (orders 10, 11, 12).

Completeness of the ear method (PROVED; both ingredients are in THM-1370):
  * Moon (1966): every strong tournament on m >= 4 vertices has a vertex v with T - v strong;
  * insertion lemma: H(T) >= H(T - v) for every vertex v.
So every strong m-tournament with H <= X is a nonconstant one-vertex ear of a strong (m-1)-tournament with H <= X.
From the exhaustive n <= 9 census (h_multiplicities.py) this script builds every strong class of order 10, 11 with
H <= LIMIT and every strong class of the next order with H = f (Busch's minimum, CITED: no strong tournament of that
order has H < f, re-checked here on the ears).

  (E0) PROVED + checked: an interval chain IC(p_1..p_k) (vertices 0..n-1, every arc forward except one back arc from
       the last to the first vertex of each block [l_i, l_i + p_i], consecutive blocks sharing an end vertex) has
       exactly prod (2^(p_i - 1) + 1) Hamiltonian paths; checked for every composition into parts 2..6, n <= 14.
  (E1) strong classes of order 10 (and 11) with H <= LIMIT; (E2) the extremal classes of the next order.
  (E3) the strong raw-count polynomials P_d, complete for every odd d < f(top order + 1).
  (E4) Busch-extremal classes E(m) against the composition count (parts 2, 3, 4 with optimal product), and the
       identification: the extremal classes are exactly the optimal interval chains, distinct compositions give distinct
       classes, self-converse iff palindromic, |Aut| = 1 (m >= 4), scores (1,1,2,...,m-2,m-2), m - 2 cyclic triples.
  (E5) exact multiplicity polynomials N(n,h) for every odd h below the completeness bound (inherited class series of
       spectrum_basins_20261020e/spectrum.md eq. (4)), cross-checked against the census for n <= 9.
Canonical forms: the S36 refined-invariant minimum code (canon_refined) for n <= 11 (55-bit codes; invariant blocks of
size <= 9 asserted); at n = 12 an exact Python individualisation-refinement minimum (big-integer codes).
Prints ALL CHECKS PASSED.  Runtime about 1 minute (LIMIT 225) or a few minutes (LIMIT 375)."""
import sys, os, time, itertools
from math import factorial, comb
from collections import defaultdict, Counter
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 's36_audit_20261014'))
import h_multiplicities as HM                                            # noqa: E402
from h_multiplicities import ham_one                                     # noqa: E402
from audit_tournaments_classes import pair_index, code_to_adj           # noqa: E402
from audit_tournaments_classes_n9 import canon_refined, vertex_invariants  # noqa: E402

FAIL = []
T0 = time.time()
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 225
assert LIMIT in (225, 375)
TOP = 11 if LIMIT == 225 else 12          # the order at which only H = f(TOP) = LIMIT is collected


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def busch(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]


def comps(total, allowed):
    if total == 0:
        yield ()
        return
    for p in allowed:
        if p <= total:
            for rest in comps(total - p, allowed):
                yield (p,) + rest


def prodH(c):
    out = 1
    for p in c:
        out *= 2 ** (p - 1) + 1
    return out


def optimal(m):
    return [c for c in comps(m - 1, (2, 3, 4)) if prodH(c) == busch(m)]


def chain(c):
    n = 1 + sum(c)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    l = 0
    for p in c:
        r = l + p
        adj[l] &= ~(1 << r)
        adj[r] |= 1 << l
        l = r
    return adj


def Hof(adj):
    return int(ham_one(len(adj), np.array(adj, dtype=np.int64)))


fact = np.array([factorial(k) for k in range(10)], dtype=np.int64)
permoff = np.zeros(10, dtype=np.int64); chunks = []; off = 0
for k in range(10):
    permoff[k] = off
    P_ = np.array(list(itertools.permutations(range(k))), dtype=np.int64).reshape(-1) if k > 0 else np.zeros(0, dtype=np.int64)
    chunks.append(P_); off += len(P_)
permtab = np.concatenate(chunks)


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
            if s[i] == s[i - 1]:
                run += 1
            else:
                run = 1
            if run > best:
                best = run
        out[t] = best
    return out


@nb.njit(parallel=True, cache=True)
def ears(n, parents):
    p = parents.shape[0]
    full = (1 << (n - 1)) - 1
    Hs = np.zeros((p, 1 << (n - 1)), dtype=np.int64)
    for t in nb.prange(p):
        adj = np.zeros(n, dtype=np.int64)
        for S in range(1, full):
            for v in range(n - 1):
                a = parents[t, v]
                if not (S >> v) & 1:
                    a |= 1 << (n - 1)
                adj[v] = a
            adj[n - 1] = S
            Hs[t, S] = ham_one(n, adj)
    return Hs


# ---------------------------------------------------------------- exact Python canonical form (any n)
def canon_py(adj):
    """Minimum code over all relabelings that respect an iterated colour refinement and individualisation (exact:
    the search tree explores every individualisation of the first non-singleton cell)."""
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

    def code_of(order):
        pos = {v: i for i, v in enumerate(order)}
        c = 0; e = 0
        for i in range(n):
            for j in range(i + 1, n):
                if outn[order[i]][order[j]]:
                    c |= 1 << e
                e += 1
        return c

    def search(col):
        col = refine(col)
        if len(set(col)) == n:
            order = sorted(range(n), key=lambda v: col[v])
            c = code_of(order)
            if best[0] is None or c < best[0]:
                best[0] = c
            return
        cells = defaultdict(list)
        for v in range(n):
            cells[col[v]].append(v)
        target = min((len(vs), k) for k, vs in cells.items() if len(vs) > 1)[1]
        for v in cells[target]:
            new = [2 * c + (1 if c > target or (c == target and u != v) else 0) for u, c in enumerate(col)]
            search(new)

    search([bin(a).count('1') for a in adj])
    return best[0]


MINH = {}


def build_children(n, parents, lim, exact=None, py=False):
    idxp = pair_index(n) if not py else None
    res = {}
    B = 400
    for b0 in range(0, len(parents), B):
        par = np.array(parents[b0:b0 + B], dtype=np.int64)
        Hs = ears(n, par)
        mh = int(Hs[:, 1:-1].min())
        MINH[n] = min(MINH.get(n, mh), mh)
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
        if py:
            for row, h in zip(cand, hv):
                c = canon_py([int(x) for x in row])
                res.setdefault(c, (int(h), [int(x) for x in row]))
        else:
            assert max_block(n, cand).max() <= 9, 'invariant block larger than 9'
            codes, au = canon_refined(n, cand, idxp, permtab, permoff, fact)
            for c, h, a in zip(codes, hv, au):
                res.setdefault(int(c), (int(h), int(a)))
    return res


def finish(n, res, py=False):
    """Add |Aut| and self-converse flags. Returns dict code -> (H, aut, sc)."""
    out = {}
    if py:
        for c, (h, row) in res.items():
            fullm = (1 << n) - 1
            conv = [fullm ^ row[v] ^ (1 << v) for v in range(n)]
            cc = canon_py(conv)
            # |Aut| by counting relabelings achieving the canonical code is expensive in Python; the extremal ones
            # are verified rigid below by direct automorphism search
            out[c] = (h, None, cc == c)
        return out
    codes = np.array(list(res.keys()), dtype=np.int64)
    A = np.array([code_to_adj(n, int(c)) for c in codes], dtype=np.int64)
    fullm = (1 << n) - 1
    conv = fullm ^ A ^ (np.int64(1) << np.arange(n, dtype=np.int64))[None, :]
    assert max_block(n, conv).max() <= 9
    cc, _ = canon_refined(n, conv, pair_index(n), permtab, permoff, fact)
    for c, x in zip(codes, cc):
        h, a = res[int(c)]
        out[int(c)] = (h, a, int(x) == int(c))
    return out


def canon_n(adj):
    n = len(adj)
    if n <= 11:
        A = np.array([adj], dtype=np.int64)
        assert max_block(n, A).max() <= 9
        codes, au = canon_refined(n, A, pair_index(n), permtab, permoff, fact)
        return int(codes[0]), int(au[0])
    return canon_py(adj), None


def aut_count(adj):
    """Exact |Aut| by backtracking (small n)."""
    n = len(adj)
    outn = [[(adj[u] >> v) & 1 for v in range(n)] for u in range(n)]
    sc = [bin(a).count('1') for a in adj]
    cnt = 0
    img = [-1] * n; used = [False] * n

    def bt(i):
        nonlocal cnt
        if i == n:
            cnt += 1
            return
        for w in range(n):
            if not used[w] and sc[w] == sc[i] and all(outn[i][j] == outn[w][img[j]] for j in range(i)):
                img[i] = w; used[w] = True
                bt(i + 1)
                used[w] = False
        img[i] = -1

    bt(0)
    return cnt


def c3count(adj):
    n = len(adj)
    return sum(1 for a, b, d in itertools.combinations(range(n), 3)
               if ((adj[a] >> b) & 1 and (adj[b] >> d) & 1 and (adj[d] >> a) & 1)
               or ((adj[a] >> d) & 1 and (adj[d] >> b) & 1 and (adj[b] >> a) & 1))


# ---------------------------------------------------------------- (E0)
print(f'(E0) interval chains [LIMIT = {LIMIT}]')
tested = 0; bad = 0
for n in range(3, 15):
    for c in comps(n - 1, (2, 3, 4, 5, 6)):
        tested += 1
        if Hof(chain(c)) != prodH(c):
            bad += 1
check(tested > 0 and bad == 0, f'H(IC(p_1..p_k)) = prod (2^(p_i - 1) + 1) for all {tested} compositions of n - 1 <= 13 into parts 2..6 (PROVED in the note; checked)')
check(all(prodH(c) >= busch(1 + sum(c)) for n in range(3, 18) for c in comps(n - 1, (2, 3, 4, 5, 6))),
      'no interval chain beats Busch (consistency); the optimal ones use parts 2, 3 and at most one 4 (or two 2s)')

# ---------------------------------------------------------------- census n <= 9
data = HM.classes_upto(9)
strong = {}
for m in range(3, 10):
    A = data[m]['A']
    H, st = HM.ham_all(m, A)
    strong[m] = [(A[i].copy(), int(H[i]), int(data[m]['aut'][i]), bool(data[m]['sc'][i])) for i in range(len(H)) if st[i]]
print(f'census n <= 9 loaded [{time.time() - T0:.0f}s]; strong 9-classes with H <= {LIMIT}: {sum(1 for x in strong[9] if x[1] <= LIMIT)}', flush=True)

# ---------------------------------------------------------------- (E1)
print(f'(E1) strong classes of order 10{" and 11" if TOP == 12 else ""} with H <= {LIMIT}')
levels = {}
par = [x[0] for x in strong[9] if x[1] <= LIMIT]
for n in range(10, TOP):
    res = build_children(n, par, LIMIT)
    levels[n] = finish(n, res)
    Sn = Counter(h for h, a, s in levels[n].values())
    print(f'     order {n}: {len(levels[n])} classes; least H = {min(Sn)}; first values {sorted(Sn)[:10]} [{time.time() - T0:.0f}s]', flush=True)
    check(min(Sn) == busch(n), f'the least H of a strong {n}-tournament is f({n}) = {busch(n)} (Busch, CITED; exact here)')
    check(all(h % a == 0 for h, a, s in levels[n].values()), f'order {n}: |Aut| divides H for every class found')
    check(all((Sn[h] - sum(1 for hh, a, s in levels[n].values() if hh == h and s)) % 2 == 0 for h in Sn),
          f'order {n}: converse pairing P_d[{n}] = SC_d[{n}] (mod 2)')
    if n % 2 == 0:
        check(all(h % 4 == 1 for h, a, s in levels[n].values() if s), f'order {n} (even): every self-converse class found has H = 1 (mod 4)')
    par = [code_to_adj(n, c) for c in levels[n]]

# ---------------------------------------------------------------- (E2)
print(f'(E2) strong classes of order {TOP} with H = f({TOP}) = {busch(TOP)}')
resT = build_children(TOP, par, LIMIT, exact=busch(TOP), py=(TOP >= 12))
top = finish(TOP, resT, py=(TOP >= 12))
lowest = MINH[TOP]
check(lowest == busch(TOP), f'no strong {TOP}-tournament has H < {busch(TOP)} (least ear value {lowest})')
print(f'     {len(top)} classes [{time.time() - T0:.0f}s]', flush=True)

# ---------------------------------------------------------------- (E3)
print(f'(E3) strong polynomials P_d for every odd d < {LIMIT}')
P = defaultdict(dict)
for m in range(3, 10):
    for adj, h, a, s in strong[m]:
        if h < LIMIT:
            P[h][m] = P[h].get(m, 0) + 1
for n, lev in levels.items():
    for h, a, s in lev.values():
        if h < LIMIT:
            P[h][n] = P[h].get(n, 0) + 1
zeros = [d for d in range(1, LIMIT, 2) if not P.get(d)]
check(zeros == [1, 7, 21], f'odd d < {LIMIT} with P_d = 0: {zeros}')
for d in range(125, LIMIT, 2):
    print(f'     P_{d} = ' + ' + '.join(f'{c}z^{m}' for m, c in sorted(P[d].items())))
print(f'     P_d(1), odd d = 3..{LIMIT - 2}:', [sum(P[d].values()) for d in range(3, LIMIT, 2)])

# ---------------------------------------------------------------- (E4)
print('(E4) Busch-extremal strong classes against the optimal interval chains')
ext = {}
for m in range(3, 10):
    ext[m] = [(adj, a, s) for adj, h, a, s in strong[m] if h == busch(m)]
for n, lev in levels.items():
    ext[n] = [(code_to_adj(n, c), a, s) for c, (h, a, s) in lev.items() if h == busch(n)]
ext[TOP] = [(None, None, s) for c, (h, a, s) in top.items()]
E = {m: len(v) for m, v in ext.items()}
opt = {m: optimal(m) for m in range(3, TOP + 1)}
print('     (m, f(m), E(m), #optimal compositions):', [(m, busch(m), E[m], len(opt[m])) for m in range(3, TOP + 1)])
check(all(E[m] == len(opt[m]) for m in range(3, TOP + 1)), f'E(m) = number of optimal compositions for m = 3..{TOP}')
# identification
ok_id = True; ok_rigid = True; ok_pal = True; ok_sc = True; ok_c3 = True
for m in range(3, TOP + 1):
    ext_codes = set()
    if m < TOP:
        for adj, a, s in ext[m]:
            ext_codes.add(canon_n([int(x) for x in adj])[0])
    else:
        ext_codes = set(top.keys())
    chain_codes = {}
    for c in opt[m]:
        adj = chain(c)
        code = canon_n(adj)[0] if m <= 11 else canon_py(adj)
        chain_codes[c] = code
        ok_rigid &= (m < 4) or aut_count(adj) == 1
        conv = [((1 << m) - 1) ^ adj[v] ^ (1 << v) for v in range(m)]
        is_sc = (canon_n(conv)[0] if m <= 11 else canon_py(conv)) == code
        ok_pal &= (is_sc == (tuple(reversed(c)) == c))
        sc_ = sorted(bin(x).count('1') for x in adj)
        ok_sc &= sc_ == ([1, 1] + list(range(2, m - 2)) + [m - 2, m - 2] if m >= 4 else [1, 1, 1])
        ok_c3 &= c3count(adj) == m - 2
    ok_id &= (set(chain_codes.values()) == ext_codes) and len(set(chain_codes.values())) == len(opt[m])
check(ok_id, f'for every m = 3..{TOP} the Busch-extremal classes are exactly the optimal interval chains, and distinct optimal compositions give distinct classes')
check(ok_pal, 'an optimal interval chain is self-converse iff its composition is a palindrome')
check(ok_rigid, 'every optimal interval chain of order >= 4 has trivial automorphism group')
check(ok_sc and ok_c3, 'every optimal interval chain of order m >= 4 has scores (1,1,2,3,...,m-3,m-2,m-2) and exactly m - 2 cyclic triples')
nsc = {m: sum(1 for c in opt[m] if tuple(reversed(c)) == c) for m in range(3, TOP + 1)}
print('     self-converse extremal classes (palindromes), m = 3..%d:' % TOP, [nsc[m] for m in range(3, TOP + 1)])

# ---------------------------------------------------------------- (E5)
print(f'(E5) exact multiplicities N(n,h) for every odd h < {LIMIT}')
strongc = defaultdict(lambda: defaultdict(int))
for d, row in P.items():
    for m, c in row.items():
        strongc[m][d] += c
Aseq = HM.seq_counts(strongc, 30, hmax=LIMIT - 1)
for n in (10, 11, 12):
    pr = HM.predict(Aseq, n)
    print(f'     n = {n}: N(n,h), h in (63, 75, 105, 125, 135, 147, 189, 195, 223): ' +
          str({h: pr.get(h, 0) for h in (63, 75, 105, 125, 135, 147, 189, 195, 223)}))
check(all(HM.predict(Aseq, n).get(h, 0) == 0 for n in range(1, 31) for h in (7, 21)), 'N(n,7) = N(n,21) = 0 for every n <= 30 (consistency)')
ok = True
for n in range(3, 10):
    pr = HM.predict(Aseq, n)
    cen = Counter(int(h) for h in HM.ham_all(n, data[n]['A'])[0])
    ok &= all(pr.get(h, 0) == cen.get(h, 0) for h in range(1, LIMIT, 2))
check(ok, f'the class series reproduces the census N(n,h) for every n <= 9 and odd h < {LIMIT}')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
