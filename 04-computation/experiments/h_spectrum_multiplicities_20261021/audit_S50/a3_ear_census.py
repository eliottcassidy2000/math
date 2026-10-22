"""Audit S50 - a3: independent inventory of strong tournaments with few Hamiltonian paths, from C_3 upward.

Method (independent of the session's census and code):
  * completeness: Moon (every strong m-tournament, m >= 4, has v with T - v strong) + insertion lemma (H(T-v) <= H(T)),
    so every strong m-class with H <= X is a nonconstant one-vertex extension of a strong (m-1)-class with H <= X;
  * the H value of every extension of a parent is read off one quadratic form (alib.ear_values, self-tested in a0);
  * exact canonical forms by my own refinement + individualisation code (alib.canon, self-tested against networkx);
  * start: the single strong 3-class C_3.  Orders 4..11: all strong classes with H <= 375.  Order 12: all extensions are
    evaluated; the least value is recorded and the classes with H = 375 are kept.
Sanity: with no H cap, orders 4..8 must give the known strong-class counts 1, 6, 35, 353, 6008 (OEIS A051337).
Writes inventory_375.json (codes, H, self-converse flag, |Aut| for orders 3..11 and the order-12 extremal classes)."""
import sys, os, time, json
from multiprocessing import Pool
import numpy as np
import alib as L

LIMIT = 375


def work(args):
    n1, code, cap, exact = args
    adj = L.decode(n1, code)
    vals = L.ear_values(n1, np.array(adj, dtype=np.int64), np.array(L.converse(adj), dtype=np.int64))
    full = (1 << n1) - 1
    v = vals[1:full]
    mn = int(v.min())
    out = {}
    if exact is not None:
        idx = np.nonzero(v == exact)[0] + 1
    elif cap is None:
        idx = np.arange(1, full)
    else:
        idx = np.nonzero(v <= cap)[0] + 1
    for x in idx:
        x = int(x)
        child = L.extend(adj, x)
        c = L.canon(child)
        if c not in out:
            out[c] = int(vals[x])
    return mn, out


def grow(pool, n, parents, cap=None, exact=None):
    tasks = [(n - 1, c, cap, exact) for c in parents]
    res = {}; mn = None
    for m, out in pool.imap_unordered(work, tasks, chunksize=8):
        mn = m if mn is None else min(mn, m)
        for c, h in out.items():
            if c in res:
                assert res[c] == h
            else:
                res[c] = h
    return res, mn


def finish(n, codes):
    """self-converse flag and |Aut| for each class."""
    out = {}
    for c in codes:
        adj = L.decode(n, c)
        cc, aut = L.canon(adj, want_aut=True)
        assert cc == c
        sc = L.canon(L.converse(adj)) == c
        out[c] = (sc, aut)
    return out


def fin_work(args):
    n, c = args
    return c, finish(n, [c])[c]


if __name__ == '__main__':
    T0 = time.time()
    FAIL = []

    def check(c, m):
        print(('  ok   ' if c else '  FAIL ') + m, flush=True)
        if not c:
            FAIL.append(m)

    c3 = L.canon(L.chain((2,)))
    with Pool(12) as pool:
        # sanity: uncapped strong census to order 8
        lev = {3: {c3: 3}}
        for n in range(4, 9):
            lev[n], _ = grow(pool, n, list(lev[n - 1]), cap=None)
        cnts = [len(lev[n]) for n in range(3, 9)]
        check(cnts == [1, 1, 6, 35, 353, 6008], f'uncapped ear census reproduces the strong-class counts 1, 1, 6, 35, 353, 6008 (A051337) for orders 3..8: {cnts} [{time.time() - T0:.0f}s]')
        full8 = lev[8]
        # capped inventory
        inv = {3: {c3: 3}}
        mins = {}
        for n in range(4, 12):
            par = [c for c, h in inv[n - 1].items() if h <= LIMIT]
            inv[n], mins[n] = grow(pool, n, par, cap=LIMIT)
            print(f'     order {n}: {len(inv[n])} strong classes with H <= {LIMIT}; least H {min(inv[n].values())}; '
                  f'least ear value {mins[n]} [{time.time() - T0:.0f}s]', flush=True)
        check(all(inv[n] == {c: h for c, h in lev[n].items() if h <= LIMIT} for n in range(4, 9)),
              'capped inventory = uncapped census restricted to H <= 375 at orders 4..8')
        par = list(inv[11])
        top, mins[12] = grow(pool, 12, par, exact=LIMIT)
        print(f'     order 12: least ear value {mins[12]}; classes with H = {LIMIT}: {len(top)} [{time.time() - T0:.0f}s]', flush=True)
        # symmetry data
        fin = {}
        for n in list(range(3, 12)) + [12]:
            src = inv[n] if n < 12 else top
            fin[n] = dict(pool.map(fin_work, [(n, c) for c in src], chunksize=16))
        print(f'     symmetry data done [{time.time() - T0:.0f}s]', flush=True)

    data = {'limit': LIMIT, 'orders': {}}
    for n in list(range(3, 12)) + [12]:
        src = inv[n] if n < 12 else top
        data['orders'][str(n)] = [[str(c), h, fin[n][c][0], fin[n][c][1]] for c, h in sorted(src.items(), key=lambda t: (t[1], t[0]))]
    data['uncapped8'] = [[str(c), h] for c, h in full8.items()]
    data['least_ear_value'] = {str(k): int(v) for k, v in mins.items()}
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inventory_375.json'), 'w') as fh:
        json.dump(data, fh)

    # ---------------------------------------------------------------- reports
    cnt = {n: sum(1 for h in inv[n].values() if h <= LIMIT) for n in range(3, 12)}
    print('     strong classes with H <= 375 by order 3..11:', [cnt[n] for n in range(3, 12)])
    n9_225 = sum(1 for h in inv[9].values() if h <= 225); n10_225 = sum(1 for h in inv[10].values() if h <= 225)
    check(cnt[9] == 3492 and n9_225 == 816, f'order 9: {cnt[9]} classes with H <= 375 (session 3492), {n9_225} with H <= 225 (session 816)')
    check(cnt[10] == 1739 and n10_225 == 239, f'order 10: {cnt[10]} classes with H <= 375 (session 1739), {n10_225} with H <= 225 (session 239)')
    check(cnt[11] == 411, f'order 11: {cnt[11]} classes with H <= 375 (session 411)')
    for n in range(4, 13):
        check(mins[n] == L.busch(n) if n >= 4 else True, f'order {n}: least H over all strong {n}-classes (ear minimum) = {mins[n]} = f({n}) = {L.busch(n)}')
    # extremal classes
    E = {}; SCE = {}; ident = True; rigid = True
    for m in range(3, 13):
        src = inv[m] if m < 12 else top
        ext = {c for c, h in src.items() if h == L.busch(m)}
        opt = [p for p in L.compositions(m - 1, 2, 4) if L.prod_law(p) == L.busch(m)]
        chain_codes = {L.canon(L.chain(p)) for p in opt}
        E[m] = len(ext)
        ident &= ext == chain_codes and len(chain_codes) == len(opt)
        SCE[m] = sum(1 for c in ext if fin[m][c][0])
        rigid &= all(fin[m][c][1] == (3 if m == 3 else 1) for c in ext)
    check(ident, f'for m = 3..12 the strong classes with H = f(m) are exactly the optimal chains: E(m) = {[E[m] for m in range(3, 13)]}')
    check([SCE[m] for m in range(3, 13)] == [1, 1, 2, 0, 1, 1, 1, 1, 3, 0] and rigid,
          f'self-converse extremal classes {[SCE[m] for m in range(3, 13)]}; all extremal classes of order >= 4 rigid')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
