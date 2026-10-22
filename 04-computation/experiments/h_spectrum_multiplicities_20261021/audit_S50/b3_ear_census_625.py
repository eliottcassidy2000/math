"""Audit S50, second pass - b3: independent inventory of strong tournaments with H <= 625 (own code; same method as
a3_ear_census.py: Moon + insertion lemma completeness, the quadratic-form extension values of alib.ear_values,
own refinement/individualisation canonical forms), started from C_3.
Checks THM-4744 (3), (4) as retyped: 13,439 / 8,769 / 3,536 / 658 strong classes with H <= 625 at orders 9-12;
no strong 13-tournament has H < 625; E(13) = 1, the chain IC(3,3,3,3); odd d < 625 missing from every order: 1, 7, 21;
and the session's printed list P_d(1), odd d = 375..623 (strong_ears_625.out, read as data)."""
import os, re, time
from multiprocessing import Pool
from collections import Counter
import numpy as np
import alib as L

LIMIT = 625
HERE = os.path.dirname(os.path.abspath(__file__))


def work(args):
    n1, code, cap, exact = args
    adj = L.decode(n1, code)
    vals = L.ear_values(n1, np.array(adj, dtype=np.int64), np.array(L.converse(adj), dtype=np.int64))
    full = (1 << n1) - 1
    v = vals[1:full]
    mn = int(v.min())
    out = {}
    idx = (np.nonzero(v == exact)[0] + 1) if exact is not None else (np.nonzero(v <= cap)[0] + 1)
    for x in idx:
        x = int(x)
        c = L.canon(L.extend(adj, x))
        if c not in out:
            out[c] = int(vals[x])
    return mn, out


def grow(pool, n, parents, cap=None, exact=None):
    res = {}; mn = None
    for m, out in pool.imap_unordered(work, [(n - 1, c, cap, exact) for c in parents], chunksize=4):
        mn = m if mn is None else min(mn, m)
        for c, h in out.items():
            assert res.setdefault(c, h) == h
    return res, mn


if __name__ == '__main__':
    T0 = time.time()
    FAIL = []

    def check(c, m):
        print(('  ok   ' if c else '  FAIL ') + m, flush=True)
        if not c:
            FAIL.append(m)

    inv = {3: {L.canon(L.chain((2,))): 3}}
    mins = {}
    with Pool(12) as pool:
        for n in range(4, 13):
            par = [c for c, h in inv[n - 1].items() if h <= LIMIT]
            inv[n], mins[n] = grow(pool, n, par, cap=LIMIT)
            print(f'     order {n}: {len(inv[n])} strong classes with H <= {LIMIT}; least ear value {mins[n]} [{time.time() - T0:.0f}s]', flush=True)
        top, mins[13] = grow(pool, 13, list(inv[12]), exact=LIMIT)
    print(f'     order 13: least ear value {mins[13]}; classes with H = 625: {len(top)} [{time.time() - T0:.0f}s]')
    cnt = [len(inv[n]) for n in (9, 10, 11, 12)]
    check(cnt == [13439, 8769, 3536, 658], f'strong classes with H <= 625 at orders 9..12: {cnt} (THM-4744 (3): 13,439 / 8,769 / 3,536 / 658)')
    check(all(mins[n] == L.busch(n) for n in range(4, 14)), f'least ear value = f(n) for n = 4..13: {[mins[n] for n in range(4, 14)]}')
    check(set(top) == {L.canon(L.chain((3, 3, 3, 3)))}, 'E(13) = 1: the only strong 13-class with H = 625 is IC(3,3,3,3)')
    P = Counter()
    for n, d in inv.items():
        for c, h in d.items():
            if h < LIMIT:
                P[h] += 1
    zeros = [d for d in range(1, LIMIT, 2) if P[d] == 0]
    check(zeros == [1, 7, 21], f'odd d < 625 attained by no strong tournament: {zeros}')
    mine = [P[d] for d in range(375, LIMIT, 2)]
    txt = open(os.path.join(os.path.dirname(HERE), 'strong_ears_625.out')).read()
    m = re.search(r'P_d\(1\) for odd d = 375\.\.623: \[([^\]]*)\]', txt)
    sess = [int(x) for x in m.group(1).split(',')] if m else None
    check(sess == mine, f'P_d(1) for odd d = 375..623 agrees with strong_ears_625.out entry by entry ({len(mine)} values)')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
