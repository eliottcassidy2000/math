"""S51 audit c4: the auditor's own complete inventory of strong tournament classes with H <= 2100 (all orders), by exact
ears from C_3 with the auditor's own ear values (alib51.ear_all) and canonical form (alib51.canon); nothing from the
session or from the S50 audit is imported.  Completeness: every strong n-tournament (n >= 4) has a vertex v with T - v
strong (Moon 1966), and H(T - v) <= H(T) (insert v into a Hamiltonian path of T - v), so every strong class with
H <= 2100 is a non-constant ear of a strong class with H <= 2100.
Checks: class counts by order against repunit_rarity.out (I); the THM-4744 counts below 625; every H re-verified by an
exact DP on the decoded representative; strongness; |Aut| divides H; every odd 3 <= d < 2100 except 7, 21 is a strong
count; the S(d) values printed by the session; and (if it is plain data) the session's own pickled inventory
inv_2100.pkl, compared value-multiset by value-multiset per order (read with a restricted unpickler: only built-in
dict/int/tuple opcodes are accepted).
Writes inventory_2100.npz to the scratchpad (cache for c5, c6).  Prints ALL CHECKS PASSED or the failures."""
import os, sys, time, math, pickle, pickletools, io
from collections import Counter, defaultdict
import numpy as np
import numba as nb
import alib51 as L

FAIL = []
T0 = time.time()
LIMIT = 2100
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'
SESS = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\s51\inv_2100.pkl'
os.makedirs(SCR, exist_ok=True)


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(parallel=True, cache=True)
def children(n1, parents, cap):
    """for each parent (order n1), every non-constant ear x with H(child) <= cap: rows (w1, w0, H) in a per-parent slot"""
    P = parents.shape[0]
    X = 1 << n1
    out = np.full((P * X, 3), -1, dtype=np.int64)
    mins = np.zeros(P, dtype=np.int64)
    for t in nb.prange(P):
        ev = L.ear_all(n1, parents[t])
        mn = 1 << 62
        k = 0
        for x in range(1, X - 1):
            h = ev[x]
            if h < mn:
                mn = h
            if h <= cap:
                ch = L.extend(n1, parents[t], x)
                w1, w0, a = L.canon(n1 + 1, ch)
                out[t * X + k, 0] = w1
                out[t * X + k, 1] = w0
                out[t * X + k, 2] = h
                k += 1
        mins[t] = mn
    return out, mins


@nb.njit(parallel=True, cache=True)
def verify(n, codes):
    K = codes.shape[0]
    H = np.zeros(K, dtype=np.int64); st = np.zeros(K, dtype=np.int64); aut = np.zeros(K, dtype=np.int64)
    A = np.zeros((K, n), dtype=np.int64)
    full = (1 << n) - 1
    for t in nb.prange(K):
        adj = L.decode(n, codes[t, 0], codes[t, 1])
        for v in range(n):
            A[t, v] = adj[v]
        H[t] = L.ham(n, adj)
        st[t] = 1 if L.is_strong_mask(n, adj, full) else 0
        w1, w0, a = L.canon(n, adj)
        aut[t] = a
        if w1 != codes[t, 0] or w0 != codes[t, 1]:
            aut[t] = -1
    return H, st, aut, A


def build_inventory():
    c3 = L.canon(3, L.reverse_arcs(3, [(2, 0)]))
    inv = {3: dict(codes=np.array([[c3[0], c3[1]]], dtype=np.int64), H=np.array([3], dtype=np.int64))}
    mins = {}
    n = 4
    while True:
        par = np.array([L.decode(n - 1, int(w1), int(w0)) for w1, w0 in inv[n - 1]['codes']], dtype=np.int64)
        X = 1 << (n - 1)
        chunk = max(16, (1 << 22) // X)
        rows = []
        mn = None
        for s in range(0, len(par), chunk):
            out, mm = children(n - 1, par[s:s + chunk], LIMIT)
            out = out[out[:, 0] >= 0]
            rows.append(np.unique(out, axis=0))
            m0 = int(mm.min())
            mn = m0 if mn is None else min(mn, m0)
        allr = np.unique(np.concatenate(rows), axis=0) if rows else np.zeros((0, 3), dtype=np.int64)
        codes, idx = np.unique(allr[:, :2], axis=0, return_index=True)
        assert len(codes) == len(allr), 'one class received two different H values'
        mins[n] = mn
        print(f'     order {n}: {len(codes)} strong classes with H <= {LIMIT}; least ear value {mn} [{time.time() - T0:.0f}s]', flush=True)
        if len(codes) == 0:
            break
        inv[n] = dict(codes=codes, H=allr[idx, 2])
        n += 1
    return inv, mins


class SafeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f'refused global {module}.{name}')


if __name__ == '__main__':
    print('(I) strong inventory below', LIMIT, 'by the auditor\'s own ear census')
    inv, mins = build_inventory()
    counts = {n: len(inv[n]['H']) for n in sorted(inv)}
    sess = {3: 1, 4: 1, 5: 6, 6: 35, 7: 353, 8: 6008, 9: 146376, 10: 268683, 11: 180284, 12: 88379, 13: 28817, 14: 4545, 15: 56}
    check(counts == sess, f'classes by order {counts} = repunit_rarity.out (I); total {sum(counts.values())} (session: 723,544)')
    busch = lambda m: {0: 3 * 5 ** (m // 3 - 1), 1: 5 ** (m // 3), 2: 9 * 5 ** (m // 3 - 1)}[m % 3]
    check(all(mins[n] == busch(n) for n in mins), f'least ear value at every order = Busch f(n): {mins}')
    c625 = {m: int((inv[m]['H'] <= 625).sum()) for m in (9, 10, 11, 12, 13)}
    check(c625 == {9: 13439, 10: 8769, 11: 3536, 12: 658, 13: 1}, f'classes with H <= 625 at orders 9-13: {c625} (THM-4744 / S50 audit b3)')
    print('(V) verification of every class: exact DP, strongness, |Aut|')
    ok = True; okst = True; okaut = True
    for n in sorted(inv):
        H, st, aut, A = verify(n, inv[n]['codes'])
        ok &= bool(np.all(H == inv[n]['H']))
        okst &= bool(np.all(st == 1))
        okaut &= bool(np.all(aut > 0)) and bool(np.all(H % aut == 0))
        inv[n]['aut'] = aut
        inv[n]['A'] = A
    check(ok, 'for every class the exact Hamiltonian-path DP equals the ear value')
    check(okst, 'every class is strong')
    check(okaut, 'canonical codes are stable under re-canonicalisation and |Aut| divides H for every class')
    S = Counter(); Sm = defaultdict(Counter)
    for n in inv:
        for h in inv[n]['H']:
            S[int(h)] += 1; Sm[n][int(h)] += 1
    check(all(S[d] > 0 for d in range(3, LIMIT, 2) if d not in (7, 21)) and S[7] == 0 and S[21] == 0 and S[1] == 0,
          'every odd 3 <= d < 2100 except 7, 21 is a strong count (and 1, 7, 21 are not)')
    want = {31: 4, 33: 8, 63: 2, 65: 14, 127: 12, 129: 40, 255: 107, 257: 73, 511: 164, 513: 265, 1023: 663, 1025: 537, 2047: 1599, 2049: 1968, 1365: 1037}
    got = {d: S[d] for d in want}
    check(got == want, f'S(d) at the session\'s printed values: {got}')
    np.savez_compressed(os.path.join(SCR, 'inventory_2100.npz'), **{f'{k}{n}': v for n in inv for k, v in inv[n].items()})
    print(f'     cached to {os.path.join(SCR, "inventory_2100.npz")} [{time.time() - T0:.0f}s]')
    print('(P) the session\'s own pickled inventory (data only)')
    if os.path.exists(SESS):
        raw = open(SESS, 'rb').read()
        ops = Counter(op.name for op, arg, pos in pickletools.genops(raw))
        bad = [o for o in ops if 'GLOBAL' in o or o in ('REDUCE', 'BUILD', 'INST', 'OBJ', 'NEWOBJ', 'NEWOBJ_EX', 'EXT1', 'EXT2', 'EXT4', 'PERSID', 'BINPERSID')]
        check(not bad, f'inv_2100.pkl uses only plain-data opcodes ({sorted(ops)[:12]}...)')
        if not bad:
            sinv = SafeUnpickler(io.BytesIO(raw)).load()
            sc = {int(n): Counter(int(h) for h in d.values()) for n, d in sinv.items()}
            agree = all(sc.get(n, Counter()) == Sm[n] for n in set(sc) | set(Sm))
            check(agree, 'the session\'s inv_2100.pkl has exactly the auditor\'s multiset of H values at every order (class by class counts per value)')
            print(f'     session pickle mtime: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(SESS)))} local')
    else:
        print('     (session pickle not found; skipped)')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
