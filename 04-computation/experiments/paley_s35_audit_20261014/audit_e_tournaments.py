#!/usr/bin/env python3
"""Audit (e) of S35 / HYP-9277: balance of uniformly random Hamiltonian paths in tournaments.

delta(T) = max over pairs x != y of min(p, 1 - p), p = Pr[x precedes y] over the Hamiltonian paths of T.
Methods (independent of the session's forward/backward subset DP):
  n <= 6: brute force over all n! permutations for every labelled tournament;
  n = 7 : exhaustive iterative-DFS enumeration of the Hamiltonian paths of all 2^21 labelled tournaments (numba,
          parallel), cross-checked against permutation brute force on random tournaments;
  n = 8, 9, 10: the stated counterexamples by DFS enumeration (n = 8 also by 8! permutations).
Labelled tournaments: bit e of the code (pairs (i,j), i<j, lexicographic) set <=> i -> j  (any bijective encoding
gives the same labelled statistics).
"""
import itertools, math, time
from fractions import Fraction
from collections import Counter
import numpy as np
import numba as nb


def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def adj_from_code(n, code):
    adj = [0] * n
    for e, (i, j) in enumerate(pairs_of(n)):
        if code >> e & 1:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return adj


@nb.njit(cache=True)
def stats_dfs(n, adj, C):
    """enumerate Hamiltonian paths by DFS; C[x, y] += 1 when x precedes y.  Returns H."""
    for x in range(n):
        for y in range(n):
            C[x, y] = 0
    H = 0
    path = np.zeros(n, dtype=np.int64)
    nxt = np.zeros(n + 1, dtype=np.int64)
    for s in range(n):
        path[0] = s
        used = 1 << s
        d = 1
        nxt[1] = 0
        while d >= 1:
            if d == n:
                H += 1
                for i in range(n):
                    for j in range(i + 1, n):
                        C[path[i], path[j]] += 1
                d -= 1
                used &= ~(1 << path[d]) if d >= 1 else used
                continue
            found = False
            while nxt[d] < n:
                w = nxt[d]
                nxt[d] += 1
                if not (used >> w) & 1 and (adj[path[d - 1]] >> w) & 1:
                    path[d] = w
                    used |= 1 << w
                    d += 1
                    nxt[d] = 0
                    found = True
                    break
            if not found:
                d -= 1
                if d >= 1:
                    used &= ~(1 << path[d])
    return H


@nb.njit(cache=True)
def stats_perm(n, adj, perms, C):
    for x in range(n):
        for y in range(n):
            C[x, y] = 0
    H = 0
    for r in range(perms.shape[0]):
        ok = True
        for i in range(n - 1):
            if not (adj[perms[r, i]] >> perms[r, i + 1]) & 1:
                ok = False
                break
        if ok:
            H += 1
            for i in range(n):
                for j in range(i + 1, n):
                    C[perms[r, i], perms[r, j]] += 1
    return H


@nb.njit(cache=True)
def strong(n, adj):
    full = (1 << n) - 1
    seen = 1
    fr = 1
    while fr:
        nx_ = 0
        for v in range(n):
            if (fr >> v) & 1:
                nx_ |= adj[v]
        nx_ &= ~seen
        seen |= nx_
        fr = nx_
    if seen != full:
        return False
    seen = 1
    fr = 1
    while fr:
        nx_ = 0
        for v in range(n):
            if (fr >> v) & 1:
                for u in range(n):
                    if (adj[u] >> v) & 1:
                        nx_ |= 1 << u
        nx_ &= ~seen
        seen |= nx_
        fr = nx_
    return seen == full


@nb.njit(cache=True)
def comp_sizes_ok(n, adj):
    """True iff every strong component has size 1 or 3 and at least one has size 3."""
    # reach[v] = set reachable from v
    reach = np.zeros(n, dtype=np.int64)
    for v in range(n):
        seen = 1 << v
        fr = 1 << v
        while fr:
            nx_ = 0
            for u in range(n):
                if (fr >> u) & 1:
                    nx_ |= adj[u]
            nx_ &= ~seen
            seen |= nx_
            fr = nx_
        reach[v] = seen
    has3 = False
    for v in range(n):
        comp = 0
        for u in range(n):
            if (reach[v] >> u) & 1 and (reach[u] >> v) & 1:
                comp |= 1 << u
        sz = 0
        for u in range(n):
            if (comp >> u) & 1:
                sz += 1
        if sz != 1 and sz != 3:
            return False
        if sz == 3:
            has3 = True
    return has3


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk, use_perm, perms):
    E = pu.shape[0]
    total = 1 << E
    chunk = (total + nchunk - 1) // nchunk
    # per chunk: count delta==1/3, delta<1/3 (nontransitive), mismatch with composition prediction,
    # min strong delta (num, den), min arc ratio (num, den), strong tournaments without balanced arc
    at13 = np.zeros(nchunk, dtype=np.int64)
    below = np.zeros(nchunk, dtype=np.int64)
    mism = np.zeros(nchunk, dtype=np.int64)
    sn = np.ones(nchunk, dtype=np.int64)
    sd = np.ones(nchunk, dtype=np.int64)
    an = np.ones(nchunk, dtype=np.int64)
    ad = np.ones(nchunk, dtype=np.int64)
    nobal = np.zeros(nchunk, dtype=np.int64)
    nstrong = np.zeros(nchunk, dtype=np.int64)
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        C = np.zeros((n, n), dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                adj[i] = 0
            for e in range(E):
                if (code >> e) & 1:
                    adj[pu[e]] |= 1 << pv[e]
                else:
                    adj[pv[e]] |= 1 << pu[e]
            if use_perm:
                H = stats_perm(n, adj, perms, C)
            else:
                H = stats_dfs(n, adj, C)
            m = 0
            for x in range(n):
                for y in range(x + 1, n):
                    v = C[x, y]
                    if H - v < v:
                        v = H - v
                    if v > m:
                        m = v
            # arcs
            for x in range(n):
                for y in range(n):
                    if x != y and (adj[x] >> y) & 1:
                        if C[x, y] * ad[c] < an[c] * H:
                            an[c] = C[x, y]
                            ad[c] = H
            if H > 1:
                if 3 * m == H:
                    at13[c] += 1
                elif 3 * m < H:
                    below[c] += 1
                pred = comp_sizes_ok(n, adj)
                if pred != (3 * m == H):
                    mism[c] += 1
            if strong(n, adj):
                nstrong[c] += 1
                if m * sd[c] < sn[c] * H:
                    sn[c] = m
                    sd[c] = H
                bal = False
                for x in range(n):
                    for y in range(n):
                        if x != y and (adj[x] >> y) & 1:
                            if 3 * C[x, y] >= H and 3 * C[x, y] <= 2 * H:
                                bal = True
                if not bal:
                    nobal[c] += 1
    return at13, below, mism, sn, sd, an, ad, nobal, nstrong


def fr_min(nums, dens):
    return min(Fraction(int(a), int(b)) for a, b in zip(nums, dens))


def tour_stats(n, adj):
    C = np.zeros((n, n), dtype=np.int64)
    H = int(stats_dfs(n, np.array(adj, dtype=np.int64), C))
    delta = max(Fraction(min(int(C[x, y]), H - int(C[x, y])), H) for x in range(n) for y in range(x + 1, n))
    arcmin = min((Fraction(int(C[x, y]), H), x, y) for x in range(n) for y in range(n) if x != y and adj[x] >> y & 1)
    return H, C, delta, arcmin


def predicted_labelled(n):
    tot = 0
    for parts in range(n + 1):
        pass

    def rec(rem, b):
        nonlocal tot
        if rem == 0:
            if b:
                tot += math.factorial(n) // 3 ** b
            return
        rec(rem - 1, b)
        if rem >= 3:
            rec(rem - 3, b + 1)
    rec(n, 0)
    return tot


def main():
    t0 = time.time()
    for n in range(3, 8):
        pr = pairs_of(n)
        pu = np.array([p[0] for p in pr], dtype=np.int64)
        pv = np.array([p[1] for p in pr], dtype=np.int64)
        use_perm = n <= 6
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64) if use_perm else np.zeros((1, n), dtype=np.int64)
        at13, below, mism, sn, sd, an, ad, nobal, nstrong = sweep(n, pu, pv, 256, use_perm, perms)
        print(f'n = {n} ({"permutation brute force" if use_perm else "DFS path enumeration"}): delta = 1/3 at {int(at13.sum())} labelled '
              f'(composition count {predicted_labelled(n)}); non-transitive with delta < 1/3: {int(below.sum())}; '
              f'delta = 1/3 <=> components in {{1, C3}} mismatches: {int(mism.sum())}; least strong delta: {fr_min(sn, sd) if nstrong.sum() else None}; '
              f'min arc ratio p(tail first): {fr_min(an, ad)}; strong tournaments: {int(nstrong.sum())}, without a balanced arc: {int(nobal.sum())}  '
              f'[{time.time() - t0:.0f} s]', flush=True)
    # cross-check DFS against permutations on random 7-tournaments
    rng = np.random.default_rng(7)
    perms7 = np.array(list(itertools.permutations(range(7))), dtype=np.int64)
    C1 = np.zeros((7, 7), dtype=np.int64)
    C2 = np.zeros((7, 7), dtype=np.int64)
    bad = 0
    for _ in range(3000):
        adj = np.array(adj_from_code(7, int(rng.integers(0, 1 << 21))), dtype=np.int64)
        h1 = stats_dfs(7, adj, C1)
        h2 = stats_perm(7, adj, perms7, C2)
        if h1 != h2 or not np.array_equal(C1, C2):
            bad += 1
    print(f'n = 7: DFS = permutation brute force on 3000 random tournaments: {bad == 0}')
    # strong delta values by isomorphism class at n = 4, 5, 6 (via labelled sweep, distinct values)
    for n in (4, 5, 6):
        vals = Counter()
        ex = {}
        for code in range(1 << len(pairs_of(n))):
            adj = adj_from_code(n, code)
            if not strong(n, np.array(adj, dtype=np.int64)):
                continue
            H, C, d, _ = tour_stats(n, adj)
            vals[d] += 1
            if d not in ex:
                ex[d] = (H, tuple(sorted(bin(a).count('1') for a in adj)))
        ds = sorted(vals)
        print(f'n = {n}: distinct strong delta values (sorted): {[str(x) for x in ds[:10]]}{" ..." if len(ds) > 10 else ""} ({len(ds)} values); '
              f'least: {ds[0]} with (H, scores) = {ex[ds[0]]}')
    # n = 8 counterexample to Claim A
    adj8 = [190, 176, 250, 242, 64, 144, 163, 16]
    ok = all(((adj8[i] >> j) & 1) + ((adj8[j] >> i) & 1) == 1 for i in range(8) for j in range(i + 1, 8))
    perms8 = np.array(list(itertools.permutations(range(8))), dtype=np.int64)
    C = np.zeros((8, 8), dtype=np.int64)
    H = int(stats_perm(8, np.array(adj8, dtype=np.int64), perms8, C))
    H2, C2, d8, am8 = tour_stats(8, adj8)
    print(f'n = 8 Claim-A counterexample: valid tournament {ok}; arc 0 -> 4 present: {bool(adj8[0] >> 4 & 1)}; '
          f'8! brute force H = {H}, #(0 before 4) = {int(C[0, 4])} -> {Fraction(int(C[0, 4]), H)}; DFS agrees: {H == H2 and np.array_equal(C, C2)}; '
          f'delta = {d8}; min arc ratio {am8[0]} at arc {am8[1]} -> {am8[2]}; strong: {strong(8, np.array(adj8, dtype=np.int64))}')
    claims = [(8, [72, 197, 1, 6, 111, 143, 44, 93], 'strong min 64/161'),
              (9, [460, 493, 328, 256, 367, 461, 264, 92, 128], 'Claim A 17/81, delta 4/9'),
              (9, [508, 457, 498, 500, 98, 386, 416, 272, 16], 'strong min 104/257'),
              (10, [734, 220, 208, 116, 128, 983, 144, 8, 735, 222], 'Claim A 25/153, delta 8/17'),
              (10, [340, 781, 944, 741, 970, 595, 646, 291, 104, 385], 'strong min 292/693')]
    for n, adj, label in claims:
        ok = all(((adj[i] >> j) & 1) + ((adj[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n))
        H, C, d, am = tour_stats(n, adj)
        print(f'n = {n} ({label}): valid {ok}; H = {H}; delta = {d}; min arc ratio = {am[0]}; strong = {strong(n, np.array(adj, dtype=np.int64))}')
    # Paley tournaments
    for q in (3, 7, 11):
        sq = {(x * x) % q for x in range(1, q)}
        adj = [sum(1 << y for y in range(q) if (y - x) % q in sq) for x in range(q)]
        H, C, d, am = tour_stats(q, adj)
        print(f'Paley P_{q}: H = {H}, #(0 before 1) = {int(C[0, 1])}, p = {Fraction(int(C[0, 1]), H)}')
    print(f'[{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
