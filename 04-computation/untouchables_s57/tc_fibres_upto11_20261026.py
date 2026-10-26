#!/usr/bin/env python3
"""HYP-9352(2) for all tournaments on n <= 11 vertices: tc = H/|Aut| is never 7 or 21 (opus-2026-10-26-S57).

tc(T) = H(T)/|Aut(T)| is the number of marked-path presentations of T (THM-4741), i.e. the number of labelled graphs G
on [n-1] whose path tournament T_G (path 1 -> ... -> n, non-path pair {i, j+1} backward iff {i,j} in G) is isomorphic
to T.  Reduction (proved in the S57 note):
  (R1) tc is multiplicative over the strong components (THM-2450 sec. 2), and 7, 21 = 3*7 force a component with
       tc in {7, 21}; so a minimal counterexample is strong.
  (R2) H(S) >= f(n) for strong S (Busch; THM-4744): f = 1,-,3,5,9,15,25,45,75,125,225 at n = 1..11.  Hence
       |Aut(S)| >= f(n)/21, and for n >= 7 a rigid strong S has tc = H >= 25.
  (R3) |Aut| is odd; an odd-order group on n <= 11 points whose order is at least the R2 bound contains an element of
       one of the cycle types listed in TYPES[n] (group-theory lemma in the note: Cauchy for 5, 7, 11; a 3-group of
       order >= 9 contains an element with >= 2 three-cycles; at n = 11 an order >= 27 3-group contains (3)(3)(3) or a
       9-cycle, and an order 15*3^j group contains an element of order 15).
So it suffices to enumerate, for each listed cycle type, every tournament invariant under one fixed permutation of
that type (all isomorphism classes of such tournaments occur), keep the strong ones with H = 7a or 21a where a is a
multiple of the element order, and compute |Aut| exactly.  n <= 6 is done exhaustively over all labelled tournaments.
Self-tests: H and |Aut| of C3, T5, P7 and brute force on random small tournaments; Busch's f(n) exhaustively n <= 6.
Usage: python tc_fibres_upto11_20261026.py [nmax]
"""
import itertools, json, math, os, random, sys, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
BUSCH = {3: 3, 4: 5, 5: 9, 6: 15, 7: 25, 8: 45, 9: 75, 10: 125, 11: 225}
TYPES = {
    7: [(3,), (3, 3), (5,), (7,)],
    8: [(3,), (3, 3), (5,), (7,)],
    9: [(3, 3), (3, 3, 3), (5,), (7,)],
    10: [(3, 3), (3, 3, 3), (5,), (5, 5), (7,)],
    11: [(3, 3, 3), (5, 5), (5, 3), (5, 3, 3), (7,), (11,)],
}


def perm_from_type(n, cyc):
    p = list(range(n))
    start = 0
    for L in cyc:
        for i in range(L):
            p[start + i] = start + (i + 1) % L
        start += L
    return p


def pair_orbits(p, n):
    seen = set()
    orbs = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen:
                continue
            orb = []
            a, b = i, j
            while (min(a, b), max(a, b)) not in seen:
                seen.add((min(a, b), max(a, b)))
                orb.append((a, b))
                a, b = p[a], p[b]
            orbs.append(orb)
    return orbs


@nb.njit(cache=False)
def build(bits, n, oa, ob, ooff, adj):
    for i in range(n):
        for j in range(n):
            adj[i, j] = 0
    for t in range(ooff.shape[0] - 1):
        o = (bits >> t) & 1
        for e in range(ooff[t], ooff[t + 1]):
            if o:
                adj[oa[e], ob[e]] = 1
            else:
                adj[ob[e], oa[e]] = 1


@nb.njit(cache=False)
def is_strong(adj, n):
    for direction in range(2):
        seen = np.zeros(n, np.uint8)
        stack = np.zeros(n, np.int64)
        top = 1
        stack[0] = 0
        seen[0] = 1
        while top > 0:
            top -= 1
            v = stack[top]
            for w in range(n):
                e = adj[v, w] if direction == 0 else adj[w, v]
                if e and not seen[w]:
                    seen[w] = 1
                    stack[top] = w
                    top += 1
        for w in range(n):
            if not seen[w]:
                return False
    return True


@nb.njit(cache=False)
def ham(adj, n, dp):
    full = 1 << n
    for m in range(full):
        for v in range(n):
            dp[m, v] = 0
    for v in range(n):
        dp[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            c = dp[m, v]
            if c == 0:
                continue
            for w in range(n):
                if (m >> w) & 1 == 0 and adj[v, w]:
                    dp[m | (1 << w), w] += c
    s = 0
    for v in range(n):
        s += dp[full - 1, v]
    return s


@nb.njit(cache=False)
def aut_count(adj, n):
    """|Aut| by iterative backtracking with out-degree pruning."""
    score = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(n):
            score[i] += adj[i, j]
    order = np.argsort(score)
    img = -np.ones(n, np.int64)
    used = np.zeros(n, np.uint8)
    nxt = np.zeros(n + 1, np.int64)     # next candidate image to try at depth k
    k = 0
    cnt = 0
    while k >= 0:
        if k == n:
            cnt += 1
            k -= 1
            if k >= 0:
                used[img[order[k]]] = 0
                img[order[k]] = -1
            continue
        v = order[k]
        found = False
        w = nxt[k]
        while w < n:
            if used[w] == 0 and score[w] == score[v]:
                ok = True
                for t in range(k):
                    u = order[t]
                    if adj[u, v] != adj[img[u], w]:
                        ok = False
                        break
                if ok:
                    found = True
                    break
            w += 1
        if found:
            img[v] = w
            used[w] = 1
            nxt[k] = w + 1
            k += 1
            if k <= n - 1:
                nxt[k] = 0
        else:
            nxt[k] = 0
            k -= 1
            if k >= 0:
                used[img[order[k]]] = 0
                img[order[k]] = -1
    return cnt


@nb.njit(cache=False)
def scan_type(n, oa, ob, ooff, pmul, lo_aut):
    """All 2^m invariant tournaments: count strong ones; return those with H in {7a, 21a}, a = |Aut|, p | a."""
    m = ooff.shape[0] - 1
    adj = np.zeros((n, n), np.uint8)
    dp = np.zeros((1 << n, n), np.int64)
    nstrong = 0
    hmin = 1 << 62
    hits = []
    ncand = 0
    for bits in range(1 << m):
        build(bits, n, oa, ob, ooff, adj)
        if not is_strong(adj, n):
            continue
        nstrong += 1
        H = ham(adj, n, dp)
        if H < hmin:
            hmin = H
        cand = False
        for tc in (7, 21):
            if H % (tc * pmul) == 0 and H // tc >= lo_aut:
                cand = True
        if cand:
            ncand += 1
            A = aut_count(adj, n)
            if A > 0 and (H == 7 * A or H == 21 * A):
                hits.append((bits, H, A))
    return nstrong, ncand, hmin, hits


def exhaustive_small(n):
    """All labelled tournaments on n <= 6 vertices: tc spectrum, Busch minimum, strong-spectrum minimum."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    dp = np.zeros((1 << n, n), np.int64)
    tcs, stc, hmin = set(), set(), None
    for bits in range(1 << len(pairs)):
        adj = np.zeros((n, n), np.uint8)
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1:
                adj[i, j] = 1
            else:
                adj[j, i] = 1
        H = ham(adj, n, dp)
        A = aut_count(adj, n)
        assert H % A == 0
        tc = H // A
        tcs.add(tc)
        if is_strong(adj, n):
            stc.add(tc)
            hmin = H if hmin is None else min(hmin, H)
    return sorted(tcs), (min(stc) if stc else None), hmin


def self_tests():
    def circ(n, S):
        adj = np.zeros((n, n), np.uint8)
        for i in range(n):
            for s in S:
                adj[i, (i + s) % n] = 1
        return adj
    out = {}
    for name, n, S, H0, A0 in (('C3', 3, [1], 3, 3), ('T5=Z5{1,2}', 5, [1, 2], 15, 5), ('P7', 7, [1, 2, 4], 189, 21)):
        adj = circ(n, S)
        dp = np.zeros((1 << n, n), np.int64)
        H, A = int(ham(adj, n, dp)), int(aut_count(adj, n))
        assert (H, A) == (H0, A0), (name, H, A)
        out[name] = (H, A, H // A)
    rnd = random.Random(2026)
    for trial in range(60):
        n = rnd.choice([4, 5, 6, 7])
        adj = np.zeros((n, n), np.uint8)
        for i in range(n):
            for j in range(i + 1, n):
                if rnd.random() < 0.5:
                    adj[i, j] = 1
                else:
                    adj[j, i] = 1
        dp = np.zeros((1 << n, n), np.int64)
        H = int(ham(adj, n, dp))
        Hb = sum(all(adj[p[i], p[i + 1]] for i in range(n - 1)) for p in itertools.permutations(range(n)))
        Ab = sum(all(adj[i, j] == adj[p[i], p[j]] for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))
        assert H == Hb and int(aut_count(adj, n)) == Ab, (n, H, Hb)
    out['random_bruteforce_checks'] = 60
    return out


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    t0 = time.time()
    res = {'self_tests': self_tests()}
    print('self-tests passed:', res['self_tests'], flush=True)
    for n in range(3, 7):
        tcs, smin, hmin = exhaustive_small(n)
        bad = [t for t in tcs if t in (7, 21)]
        assert hmin == BUSCH[n], (n, hmin)
        res[f'n={n}'] = {'method': 'exhaustive labelled', 'tc_values': tcs, 'strong_tc_min': smin, 'busch_min_H': hmin,
                         'tc_7_or_21': bad}
        print(f'n={n}: exhaustive, tc spectrum {tcs[:12]}..., strong tc min {smin}, min strong H {hmin} '
              f'(Busch {BUSCH[n]}), tc in {{7,21}}: {bad}', flush=True)
    for n in range(7, nmax + 1):
        lo_aut = max(3, math.ceil(BUSCH[n] / 21))
        if lo_aut % 2 == 0:
            lo_aut += 1
        rows = {}
        for cyc in TYPES[n]:
            p = perm_from_type(n, cyc)
            orbs = pair_orbits(p, n)
            oa = np.array([a for o in orbs for (a, b) in o], np.int64)
            ob = np.array([b for o in orbs for (a, b) in o], np.int64)
            ooff = np.cumsum([0] + [len(o) for o in orbs]).astype(np.int64)
            pmul = math.lcm(*cyc)
            nstrong, ncand, hmin, hits = scan_type(n, oa, ob, ooff, pmul, lo_aut)
            rows[str(cyc)] = {'pair_orbits': len(orbs), 'invariant_tournaments': 1 << len(orbs), 'strong': nstrong,
                              'min_H_strong_invariant': int(hmin) if nstrong else None, 'aut_checked': ncand,
                              'tc_7_or_21': [list(h) for h in hits]}
            print(f'n={n} type {cyc}+fixed: 2^{len(orbs)} invariant, {nstrong} strong (min H {hmin if nstrong else "-"}), '
                  f'{ncand} |Aut| computations, hits {list(hits)} ({time.time() - t0:.0f}s)', flush=True)
            if nstrong:
                assert hmin >= BUSCH[n]
        res[f'n={n}'] = {'method': 'invariant enumeration', 'aut_lower_bound': lo_aut, 'types': rows}
    allhits = [h for k, v in res.items() if k.startswith('n=') for t in v.get('types', {}).values() for h in t['tc_7_or_21']]
    allhits += [b for k, v in res.items() if k.startswith('n=') for b in v.get('tc_7_or_21', [])]
    print(f'TOTAL tc in {{7,21}} for n <= {nmax}: {allhits}  ({time.time() - t0:.0f}s)', flush=True)
    res['total_hits'] = allhits
    json.dump(res, open(os.path.join(HERE, f'tc_fibres_upto{nmax}.json'), 'w'), indent=1)
