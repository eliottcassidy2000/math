#!/usr/bin/env python3
"""Independent audit of the lenses note (05-knowledge/results/lenses_tm_ca_fractran_tournaments_20261101.md).

  N1  PRIMEGAME's 14 fractions: per-state groups by the state prime of the denominator; the state-free group; whether
      PRIMEGAME is "one-hot" and whether its realisable co-applicable sets lie inside single groups (section 4 and
      results item 5 of the note).  Run from 2 for 200,000 steps, collect every applicable set.
  N2  GOLDBACHGAME (THM-4856, list read as data from goldbachgame.out): every fraction has exactly one state prime in
      numerator and denominator (so the per-state claim applies there).
  N3  H-spectrum section: H multiplicative over strong components (all tournaments n <= 6); minimum H over strong
      m-tournaments m = 3..6 equals Moon's chain value; H(strong m-tournament) >= m.
usage (from the repo root): python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit_lenses/audit_note.py
"""
import itertools
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []


def report(name, ok, detail=''):
    print(f'  [{"OK" if ok else "MISMATCH"}] {name}' + (f' -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


def factor(m):
    out, p = {}, 2
    while p * p <= m:
        while m % p == 0:
            out[p] = out.get(p, 0) + 1
            m //= p
        p += 1
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


PRIMEGAME = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13),
             (13, 11), (15, 14), (15, 2), (55, 1)]
STATES = {11, 13, 17, 19, 23, 29}


def main():
    print('N1  PRIMEGAME priority structure')
    groups = defaultdict(list)
    for i, (a, b) in enumerate(PRIMEGAME):
        st = [p for p in factor(b) if p in STATES]
        assert len(st) <= 1
        groups[st[0] if st else 'free'].append(i)
    sizes = {k: len(v) for k, v in groups.items()}
    report('group sizes 17:3, 13:2, 11:2, 19:2, 23:1, 29:1, state-free:3', sizes == {17: 3, 13: 2, 11: 2, 19: 2, 23: 1,
                                                                                     29: 1, 'free': 3}, f'{sizes}')
    n = 2
    appl_sets = Counter()
    zero_state = 0
    mixed = Counter()
    for step in range(200_000):
        app = tuple(i for i, (a, b) in enumerate(PRIMEGAME) if n % b == 0)
        appl_sets[app] += 1
        st = [p for p in STATES if n % p == 0]
        assert len(st) <= 1 and all(n % (p * p) for p in st)
        if not st:
            zero_state += 1
        gs = {('free' if i in groups['free'] else [k for k, v in groups.items() if i in v][0]) for i in app}
        if len(gs) > 1:
            mixed[app] += 1
        i = app[0]
        n = n // PRIMEGAME[i][1] * PRIMEGAME[i][0]
    always = all(13 in s for s in appl_sets)
    report('PRIMEGAME is NOT one-hot: reachable integers with no state prime occur, and 55/1 (index 13) is applicable '
           'at EVERY integer, so it is co-applicable with every state group', zero_state > 0 and always and len(mixed) > 0,
           f'{zero_state} of 200,000 visited integers have no state prime; {len(mixed)} distinct realised applicable '
           f'sets straddle a state group and the state-free group')
    ex = sorted(mixed, key=lambda s: -mixed[s])[:4]
    print('    most frequent straddling applicable sets (fraction indices, 0-based; 8 = 1/17, 12 = 15/2, 13 = 55/1):',
          [(s, mixed[s]) for s in ex])
    # the order between groups matters: putting the state-free group first changes the run immediately
    order2 = [PRIMEGAME[i] for i in groups['free']] + [PRIMEGAME[i] for i in range(14) if i not in groups['free']]
    seq1, seq2, n1, n2 = [], [], 2, 2
    for _ in range(10):
        for a, b in PRIMEGAME:
            if n1 % b == 0:
                n1 = n1 // b * a
                break
        for a, b in order2:
            if n2 % b == 0:
                n2 = n2 // b * a
                break
        seq1.append(n1)
        seq2.append(n2)
    report('the inter-group orientation matters (state-free group first gives a different run)', seq1 != seq2,
           f'PRIMEGAME {seq1[:5]} vs reordered {seq2[:5]}')

    print('N2  GOLDBACHGAME one-hot check (fractions read from goldbachgame.out)')
    txt = open(os.path.join(HERE, '..', 'goldbachgame.out')).read().split('\n')[1]
    fr = [tuple(map(int, f.split('/'))) for f in re.findall(r'\d+/\d+', txt)]
    regs = {2, 3, 5, 7, 31, 37, 41}
    ok = len(fr) == 37
    for a, b in fr:
        sa = [p for p in factor(a) if p not in regs]
        sb = [p for p in factor(b) if p not in regs]
        ok &= len(sa) == 1 and len(sb) == 1 and factor(a)[sa[0]] == 1 and factor(b)[sb[0]] == 1
    report('GOLDBACHGAME: 37 fractions, each with exactly one non-register (state) prime, exponent 1, in numerator and '
           'denominator', ok)

    print('N3  H-spectrum section: multiplicativity over strong components, strong minima')

    def ham_paths(adj, verts):
        verts = list(verts)
        k = len(verts)
        idx = {v: i for i, v in enumerate(verts)}
        dp = [[0] * k for _ in range(1 << k)]
        for i in range(k):
            dp[1 << i][i] = 1
        for mask in range(1 << k):
            for i in range(k):
                c = dp[mask][i]
                if c:
                    for j in range(k):
                        if not mask >> j & 1 and adj[verts[i]][verts[j]]:
                            dp[mask | 1 << j][j] += c
        return sum(dp[(1 << k) - 1])

    def sccs(adj, n):
        reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
        for k in range(n):
            for i in range(n):
                if reach[i][k]:
                    for j in range(n):
                        if reach[k][j]:
                            reach[i][j] = True
        comps, seen = [], set()
        for v in range(n):
            if v not in seen:
                c = [w for w in range(n) if reach[v][w] and reach[w][v]]
                seen |= set(c)
                comps.append(c)
        return comps

    ok_mult = True
    strong_min = {}
    for n in range(3, 7):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for bits in range(1 << len(pairs)):
            adj = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                if bits >> k & 1:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
            H = ham_paths(adj, range(n))
            comps = sccs(adj, n)
            prod = 1
            for c in comps:
                prod *= ham_paths(adj, c)
            ok_mult &= prod == H
            if len(comps) == 1:
                strong_min[n] = min(strong_min.get(n, H), H)
    moon = {3: 3, 4: 5, 5: 9, 6: 15}
    report('H(T) = product of H over strong components, every labelled tournament n = 3..6', ok_mult)
    report('minimum H over strong m-tournaments = Moon chain value 3, 5, 9, 15 (m = 3..6), all >= m',
           strong_min == moon and all(strong_min[m] >= m for m in moon), f'{strong_min}')
    print()
    print('ALL AUDIT CHECKS AGREE' if not FAILS else f'MISMATCHES: {FAILS}')
    return 0 if not FAILS else 1


if __name__ == '__main__':
    sys.exit(main())
