"""S51 audit c8: the note's section 8 ("refuted or withdrawn") numbers, recomputed from the auditor's census (c2 cache).
(a) ear increments over all non-constant one-vertex extensions of the strong classes of orders 5, 6, 7: the 2-adic
    valuation shares of Delta H, and the transition matrix of H mod 8 (parent residue -> child residue), which the
    note says "mixes in one step" but refuted_leads.py computes without printing;
(b) I(Omega, -1) = 1 - alpha_1 + alpha_2 on the strong classes of orders 5-8 (alpha_k from directed odd cycles of every
    vertex subset; at n <= 8 at most two disjoint odd cycles), its range and whether it vanishes;
(d) the one-apex toy model multiplicities (max <= 8 for 16 <= r < 1100).
Prints ALL CHECKS PASSED or the failures."""
import os, math, time, itertools
from collections import Counter, defaultdict
import numpy as np
import numba as nb
import alib51 as L

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


Z = np.load(os.path.join(SCR, 'census_n9.npz'))
strong = {}
for n in range(3, 10):
    An, Hn, sn = Z[f'A{n}'], Z[f'H{n}'], Z[f'st{n}']          # load each compressed array once
    strong[n] = [(An[t], int(Hn[t])) for t in range(len(Hn)) if sn[t]]


def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    return c


print('(a) ear increments')
for m in (5, 6, 7):
    dist = Counter(); trans = defaultdict(Counter)
    for adj, h in strong[m]:
        ev = L.ear_all(m, adj)
        for x in range(1, (1 << m) - 1):
            dist[min(v2(int(ev[x]) - h), 6)] += 1
            trans[h % 8][int(ev[x]) % 8] += 1
    tot = sum(dist.values())
    sh = {k: round(c / tot, 3) for k, c in sorted(dist.items())}
    print(f'     order {m} -> {m + 1}: v_2(Delta H) shares {sh}')
    rows = {r: {s: round(trans[r][s] / sum(trans[r].values()), 3) for s in (1, 3, 5, 7)} for r in sorted(trans)}
    print(f'     order {m} -> {m + 1}: P(child H mod 8 | parent H mod 8): {rows}')
    tvd = max(0.5 * sum(abs(p - 0.25) for p in row.values()) for row in rows.values())
    print(f'     max total-variation distance of a row from uniform: {tvd:.3f}')
    if m == 5:
        check(sh == {1: 0.511, 2: 0.267, 3: 0.122, 4: 0.067, 5: 0.033}, 'refuted_leads.out (a) order-5 shares reproduced')
    if m == 7:
        check(sh == {1: 0.506, 2: 0.253, 3: 0.119, 4: 0.061, 5: 0.03, 6: 0.031}, 'refuted_leads.out (a) order-7 shares reproduced')
        check(tvd < 0.15, f'"residues mix in one step": every row of the order-7 transition matrix is within TV distance {tvd:.3f} of uniform')


@nb.njit(cache=True)
def odd_cycle_alphas(n, adj):
    """alpha_1 = # directed odd cycles, alpha_2 = # unordered pairs of vertex-disjoint directed odd cycles (n <= 8)"""
    N = 1 << n
    cyc = np.zeros(N, dtype=np.int64)        # directed Hamiltonian cycles of T[S]
    for s in range(n):
        dp = np.zeros((N, n), dtype=np.int64)
        dp[1 << s, s] = 1
        for S in range(N):
            if not (S >> s) & 1 or (S & ((1 << s) - 1)):
                continue
            for v in range(n):
                c = dp[S, v]
                if c == 0:
                    continue
                if v != s and (adj[v] >> s) & 1:
                    cyc[S] += c
                for u in range(s + 1, n):
                    if not (S >> u) & 1 and (adj[v] >> u) & 1:
                        dp[S | (1 << u), u] += c
    a1 = 0
    a2 = 0
    for S in range(N):
        k = 0
        x = S
        while x:
            x &= x - 1
            k += 1
        if k >= 3 and k % 2 == 1 and cyc[S] > 0:
            a1 += cyc[S]
            comp = (N - 1) ^ S
            T = comp
            while T:
                kk = 0
                y = T
                while y:
                    y &= y - 1
                    kk += 1
                if kk >= 3 and kk % 2 == 1 and T > S:
                    a2 += cyc[S] * cyc[T]
                T = (T - 1) & comp
    return a1, a2


print('(b) I(Omega, -1) on strong classes of orders 5-8')
rng = {}
for m in (5, 6, 7, 8):
    vals = []
    for adj, h in strong[m]:
        a1, a2 = odd_cycle_alphas(m, adj)
        assert 1 + 2 * a1 + 4 * a2 == h
        vals.append(1 - a1 + a2)
    rng[m] = (min(vals), max(vals), vals.count(0))
    print(f'     order {m}: I(Omega,-1) in [{min(vals)}, {max(vals)}], zero {vals.count(0)} times ({len(vals)} classes; H = 1 + 2 a1 + 4 a2 checked)')
check(rng == {5: (-6, -3, 0), 6: (-18, -4, 0), 7: (-72, -5, 0), 8: (-162, -3, 0)}, 'refuted_leads.out (b) ranges reproduced; I(Omega,-1) never 0')

print('(d) one-apex toy model')
mu = Counter()
for Lw in range(2, 16):
    for w in itertools.product((0, 1), repeat=Lw):
        if w[0] != 0 or w[-1] != 1:
            continue
        A = {i + 1 for i in range(Lw) if w[i] == 1}       # chain positions beating the apex? use the auditor's closed form
        # the session's word: w[a] = 1 for a position that the apex loses to (r counts pairs b < a with w[b] = 0, w[a] = 1)
        r = (L.one_apex_paths(Lw + 1, {i + 1 for i in range(Lw) if w[i] == 1}) - 1) // 2
        mu[r] += 1
check(max(mu[r] for r in range(16, 1100)) <= 8, f'one-apex word multiplicities <= 8 for 16 <= r < 1100 (max {max(mu[r] for r in range(16, 1100))})')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
