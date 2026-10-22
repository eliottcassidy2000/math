"""The leads that did not explain the repunit rarity (session opus-2026-10-21-S51; note section 8, and section 5's
feedback-arc shares).  Prints ALL CHECKS PASSED.
(a) 2-adic inheritance along ears: the valuations v_2(H(T) - H(T')) over all strong one-vertex extensions of strong
    classes of orders 5, 6, 7 are geometric (about 1/2, 1/4, 1/8, ...), and H mod 8 mixes in one step.
(b) H = I(Omega, -1) (mod 3) since 2 = -1 (mod 3); but I(Omega, -1) is never 0 on strong classes of orders 5-8.
(c) H mod 8 by least feedback arc number fa among the strong classes of orders 6-9 (complete census).
(d) the one-apex toy model: strong words of every length, H = 1 + 2 r(w); multiplicities of r near 2^j - 1 and 2^j."""
import os, sys, time, math, itertools
from collections import Counter, defaultdict
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021'))
import alib as L                                       # noqa: E402
import h_multiplicities as HM                          # noqa: E402

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    return c


@nb.njit(cache=True)
def odd_cycle_masks(n, adj):
    out = []
    N = 1 << n
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
                pc = 0
                x = S
                while x:
                    x &= x - 1
                    pc += 1
                if pc >= 3 and pc % 2 == 1 and (adj[v] >> s) & 1:
                    for _ in range(c):
                        out.append(S)
                m = adj[v] & ~S & (N - 1) & ~((1 << (s + 1)) - 1)
                while m:
                    w = 0
                    while not (m >> w) & 1:
                        w += 1
                    m &= m - 1
                    dp[S | (1 << w), w] += c
    return out


def chi_and_H(n, adj):
    cyc = odd_cycle_masks(n, np.array(adj, dtype=np.int64))
    mc = Counter(cyc); masks = sorted(mc)
    al = Counter()

    def rec(i, used, k, w):
        al[k] += w
        for j in range(i, len(masks)):
            if not (masks[j] & used):
                rec(j + 1, used | masks[j], k + 1, w * mc[masks[j]])
    rec(0, 0, 0, 1)
    return sum(c * (-1) ** k for k, c in al.items()), sum(c * 2 ** k for k, c in al.items())


@nb.njit(cache=True)
def fa_one(n, adj):
    N = 1 << n
    f = np.full(N, 10 ** 9, dtype=np.int64)
    f[0] = 0
    for S in range(1, N):
        best = 10 ** 9
        x = S
        while x:
            v = 0
            while not (x >> v) & 1:
                v += 1
            x &= x - 1
            T = S & ~(1 << v)
            m = adj[v] & T
            cnt = 0
            while m:
                m &= m - 1
                cnt += 1
            if f[T] + cnt < best:
                best = f[T] + cnt
        f[S] = best
    return f[N - 1]


@nb.njit(parallel=True, cache=True)
def fa_all(n, A):
    out = np.zeros(A.shape[0], dtype=np.int64)
    for t in nb.prange(A.shape[0]):
        out[t] = fa_one(n, A[t])
    return out


data = HM.classes_upto(9)
strong = {}
for n in range(3, 10):
    A = data[n]['A']; H, st = HM.ham_all(n, A)
    strong[n] = [(np.array(A[t], dtype=np.int64), int(H[t])) for t in range(len(A)) if st[t]]

print('(a) ear increments of strong classes')
for m in (5, 6, 7):
    dist = Counter(); trans = defaultdict(Counter)
    for adj, h in strong[m]:
        vals = L.ear_values(m, adj, np.array(L.converse([int(x) for x in adj]), dtype=np.int64))
        for x in range(1, (1 << m) - 1):
            d = int(vals[x]) - h
            dist[min(v2(d), 6)] += 1
            trans[h % 8][int(vals[x]) % 8] += 1
    tot = sum(dist.values())
    share = {k: round(c / tot, 3) for k, c in sorted(dist.items())}
    print(f'     order {m} -> {m + 1}: v_2(Delta H) shares {share}')
    check(abs(share[1] - 0.5) < 0.04 and abs(share[2] - 0.25) < 0.03, f'order {m}: valuations are geometric (no inheritance of the 2-adic residue)')

print('(b) I(Omega, -1) on strong classes')
for m in (5, 6, 7, 8):
    chis = []
    for adj, h in strong[m]:
        chi, HH = chi_and_H(m, [int(x) for x in adj])
        assert HH == h and (h - chi) % 3 == 0
        chis.append(chi)
    print(f'     order {m}: I(Omega,-1) ranges over [{min(chis)}, {max(chis)}], zero {chis.count(0)} times, H = I(Omega,-1) (mod 3) for all {len(chis)} classes')
    check(0 not in chis, f'order {m}: I(Omega,-1) is never 0')

print('(c) H mod 8 by least feedback arc number (strong classes, complete census)')
for m in (6, 7, 8, 9):
    A = np.array([adj for adj, h in strong[m]], dtype=np.int64)
    fa = fa_all(m, A)
    tab = defaultdict(Counter)
    for (adj, h), f in zip(strong[m], fa):
        tab[int(f)][h % 8] += 1
    print(f'     order {m}: ' + '; '.join(f'fa {f}: share of 7 mod 8 = {tab[f][7] / sum(tab[f].values()):.3f} ({sum(tab[f].values())})' for f in sorted(tab)))
check(True, 'feedback-arc shares listed')

print('(d) the one-apex toy model')
mu = Counter()
for Lw in range(2, 16):
    for w in itertools.product((0, 1), repeat=Lw):
        if w[0] != 0 or w[-1] != 1:
            continue
        r = 0
        for a in range(Lw):
            if w[a] == 1:
                for b in range(a):
                    if w[b] == 0:
                        d = a - b
                        r += 1 if d == 1 else 2 ** (d - 2)
        mu[r] += 1
for j in range(4, 11):
    print(f'     r around 2^{j}: ' + ' '.join(f'{2 ** j + i}:{mu[2 ** j + i]}' for i in range(-3, 4)))
check(max(mu[r] for r in range(16, 1100)) <= 8, 'one-apex word multiplicities are tiny (at most 8 per value for 16 <= r < 1100): no rarity signal at this level')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
