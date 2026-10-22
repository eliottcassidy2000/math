"""Two reversed arcs (feedback arc number <= 2): exact formulas for H, and the Mersenne side needs three arcs.
Session opus-2026-10-21-S51; note 05-knowledge/results/repunit_rarity_hyp9357_20261021.md; THM-4748.

Setting: the transitive order 0 < 1 < ... < m-1 (i -> j for i < j), with the arcs (j1, i1), (j2, i2) reversed
(spans >= 2; in a minimum feedback order no back arc has span 1).  Configurations of the two intervals:
  disjoint / touching  -> H = (2^(a) + 1)(2^(b) + 1)                                    (product; touching = Moon chain)
  crossing  i1 < i2 < j1 < j2,  p = i2 - i1 - 1, r = j1 - i2 - 1, q = j2 - j1 - 1:
     H = 1 + 2^(p+r+1) + 2^(q+r+1) - [q=0] 2^p - [p=0] 2^q + 2^(p+q+1) 3^r
  strictly nested  i1 < i2 < j2 < j1,  s = j2 - i2, x = i2 - i1 - 1, y = j1 - j2 - 1:
     H = 1 + 2^(s-1) + 2^(x+y+1) (2^s - 1 + 3^(s-1))
  shared endpoint (i1 = i2 or j1 = j2), outer span L, inner span s, 2 <= s <= L-1:
     H = 1 + 2^(L-1) - 2^(L-1-s) + 2^(s-1) - [s = L-1]
Proof: a Hamiltonian path is a sequence of increasing runs joined by back arcs, each back arc used at most once, and
no run contains a reversed pair consecutively; counting the run assignments of the vertices gives the formulas
(each vertex covered by one interval has 2 run choices, by two intervals 3 choices; the subtracted terms are the
assignments that place a reversed pair consecutively).  The script checks every formula against an exact
Hamiltonian-path DP for every configuration with m <= 14.

Consequences (PROVED, checked here for m <= 14 and argued for all m in THM-4748):
  (C1) the values = 7 (mod 8) of tournaments with fa <= 2 are exactly 15, 23, 3*2^j - 1 (j >= 3), 5*2^j - 1 (j >= 3);
  (C2) hence no tournament with fa <= 2 has H = 2^k - 1 for k >= 5 (and H = 7 never occurs);
  (C3) 2^k + 1 = H(NT_(k+2)) has fa = 1.
  (FAN) the Fibonacci fan: the transitive order on {1, ..., k+1} plus a vertex 0 that loses exactly to a = floor(k/2),
        b = k + 1 - a and a + b = k + 1 (so b - a is 1 or 2) has H = 2^k - 1 (k >= 4) and three back arcs, so for
        k >= 5 the least feedback arc number of H = 2^k - 1 is exactly 3.  Via the one-apex formula (KNOWN,
        05-knowledge/variables/single-core-cycle-count.md): for an apex 0 losing exactly to the vertices of a set
        A of positions >= 2, H - 1 = sum_{s in A} 2^(s-1) - sum_{s<t in A} 2^max(t-s-1, 1); for A = {a, b, a+b} the
        pair terms 2^(b-1), 2^(a-1) cancel the singleton terms of a and b, leaving 2^(a+b-1) - 2.  Checked by DP for
        k <= 14, by the one-apex formula for k <= 400, and the pairwise formula against DP for every A with m <= 12.
Also checked against the full class census n <= 9 (minimum feedback arc number by subset DP).
Prints ALL CHECKS PASSED."""
import sys, os, time, itertools
import numpy as np
import numba as nb

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(cache=True)
def ham(n, out):
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            m = out[v] & ~S
            while m:
                w = 0
                while not (m >> w) & 1:
                    w += 1
                m &= m - 1
                dp[S | (1 << w), w] += c
    t = 0
    for v in range(n):
        t += dp[N - 1, v]
    return t


def tour(m, backs):
    out = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    for (j, i) in backs:
        out[i] &= ~(1 << j)
        out[j] |= 1 << i
    return np.array(out, dtype=np.int64)


def formula(i1, j1, i2, j2):
    """H for back arcs j1 -> i1 and j2 -> i2 with i1 <= i2 (spans >= 2)"""
    if j1 <= i2:                       # disjoint or touching
        return (2 ** (j1 - i1 - 1) + 1) * (2 ** (j2 - i2 - 1) + 1)
    if i1 < i2 and j1 < j2:            # crossing
        p, r, q = i2 - i1 - 1, j1 - i2 - 1, j2 - j1 - 1
        return 1 + 2 ** (p + r + 1) + 2 ** (q + r + 1) - (2 ** p if q == 0 else 0) - (2 ** q if p == 0 else 0) + 2 ** (p + q + 1) * 3 ** r
    # nested: interval 2 inside interval 1, or 1 inside 2 when i1 == i2
    if i1 == i2:
        L, s = max(j1, j2) - i1, min(j1, j2) - i1
        return 1 + 2 ** (L - 1) - 2 ** (L - 1 - s) + 2 ** (s - 1) - (1 if s == L - 1 else 0)
    if j2 == j1:
        L, s = j1 - i1, j2 - i2
        return 1 + 2 ** (L - 1) - 2 ** (L - 1 - s) + 2 ** (s - 1) - (1 if s == L - 1 else 0)
    if j2 < j1:                        # strictly nested
        s, x, y = j2 - i2, i2 - i1 - 1, j1 - j2 - 1
        return 1 + 2 ** (s - 1) + 2 ** (x + y + 1) * (2 ** s - 1 + 3 ** (s - 1))
    raise ValueError


print('(F) two-arc formulas against the exact Hamiltonian-path DP, every configuration with m <= 14')
cnt = 0; bad = 0
vals7 = set()
for m in range(4, 15):
    arcs = [(j, i) for i in range(m) for j in range(i + 2, m)]
    for (a, b) in itertools.combinations(arcs, 2):
        (j1, i1), (j2, i2) = sorted([a, b], key=lambda t: (t[1], t[0]))
        h = int(ham(m, tour(m, [a, b])))
        f = formula(i1, j1, i2, j2)
        cnt += 1
        if h != f:
            bad += 1
            if bad < 5:
                print('    mismatch', m, a, b, h, f)
        if h % 8 == 7:
            vals7.add(h)
check(bad == 0, f'{cnt} two-arc configurations: formula = DP in every case')

print('(C) consequences')
H7 = sorted(vals7)
thabit = {3 * 2 ** j - 1 for j in range(3, 200)} | {5 * 2 ** j - 1 for j in range(3, 200)} | {15, 23}
check(all(h in thabit for h in H7), f'values = 7 mod 8 with two reversed arcs (m <= 14) are all in {{15, 23, 3*2^j-1, 5*2^j-1}}: {H7[:12]} ...')
mers = {2 ** k - 1 for k in range(3, 200)}
check(not any(h in mers for h in H7 if h != 15), 'no 2^k - 1 (k >= 5) among them; 15 = 2^4 - 1 is the only repunit')
# formula-level check of (C1)/(C2) for large spans: scan the formulas symbolically over parameter grids
ok = True
for p in range(0, 40):
    for q in range(0, 40):
        for r in range(0, 40):
            h = 1 + 2 ** (p + r + 1) + 2 ** (q + r + 1) - (2 ** p if q == 0 else 0) - (2 ** q if p == 0 else 0) + 2 ** (p + q + 1) * 3 ** r
            if h % 8 == 7 and h not in thabit:
                ok = False
            if h in mers and h != 15:
                ok = False
for s in range(2, 60):
    for xy in range(0, 60):
        h = 1 + 2 ** (s - 1) + 2 ** (xy + 1) * (2 ** s - 1 + 3 ** (s - 1))
        if (h % 8 == 7 and h not in thabit) or (h in mers and h != 15):
            ok = False
for L in range(3, 80):
    for s in range(2, L):
        h = 1 + 2 ** (L - 1) - 2 ** (L - 1 - s) + 2 ** (s - 1) - (1 if s == L - 1 else 0)
        if (h % 8 == 7 and h not in thabit) or (h in mers and h != 15):
            ok = False
for a in range(1, 80):
    for b in range(1, 80):
        h = (2 ** a + 1) * (2 ** b + 1)
        if (h % 8 == 7 and h not in thabit) or (h in mers and h != 15):
            ok = False
check(ok, 'formula scan (crossing p,q,r < 40; nested s,x+y < 60; shared endpoint L < 80; products a,b < 80): '
          'every value = 7 mod 8 is 15, 23, 3*2^j-1 or 5*2^j-1, and no repunit 2^k-1 with k >= 5 occurs')
check(all(int(ham(k + 2, tour(k + 2, [(k + 1, 0)]))) == 2 ** k + 1 for k in range(1, 13)), '2^k + 1 = H(NT_(k+2)): one reversed arc (k <= 12)')

print('(FAN) the Fibonacci fan realizes 2^k - 1 with three back arcs')


def fan(k):
    a = k // 2; b = k + 1 - a
    m = k + 2
    out = [0] * m
    for i in range(1, m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    for j in range(1, m):
        if j in (a, b, a + b):
            out[j] |= 1 << 0
        else:
            out[0] |= 1 << j
    return m, np.array(out, dtype=np.int64), (a, b, a + b)


def one_apex_pairwise(A):
    A = sorted(A)
    return 1 + sum(2 ** (s - 1) for s in A) - sum(2 ** max(t - s - 1, 1) for s, t in itertools.combinations(A, 2))


ok = True
for k in range(4, 15):
    m, out, A = fan(k)
    ok &= int(ham(m, out)) == 2 ** k - 1 == one_apex_pairwise(A)
check(ok, 'H(fan_k) = 2^k - 1 by DP for 4 <= k <= 14 (the apex loses exactly to floor(k/2), k+1-floor(k/2), k+1)')
check(all(one_apex_pairwise((k // 2, k + 1 - k // 2, k + 1)) == 2 ** k - 1 for k in range(4, 401)), 'and by the one-apex pairwise formula for 4 <= k <= 400')
ok = True
for m in range(3, 13):
    for r in range(0, m - 1):
        for A in itertools.combinations(range(2, m), r):
            out = [0] * m
            for i in range(1, m):
                for j in range(i + 1, m):
                    out[i] |= 1 << j
            for j in range(1, m):
                if j in A:
                    out[j] |= 1
                else:
                    out[0] |= 1 << j
            ok &= int(ham(m, np.array(out, dtype=np.int64))) == one_apex_pairwise(A)
check(ok, 'one-apex pairwise formula H - 1 = sum 2^(s-1) - sum_{s<t} 2^max(t-s-1,1) equals the DP for every in-set A, m <= 12')

print('(K) the class census n <= 9: every class with H = 2^k - 1, k >= 5, has feedback arc number >= 3')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'h_spectrum_multiplicities_20261021'))
import h_multiplicities as HM                                              # noqa: E402


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


data = HM.classes_upto(9)
for n in range(5, 10):
    A = data[n]['A']
    H, st = HM.ham_all(n, A)
    low = [(int(h), int(fa_one(n, A[t]))) for t, h in enumerate(H)]
    reps = sorted(set((h, f) for h, f in low if h in mers and h >= 31))
    minfa = {}
    for h, f in low:
        if h in mers and h >= 31:
            minfa[h] = min(minfa.get(h, 99), f)
    fa2_7 = sorted(set(h for h, f in low if f <= 2 and h % 8 == 7))
    check(all(f >= 3 for f in minfa.values()) and all(h in thabit for h in fa2_7),
          f'n = {n}: least feedback arc number of classes with H = 2^k - 1 (k >= 5): {minfa}; fa <= 2 values = 7 mod 8: {fa2_7}')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
