"""S51 audit c1: THM-4748 statements 1-5 by the auditor's own code (alib51; no session code imported).
(i)   every pair of reversed arcs of the transitive order, ALL spans (span 1 included), m <= 16: exact H by DP; the
      formulas of statement 1 (typed from the theorem text) against the DP for every pair with spans >= 2;
(ii)  statement 2: the values = 7 (mod 8) over ALL tournaments with fa <= 2 on m <= 16 vertices (every reversal of at
      most two arcs, any spans); every listed value occurs; exact-integer formula scan far beyond m = 16;
(iii) statement 4: the Fibonacci fan by DP (k <= 18: H, strong, three back arcs, fa = 3) and by the auditor's own
      closed form for one-apex tournaments (k <= 3000); the one-apex pairwise formula against DP and the closed form;
      which one-apex in-sets give 2^k - 1;
(iv)  statement 5 (NT_(k+2)) and bof's 2016 crossing value 3*2^(n-4) + 3 (MathOverflow 234029) as a special case.
Prints ALL CHECKS PASSED or the failures."""
import itertools, time
import numpy as np
import numba as nb
import alib51 as L

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(parallel=True, cache=True)
def ham_batch(m, A):
    out = np.zeros(A.shape[0], dtype=np.int64)
    for t in nb.prange(A.shape[0]):
        out[t] = L.ham(m, A[t])
    return out


def build(m, backs):
    return L.reverse_arcs(m, backs)


# ---------------------------------------------------------------- the formulas, typed from the THM-4748 text
def formula_text(m, a1, a2):
    """a = (y, x) back arc y -> x, spans >= 2.  Returns (configuration, H) per THM-4748 statement 1."""
    (y1, x1), (y2, x2) = sorted([a1, a2], key=lambda t: (t[1], t[0]))   # x1 <= x2
    if x1 == x2 or y1 == y2:                                            # shared endpoint
        L_ = max(y1 - x1, y2 - x2); s = min(y1 - x1, y2 - x2)
        return 'shared', 1 + 2 ** (L_ - 1) - 2 ** (L_ - 1 - s) + 2 ** (s - 1) - (1 if s == L_ - 1 else 0)
    if y1 <= x2:                                                        # disjoint or touching
        a = (y1 - x1) - 1; b = (y2 - x2) - 1                            # a + 1, b + 1 are the spans
        return ('touching' if y1 == x2 else 'disjoint'), (2 ** a + 1) * (2 ** b + 1)
    if x1 < x2 < y1 < y2:                                               # crossing
        p = x2 - x1 - 1; r = y1 - x2 - 1; q = y2 - y1 - 1
        return 'crossing', (1 + 2 ** (p + r + 1) + 2 ** (q + r + 1) - (2 ** p if q == 0 else 0)
                            - (2 ** q if p == 0 else 0) + 2 ** (p + q + 1) * 3 ** r)
    if x1 < x2 < y2 < y1:                                               # strictly nested
        s = y2 - x2; x = x2 - x1 - 1; y = y1 - y2 - 1
        return 'nested', 1 + 2 ** (s - 1) + 2 ** (x + y + 1) * (2 ** s - 1 + 3 ** (s - 1))
    raise ValueError((a1, a2))


thabit = lambda h: any(h == 3 * 2 ** j - 1 or h == 5 * 2 ** j - 1 for j in range(3, 400))
listed = lambda h: h in (15, 23) or thabit(h)
mers = {2 ** k - 1 for k in range(2, 400)}

print('(i) two reversed arcs, every pair, all spans, m <= 16: DP against the statement-1 formulas (spans >= 2)')
stats = {}
vals_fa2 = {}           # m -> set of H over all reversals of <= 2 arcs (= all tournaments with fa <= 2 on m vertices)
bad = []
nconf = 0
for m in range(2, 17):
    arcs = [(y, x) for x in range(m) for y in range(x + 1, m)]
    pairs = list(itertools.combinations(arcs, 2))
    A = np.array([build(m, list(pr)) for pr in pairs] + [build(m, [a]) for a in arcs] + [L.transitive(m)], dtype=np.int64)
    Hs = ham_batch(m, A)
    vals_fa2[m] = set(int(h) for h in Hs)
    for idx, pr in enumerate(pairs):
        (y1, x1), (y2, x2) = pr
        if y1 - x1 < 2 or y2 - x2 < 2:
            continue
        conf, f = formula_text(m, pr[0], pr[1])
        nconf += 1
        stats[conf] = stats.get(conf, 0) + 1
        if f != int(Hs[idx]):
            bad.append((m, pr, conf, int(Hs[idx]), f))
    # single arcs: H = 2^(span-1) + 1
    for j, a in enumerate(arcs):
        h1 = int(Hs[len(pairs) + j])
        if a[0] - a[1] >= 2 and h1 != 2 ** (a[0] - a[1] - 1) + 1:
            bad.append((m, a, 'single', h1))
    print(f'     m = {m:2d}: {len(pairs)} pairs ({sum(1 for p in pairs if p[0][0]-p[0][1] >= 2 and p[1][0]-p[1][1] >= 2)} with spans >= 2) [{time.time() - T0:.0f}s]', flush=True)
check(not bad, f'{nconf} two-arc configurations with spans >= 2, m <= 16: the THM-4748 (1) formula equals the DP in every case '
               f'(by configuration: {dict(sorted(stats.items()))}); single arcs give 2^(span-1) + 1; mismatches: {bad[:5]}')
check(nconf == sum(len(list(itertools.combinations([(y, x) for x in range(m) for y in range(x + 2, m)], 2))) for m in range(2, 17)),
      'the five configurations exhaust all pairs of back arcs with spans >= 2 (no ValueError, every pair classified)')
bnd = [(m, pr) for (m, pr, *_ ) in bad]
# boundary cases explicitly present
def present(cond):
    for m in range(4, 17):
        arcs = [(y, x) for x in range(m) for y in range(x + 2, m)]
        for pr in itertools.combinations(arcs, 2):
            (y1, x1), (y2, x2) = sorted(pr, key=lambda t: (t[1], t[0]))
            if cond(x1, y1, x2, y2):
                return True
    return False
check(present(lambda x1, y1, x2, y2: x1 < x2 < y1 < y2 and x2 - x1 - 1 == 0 and y2 - y1 - 1 == 0) and
      present(lambda x1, y1, x2, y2: x1 < x2 < y1 < y2 and x2 - x1 - 1 == 0 and y2 - y1 - 1 > 0) and
      present(lambda x1, y1, x2, y2: x1 < x2 < y1 < y2 and y1 - x2 - 1 == 0) and
      present(lambda x1, y1, x2, y2: (x1 == x2 or y1 == y2) and min(y1 - x1, y2 - x2) == max(y1 - x1, y2 - x2) - 1) and
      present(lambda x1, y1, x2, y2: x1 < x2 < y2 < y1 and y2 - x2 == 2 and x2 - x1 == 1 and y1 - y2 == 1),
      'the boundary cases p = q = 0, p = 0 < q, r = 0, s = L - 1, nested s = 2 with x = y = 0 are all among the checked configurations')

print('(ii) statement 2: H = 7 (mod 8) over ALL tournaments with fa <= 2 (every reversal of <= 2 arcs, any spans), m <= 16')
allv = set().union(*vals_fa2.values())
v7 = sorted(h for h in allv if h % 8 == 7)
check(all(listed(h) for h in v7), f'every value = 7 mod 8 is 15, 23, 3*2^j-1 or 5*2^j-1 (j >= 3): {v7}')
check(not any(h in mers and h != 15 for h in v7) and 7 not in allv, 'no 2^k - 1 with k >= 3 except 15 occurs with fa <= 2 (m <= 16)')
need = [15, 23] + [3 * 2 ** j - 1 for j in range(3, 14)] + [5 * 2 ** j - 1 for j in range(3, 13)]
miss = [h for h in need if h not in allv]
check(not miss, f'every listed value 3*2^j-1 (3 <= j <= 13) and 5*2^j-1 (3 <= j <= 12) occurs at m <= 16 (missing: {miss})')
check(all(3 * 2 ** j - 1 in vals_fa2[j + 3] and 3 * 2 ** j - 1 not in vals_fa2[j + 2] for j in range(3, 14)) and all(5 * 2 ** j - 1 in vals_fa2[j + 4] and 5 * 2 ** j - 1 not in vals_fa2[j + 3] for j in range(3, 13)),
      'the lines first occur at the predicted orders: 3*2^j-1 at m = j + 3 (shared endpoint, L = m - 1, s = L-1), 5*2^j-1 at m = j + 4 (s = L-2)')
# exact-integer scan of the formulas far beyond m = 16
ok = True; worst = []
for p in range(0, 61):
    for q in range(0, 61):
        for r in range(0, 61):
            h = 1 + 2 ** (p + r + 1) + 2 ** (q + r + 1) - (2 ** p if q == 0 else 0) - (2 ** q if p == 0 else 0) + 2 ** (p + q + 1) * 3 ** r
            if h % 8 == 7 and not listed(h):
                ok = False; worst.append(('crossing', p, q, r))
for s in range(2, 101):
    for xy in range(0, 101):
        h = 1 + 2 ** (s - 1) + 2 ** (xy + 1) * (2 ** s - 1 + 3 ** (s - 1))
        if h % 8 == 7 and not listed(h):
            ok = False; worst.append(('nested', s, xy))
for Lh in range(3, 301):
    for s in range(2, Lh):
        h = 1 + 2 ** (Lh - 1) - 2 ** (Lh - 1 - s) + 2 ** (s - 1) - (1 if s == Lh - 1 else 0)
        if h % 8 == 7 and not listed(h):
            ok = False; worst.append(('shared', Lh, s))
for a in range(1, 301):
    for b in range(1, 301):
        h = (2 ** a + 1) * (2 ** b + 1)
        if h % 8 == 7 and h != 15:
            ok = False; worst.append(('product', a, b))
check(ok, f'exact formula scan (crossing p,q,r <= 60; nested s, x+y <= 100; shared L <= 300; products a,b <= 300): '
          f'H = 7 mod 8 only on 15, 23, 3*2^j-1, 5*2^j-1 (counterexamples: {worst[:5]})')
# the mod-8 case analysis exactly as written in the proof of statement 2
ok = True
for q in range(1, 61):
    for r in range(0, 61):
        hm1 = 2 ** (r + 1) + 2 ** (q + r + 1) - 2 ** q + 2 ** (q + 1) * 3 ** r
        if (hm1 % 8 == 6) != ((q, r) in ((2, 0), (1, 1))):
            ok = False
for r in range(0, 61):
    if (2 ** (r + 2) - 2 + 2 * 3 ** r) % 8 not in (0, 4):
        ok = False
for xy in range(0, 61):
    if ((2 + 3 * 2 ** (xy + 2)) % 8 == 6) != (xy == 0):
        ok = False
for Lh in range(5, 301):
    if (3 * 2 ** (Lh - 3) + 3) % 8 == 7 and Lh != 5:
        ok = False
check(ok, 'each step of the written mod-8 case analysis (p = 0 < q: only (q,r) = (2,0), (1,1); p = q = 0: H-1 = 0,4 mod 8; '
          'nested s = 2: only x = y = 0; shared s = 2, t >= 2: only L = 5) holds on the scanned range')
check(23 == 3 * 2 ** 3 - 1, 'remark: 23 = 3*2^3 - 1 is already the j = 3 point of the Thabit line (listing it separately is redundant)')

print('(iii) statement 4: the Fibonacci fan')


def fan(k):
    a = k // 2; b = k + 1 - a
    return a, b, L.one_apex(k + 2, {a, b, a + b})


ok = True
rows = []
for k in range(4, 19):
    a, b, adj = fan(k)
    m = k + 2
    h = int(L.ham(m, adj))
    st = bool(L.is_strong_mask(m, adj, (1 << m) - 1))
    nback = sum(1 for i in range(m) for j in range(i + 1, m) if (int(adj[j]) >> i) & 1)
    f = int(L.fas(m, adj))
    rows.append((k, h == 2 ** k - 1, st, nback, f))
    ok &= h == 2 ** k - 1 and st and nback == 3 and (f == 3 if k >= 5 else f <= 3) and b - a in (1, 2) and a >= 2
check(ok, f'F_k for 4 <= k <= 18 by DP: H = 2^k - 1, strong, exactly 3 back arcs in the order 0 < 1 < ... < k+1, b - a in {{1, 2}}, '
          f'and fa(F_k) = 3 for k >= 5 (fa(F_4) = {rows[0][4]})')
check(all(L.one_apex_paths(k + 2, {k // 2, k + 1 - k // 2, k + 1}) == 2 ** k - 1 for k in range(4, 3001)),
      'H(F_k) = 2^k - 1 for 4 <= k <= 3000 by the auditor\'s path-decomposition closed form (independent of OCF)')
check(all(L.one_apex_pairwise({k // 2, k + 1 - k // 2, k + 1}) == 2 ** k - 1 for k in range(4, 3001)),
      'and by the THM-4748 (4) pairwise formula, 4 <= k <= 3000')
ok = True
for m in range(3, 14):
    As = [A for r in range(0, m - 1) for A in itertools.combinations(range(2, m), r)]
    arr = np.array([L.one_apex(m, set(A)) for A in As], dtype=np.int64)
    Hs = ham_batch(m, arr)
    ok &= all(int(h) == L.one_apex_pairwise(A) for h, A in zip(Hs, As))
check(ok, 'one-apex pairwise formula = DP for every in-set A subset of {2, ..., m-1}, m <= 13')
ok = True
for m in range(3, 17):
    for r in range(0, m - 1):
        for A in itertools.combinations(range(2, m), r):
            ok &= L.one_apex_paths(m, A) == L.one_apex_pairwise(A)
check(ok, 'pairwise formula = path-decomposition closed form for every A subset of {2, ..., m-1}, m <= 16')
off = [A for m in range(3, 9) for r in range(1, m) for A in itertools.combinations(range(1, m), r)
       if 1 in A and L.one_apex_paths(m, A) != L.one_apex_pairwise(A) - 1]
check(not off, 'remark: with 1 in A the pairwise formula is off by exactly +1 (so the restriction A subset of {2..m-1} in the text is needed)')
# which one-apex in-sets realize 2^k - 1 (m <= 18), by size of A
sols = {}
for m in range(3, 19):
    for r in range(1, m - 1):
        for A in itertools.combinations(range(2, m), r):
            h = L.one_apex_pairwise(A)
            if h in mers and h >= 31:
                sols.setdefault(h, []).append(A)
sizes = {h: sorted(set(len(A) for A in v)) for h, v in sorted(sols.items())}
three = {h: sorted(set(A for A in v if len(A) == 3)) for h, v in sorted(sols.items())}
fib = all(all(A[0] + A[1] == A[2] and A[1] - A[0] in (1, 2) for A in v) for v in three.values())
print(f'     one-apex in-sets A (subset of {{2..m-1}}, m <= 18) with H = 2^k - 1, k >= 5: sizes of A by H: {sizes}')
print(f'     the distinct |A| = 3 solutions: {three}')
check(all(min(s) >= 3 for s in sizes.values()), 'no one-apex in-set of size <= 2 gives 2^k - 1 (k >= 5), m <= 18 (consistent with statement 3)')
check(fib, 'every |A| = 3 one-apex solution of H = 2^k - 1 (m <= 18) is a Fibonacci triple {a, b, a+b} with b - a in {1, 2} (remark "exactly the Fibonacci relation")')

print('(iv) statement 5 and a 2016 special case of the crossing formula')
check(all(int(L.ham(k + 2, L.reverse_arcs(k + 2, [(k + 1, 0)]))) == 2 ** k + 1 for k in range(1, 19)), 'H(NT_(k+2)) = 2^k + 1, k <= 18 (one back arc)')
ok = all(int(L.ham(n, L.reverse_arcs(n, [(2, 0), (n - 1, 1)]))) == 3 * 2 ** (n - 4) + 3 == formula_text(n, (2, 0), (n - 1, 1))[1] for n in range(5, 17))
check(ok, 'bof (MathOverflow 234029, 2016): reversing v1v3 and v2vn gives 3*2^(n-4) + 3 = the crossing formula with p = r = 0, q = n - 4 (n = 5..16)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
