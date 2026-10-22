"""S51 audit c2: the auditor's own class census n <= 9 (one-vertex extensions of every (n-1)-class, own canonical form),
with H, strongness, least feedback arc number fa, |Aut| and self-converseness for every class.  Checks:
  THM-4748 (6) least fa of the classes with H = 31, 63, 127, 255 (and 33, 65); fa <= 2 classes reach 7 mod 8 only on
  the listed lines; in every minimum ordering of every class with fa <= 3 (n <= 8) no back arc has span 1;
  THM-4749 (6) the n = 8 row of the table; the S51 outputs that come from the census: repunit_rarity.out (O) octave
  densities and the n = 8 multiplicity/fa regression, refuted_leads.out (c) fa shares, class_residues.out mod-3 counts
  and labelled shares, the note's order-7 spectrum claims (holes, thin self-converse values, S_7(63) = 0, S_7(127) = 1).
Writes census_n9.npz to the scratchpad (cache for c4/c6).  Prints ALL CHECKS PASSED or the failures."""
import os, sys, time, math, itertools
from collections import Counter, defaultdict
import numpy as np
import numba as nb
import alib51 as L

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'
os.makedirs(SCR, exist_ok=True)


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(parallel=True, cache=True)
def extend_codes(n1, parents):
    P = parents.shape[0]
    X = 1 << n1
    out = np.empty((P * X, 2), dtype=np.int64)
    for t in nb.prange(P):
        for x in range(X):
            ch = L.extend(n1, parents[t], x)
            w1, w0, a = L.canon(n1 + 1, ch)
            out[t * X + x, 0] = w1
            out[t * X + x, 1] = w0
    return out


@nb.njit(parallel=True, cache=True)
def props(n, codes):
    K = codes.shape[0]
    H = np.zeros(K, dtype=np.int64); st = np.zeros(K, dtype=np.int64); fa = np.zeros(K, dtype=np.int64)
    aut = np.zeros(K, dtype=np.int64); sc = np.zeros(K, dtype=np.int64)
    A = np.zeros((K, n), dtype=np.int64)
    full = (1 << n) - 1
    for t in nb.prange(K):
        adj = L.decode(n, codes[t, 0], codes[t, 1])
        for v in range(n):
            A[t, v] = adj[v]
        H[t] = L.ham(n, adj)
        st[t] = 1 if (n == 1 or L.is_strong_mask(n, adj, full)) else 0
        fa[t] = L.fas(n, adj)
        w1, w0, a = L.canon(n, adj)
        aut[t] = a
        cv = np.empty(n, dtype=np.int64)
        for v in range(n):
            cv[v] = full ^ adj[v] ^ (1 << v)
        c1, c0, _ = L.canon(n, cv)
        sc[t] = 1 if (c1 == w1 and c0 == w0) else 0
    return H, st, fa, aut, sc, A


@nb.njit(cache=True)
def min_orders_span1(n, adj, f):
    """number of orderings attaining fa = f that contain a back arc of span 1 (Heap's algorithm over all n! orders)"""
    p = np.arange(n)
    c = np.zeros(n, dtype=np.int64)
    bad = 0
    tot = 0
    i = 0
    while True:
        # evaluate p
        nb_ = 0
        span1 = False
        for a in range(n):
            for b in range(a + 1, n):
                if (adj[p[b]] >> p[a]) & 1:
                    nb_ += 1
                    if b == a + 1:
                        span1 = True
        if nb_ == f:
            tot += 1
            if span1:
                bad += 1
        # next permutation (Heap)
        while i < n:
            if c[i] < i:
                if i % 2 == 0:
                    t = p[0]; p[0] = p[i]; p[i] = t
                else:
                    t = p[c[i]]; p[c[i]] = p[i]; p[i] = t
                c[i] += 1
                i = 0
                break
            else:
                c[i] = 0
                i += 1
        if i >= n:
            break
    return bad, tot


print('(C) class census n <= 9 by one-vertex extension and the auditor\'s canonical form')
codes = {1: np.zeros((1, 2), dtype=np.int64)}
c1 = L.canon(1, np.zeros(1, dtype=np.int64)); codes[1][0] = (c1[0], c1[1])
data = {}
for n in range(1, 10):
    if n > 1:
        par = data[n - 1]['A']
        ext = extend_codes(n - 1, par)
        codes[n] = np.unique(ext, axis=0)
    H, st, fa, aut, sc, A = props(n, codes[n])
    data[n] = dict(H=H, st=st, fa=fa, aut=aut, sc=sc, A=A)
    print(f'     n = {n}: {len(H)} classes, {int(st.sum())} strong, {int(sc.sum())} self-converse [{time.time() - T0:.0f}s]', flush=True)
cnt = [len(data[n]['H']) for n in range(1, 10)]
check(cnt == [1, 1, 2, 4, 12, 56, 456, 6880, 191536], f'class counts n = 1..9: {cnt} (A000568)')
scnt = [int(data[n]['st'].sum()) for n in range(3, 10)]
check(scnt == [1, 1, 6, 35, 353, 6008, 178133], f'strong class counts n = 3..9: {scnt} (A051337: 1, 1, 6, 35, 353, 6008, 178133)')
lab = [sum(math.factorial(n) // int(a) for a in data[n]['aut']) for n in range(1, 10)]
check(lab == [2 ** (n * (n - 1) // 2) for n in range(1, 10)], 'sum of n!/|Aut| = 2^C(n,2) for every n <= 9 (orbit-stabilizer: |Aut| exact)')
check(all(int(h) % int(a) == 0 for n in range(1, 10) for h, a in zip(data[n]['H'], data[n]['aut'])), '|Aut| divides H for every class n <= 9')
np.savez_compressed(os.path.join(SCR, 'census_n9.npz'), **{f'{k}{n}': v for n in data for k, v in data[n].items()})

print('(K) THM-4748 (6) and statement 2 on the census')
mers = [2 ** k - 1 for k in range(2, 13)]
minfa = defaultdict(lambda: 99)
for n in range(1, 10):
    for h, f in zip(data[n]['H'], data[n]['fa']):
        minfa[int(h)] = min(minfa[int(h)], int(f))
print(f'     least fa over all classes n <= 9: ' + ', '.join(f'{h}: {minfa[h]}' for h in (15, 17, 31, 33, 63, 65, 127, 129, 255, 257, 511, 513, 1023, 2047) if minfa[h] < 99))
check(all(minfa[h] == 3 for h in (31, 63, 127, 255)), 'least fa = 3 for H = 31, 63, 127, 255 (classes n <= 9)')
check(minfa[33] == 1 and minfa[65] == 1, 'least fa = 1 for H = 33, 65')
perm = {}
for n in range(1, 10):
    for h, f in zip(data[n]['H'], data[n]['fa']):
        if int(h) in (31, 63, 127, 255, 511, 1023, 2047):
            perm.setdefault((n, int(h)), Counter())[int(f)] += 1
print('     fa distribution of the classes with H = 2^k - 1 (k >= 5), by n: ' + '; '.join(f'n={n},H={h}: {dict(sorted(c.items()))}' for (n, h), c in sorted(perm.items())))
v7 = sorted(set(int(h) for n in range(1, 10) for h, f in zip(data[n]['H'], data[n]['fa']) if f <= 2 and h % 8 == 7))
lines = lambda h: h in (15, 23) or any(h == 3 * 2 ** j - 1 or h == 5 * 2 ** j - 1 for j in range(3, 40))
check(all(lines(h) for h in v7), f'classes with fa <= 2 (n <= 9) have H = 7 mod 8 only on the lines: {v7}')
bad = 0; tot = 0; ncls = 0
for n in range(3, 9):
    for t in range(len(data[n]['H'])):
        f = int(data[n]['fa'][t])
        if 1 <= f <= 3:
            b, tt = min_orders_span1(n, data[n]['A'][t], f)
            bad += b; tot += tt; ncls += 1
check(bad == 0 and tot > 0, f'no minimum-feedback ordering has a back arc of span 1: {tot} minimum orderings of the {ncls} classes with 1 <= fa <= 3, n <= 8')

print('(T) THM-4749 (6), n = 8 row, from the full census (not from the strong inventory)')
N8 = Counter(int(h) for h in data[8]['H'])
row = [(N8[2 ** k - 1], N8[2 ** k + 1]) for k in range(4, 12)]
check(row == [(20, 6), (10, 23), (2, 25), (10, 40), (28, 30), (32, 35), (0, 0), (0, 0)], f'n = 8: {row}')
N9 = Counter(int(h) for h in data[9]['H'])
row9 = [(N9[2 ** k - 1], N9[2 ** k + 1]) for k in range(4, 12)]
check(row9 == [(30, 8), (14, 35), (4, 41), (22, 84), (91, 90), (136, 155), (193, 168), (149, 174)], f'n = 9 (repunit_rarity.out line 22): {row9}')

print('(O) census-derived numbers printed by the session')
fmin = {6: 15, 7: 25, 8: 45, 9: 75}
dens = {}
for m in (6, 7, 8, 9):
    bins = defaultdict(Counter)
    for h, s in zip(data[m]['H'], data[m]['st']):
        if not s:
            continue
        r = int(h) / fmin[m]
        key = '<2' if r < 2 else '2-4' if r < 4 else '4-8' if r < 8 else '>=8'
        bins[key][int(h) % 8] += 1
    dens[m] = {k: {rr: round(bins[k][rr] * 4 / sum(bins[k].values()), 2) for rr in (1, 3, 5, 7)} for k in bins}
print('     relative density of H mod 8 by octave H/f(m), strong classes: ' + str(dens))
check([dens[m]['2-4'][7] for m in (6, 7, 8, 9)] == [0.47, 0.53, 0.78, 0.89], 'second octave [2f, 4f): 7 mod 8 relative density 0.47, 0.53, 0.78, 0.89 at orders 6-9 (note section 5)')
check([dens[m]['4-8'][7] for m in (7, 8, 9)] == [1.13, 1.01, 0.98], 'octave [4f, 8f): 1.13, 1.01, 0.98 at orders 7-9 ("about 1")')
# n = 8 multiplicity residual vs least fa
mf = defaultdict(lambda: 99)
for h, f in zip(data[8]['H'], data[8]['fa']):
    mf[int(h)] = min(mf[int(h)], int(f))
hs = [h for h in range(25, 331, 2) if N8[h] > 0]
X = np.array([[1, math.log(h), math.log(h) ** 2, h % 3 == 0, h % 5 == 0] for h in hs], dtype=float)
y = np.log([N8[h] for h in hs])
cc, *_ = np.linalg.lstsq(X, y, rcond=None); rr = y - X @ cc
corr = float(np.corrcoef(rr, [mf[h] for h in hs])[0, 1])
byf = defaultdict(list)
for h, x in zip(hs, rr):
    byf[mf[h]].append(x)
mean1 = float(np.mean(byf[1]))
print(f'     n = 8: corr(residual, least fa) = {corr:+.3f}; mean residual by least fa ' + str({k: (round(float(np.mean(v)), 3), len(v)) for k, v in sorted(byf.items())}))
check(abs(corr + 0.194) < 0.0005 and abs(mean1 - 0.645) < 0.0005, 'n = 8 regression: corr = -0.194 and mean residual +0.645 for the two least-fa-1 values (repunit_rarity.out (O))')
# refuted_leads (c): share of 7 mod 8 by fa among strong classes
shares = {}
for m in (6, 7, 8, 9):
    tab = defaultdict(Counter)
    for h, s, f in zip(data[m]['H'], data[m]['st'], data[m]['fa']):
        if s:
            tab[int(f)][int(h) % 8] += 1
    shares[m] = {f: (round(tab[f][7] / sum(tab[f].values()), 3), sum(tab[f].values())) for f in sorted(tab)}
print('     share of 7 mod 8 by least fa (strong classes): ' + str(shares))
exp = {6: {1: (0.0, 1), 2: (0.462, 13), 3: (0.118, 17), 4: (0.0, 4)},
       7: {1: (0.0, 1), 2: (0.174, 23), 3: (0.173, 104), 4: (0.15, 147), 5: (0.313, 67), 6: (0.556, 9), 7: (0.0, 2)}}
check(shares[6] == exp[6] and shares[7] == exp[7] and shares[8][2] == (0.114, 35) and shares[8][5] == (0.251, 2007)
      and shares[9][3] == (0.332, 695) and shares[9][2] == (0.082, 49), 'refuted_leads.out (c) fa shares reproduced (orders 6-9)')
# class_residues: class counts by H mod 3 and labelled shares at n = 9
c3 = Counter(int(h) % 3 for h in data[9]['H'])
labs = Counter()
for h, a in zip(data[9]['H'], data[9]['aut']):
    labs[int(h) % 3] += math.factorial(9) // int(a)
tot9 = 2 ** 36
check((c3[0], c3[1], c3[2]) == (66320, 62129, 63087), f'n = 9 classes by H mod 3: {(c3[0], c3[1], c3[2])} (66,320 / 62,129 / 63,087)')
sh = [3 * labs[r] / tot9 for r in range(3)]
check(max(abs(x - 1) for x in sh) < 0.02, f'n = 9 labelled relative shares by H mod 3 within 2% of uniform: {[round(x, 4) for x in sh]}')
lab7 = []
for n in range(5, 10):
    w7 = sum(math.factorial(n) // int(a) for h, a in zip(data[n]['H'], data[n]['aut']) if int(h) % 8 == 7)
    lab7.append(round(4 * w7 / 2 ** (n * (n - 1) // 2), 4))
check(lab7 == [0.25, 0.7969, 0.8069, 0.9383, 0.993], f'labelled relative density of H = 7 mod 8 at n = 5..9: {lab7} (note: "0.25, 0.80, 0.81, 0.94, 0.993"; these are densities relative to 1/4, not shares)')
# order-7 strong spectrum
S7 = Counter(); sc7 = defaultdict(list)
for h, s, scv in zip(data[7]['H'], data[7]['st'], data[7]['sc']):
    if s:
        S7[int(h)] += 1; sc7[int(h)].append(int(scv))
holes = [h for h in range(47, 150, 2) if S7[h] == 0]
single = [h for h in range(47, 150, 2) if S7[h] == 1]
print(f'     order-7 strong spectrum in (45, 150): holes {holes}; values with one class {single} (self-converse: {[sc7[h][0] for h in single]})')
check(holes == [63, 107, 119, 149], 'order-7 holes in (45,150) are exactly 63, 107, 119, 149')
check(S7[63] == 0 and S7[127] == 1, 'S_7(63) = 0 and S_7(127) = 1')
check(all(S7[h] == 1 and sc7[h] == [1] for h in (61, 85, 87, 125, 127)), 'the thin values 61, 85, 87, 125, 127 are each attained by exactly one strong 7-class, and it is self-converse')
print(f'     all order-7 strong values attained by exactly one class: {sorted(h for h in S7 if S7[h] == 1)}')
# N(8,h) = 2 S_7(h) + S_8(h) for h > 45
S8 = Counter(int(h) for h, s in zip(data[8]['H'], data[8]['st']) if s)
check(all(N8[h] == 2 * S7[h] + S8[h] for h in range(47, 700, 2)), 'N(8, h) = 2 S_7(h) + S_8(h) for every odd h > 45 (note section 5)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
