#!/usr/bin/env python3
"""Trailing-ones switches: the parity law for every odd source, the sibling form of the reset switch,
{7, 21} as the two binary-append maps generating F_21 = Aut(P_7), CRT independence of the Mersenne
mirror clocks, and the collision (birthday) statistics of the word map at -1.

Session opus-2026-10-06-S18.  Note: 05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md
Complements mac-mini's THM-4555 / THM-4556 / HYP-9213 / HYP-9214 (same owner prompt); reproduces their key
numbers independently where used (sections A, E).

Sections:
 A. THM-4555 (vi) census reproduced: 37 of the 60 odd a in [3,121]; least lag odd; lag sets are pairs {D, D+1}.
 B. Parity law for every odd n < 2*10^4 (run r >= 2): predicted good parity (-1)^(r-D+1) = t (mod 4);
    L+(n) = L(n) u {0 if reset(n) >= 3} is a union of pairs {D, D+1}, D good, except an unpaired good top lag r-1;
    least lag odd for every reset-2 source; also in THM-4555's wider convention (runs >= 1, lags D <= r).
 C. Ahmed's identity U^r(2m+1) = 2^(e(m)) U^r(m) + 1, and its sibling case U^r(2m+1) = 4 U^r(m) + 1 iff reset(2m+1) >= 3.
 D. {7, 21}: A_1(x) = 2x+1 (append 1) and R(x) = 4x+1 (append 01) generate F_21 = Aut(P_7) mod 7;
    mod 2^k - 1 they generate Z/p x| <2> (order kp), equal to Aut(P_p) only at k = 3;
    (4^k-1)/3 = k(2^k-1) iff 2^k + 1 = 3k iff k in {1, 3}; mod 63: R^3 = x + 21, A_1^6 = id.
 E. Equal-time switches = level sets of the odd-step time sigma; classes N(A) for A <= 4000, windows, least-squares
    fit; least-lag merge values (odd a); collisions (equal totals); template totals.
 F. CRT independence: forward switching vs backward-minimality (depth 30) of 2^a - 1, odd a <= 401.
 G. Word map at -1: reduced words of total A, distinct values, colliding pairs, sporadic pairs, Renyi-2 entropy.
 H. Exact 2-adic certified switching density (template total <= K), K <= 18; the a = 95 (mod 128) family
    (words (2,6,1) and (4,1,1,3), value 125/256) checked at large a.

Reproduce: python3 04-computation/experiments/mersenne_switch_parity_f21_compression_20261006.py  (a few minutes)
"""
import math
from collections import Counter, defaultdict
from fractions import Fraction as Fr

FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def U(x):
    y = 3 * x + 1
    return y >> ((y & -y).bit_length() - 1)


def Uk(x, k):
    for _ in range(k):
        x = U(x)
    return x


def orbit(n):
    out = [n]
    while n != 1:
        n = U(n)
        out.append(n)
    return out


def sigma(n):
    s = 0
    while n != 1:
        n = U(n)
        s += 1
    return s


def run_t(n):
    """n = 2^(r+1) t - 1 with t odd -> (r, t)"""
    y = n + 1
    e = (y & -y).bit_length() - 1
    return e - 1, y >> e


def reset_exp(n):
    r, t = run_t(n)
    y = 3 * Uk(n, r) + 1
    return (y & -y).bit_length() - 1


def meet_index(oa, ob):
    pb = {v: i for i, v in enumerate(ob)}
    for i, v in enumerate(oa):
        if pb.get(v) == i:
            return i
    return None


def word(x, steps):
    out = []
    for _ in range(steps):
        y = 3 * x + 1
        k = (y & -y).bit_length() - 1
        out.append(k)
        x = y >> k
    return out


def fval(u):
    x = Fr(-1)
    for c in u:
        x = (3 * x + 1) / 2 ** c
    return x


# ---------------------------------------------------------------- A
print('== A. THM-4555 (vi) census reproduced')
orbM = {b: orbit(2 ** b - 1) for b in range(2, 122)}
rowsA = []
for a in range(3, 122, 2):
    lags = [D for D in range(1, a - 1) if meet_index(orbM[a], orbM[a - D]) is not None]
    rowsA.append((a, lags))
swA = [a for a, l in rowsA if l]
print(f'   odd a in [3,121]: {len(rowsA)}; with an equal-time Mersenne partner: {len(swA)}', check('37/60', len(swA) == 37 and len(rowsA) == 60))
print('   least lag odd for all:', check('least odd', all(min(l) % 2 == 1 for a, l in rowsA if l)))
pairs_ok = all(((D % 2 == 0 and D - 1 in l) or (D % 2 == 1 and (a - D - 1 < 2 or D + 1 in l))) for a, l in rowsA for D in l)
print('   lag sets are unions of pairs {D, D+1}, D odd:', check('pairs', pairs_ok))
print('   non-switching odd a:', [a for a, l in rowsA if not l])

# ---------------------------------------------------------------- B
print('\n== B. Parity law for every odd source n < 2*10^4 with run r >= 2')
cache = {}


def orb(n):
    o = cache.get(n)
    if o is None:
        o = cache[n] = orbit(n)
    return o


pred_fail = pair_viol = 0
n_src = n_r2 = n_r2_sw = 0
least_odd = True
for n in range(3, 20000, 2):
    r, t = run_t(n)
    if r < 2:
        continue
    n_src += 1
    on = orb(n)
    L = set()
    for D in range(1, r):
        m = (n + 1) // 2 ** D - 1
        if m > 1 and meet_index(on, orb(m)) is not None:
            L.add(D)
    for D in range(0, r):
        m = (n + 1) // 2 ** D - 1
        if m > 1 and (reset_exp(m) >= 3) != ((((-1) ** (r - D + 1)) - t) % 4 == 0):
            pred_fail += 1
    for D in L:
        good = reset_exp((n + 1) // 2 ** D - 1) >= 3
        P = D + 1 if good else D - 1
        if 1 <= P < r and (n + 1) // 2 ** P - 1 > 1 and P not in L:
            pair_viol += 1
    if reset_exp(n) == 2:
        n_r2 += 1
        if L:
            n_r2_sw += 1
            least_odd &= min(L) % 2 == 1
print(f'   sources: {n_src}; good-parity prediction failures: {pred_fail}', check('B pred', pred_fail == 0),
      f'; pair violations (closure inside 1..r-1): {pair_viol}', check('B pairs', pair_viol == 0))
print(f'   reset-2 sources: {n_r2}, with a trailing-ones equal-time partner: {n_r2_sw}; least lag odd for all: {least_odd}', check('B odd', least_odd))
# corrected form: L+ = L u {0 if reset >= 3}; pairs {D, D+1} with D good, except an unpaired good top lag r-1
top_exc = bad_form = 0
for n in range(3, 20000, 2):
    r, t = run_t(n)
    if r < 2:
        continue
    on = orb(n)
    L = {D for D in range(1, r) if (n + 1) // 2 ** D - 1 > 1 and meet_index(on, orb((n + 1) // 2 ** D - 1)) is not None}
    Lp = set(L) | ({0} if reset_exp(n) >= 3 else set())
    for D in Lp:
        good = D == 0 or reset_exp((n + 1) // 2 ** D - 1) >= 3
        if good:
            if D + 1 <= r - 1 and (n + 1) // 2 ** (D + 1) - 1 > 1:
                bad_form += (D + 1) not in Lp
            elif D == r - 1:
                top_exc += 1
        else:
            bad_form += (D - 1) not in Lp
print(f'   corrected form (L+ = L u {{0 if reset >= 3}}): violations {bad_form}', check('B form', bad_form == 0),
      f'; unpaired good top lags r-1: {top_exc} sources')
# THM-4555's wider convention: run r >= 1, lags 1 <= D <= r, partner > 1, reset-2 sources n < 2*10^4
wide_src = wide_sw = 0
wide_odd = True
for n in range(3, 20000, 2):
    r, t = run_t(n)
    if r < 1 or reset_exp(n) != 2:
        continue
    wide_src += 1
    on = orb(n)
    Lw = [D for D in range(1, r + 1) if (n + 1) // 2 ** D - 1 > 1 and meet_index(on, orb((n + 1) // 2 ** D - 1)) is not None]
    if Lw:
        wide_sw += 1
        wide_odd &= min(Lw) % 2 == 1
print(f'   wider convention (runs >= 1, D <= r): reset-2 sources {wide_src}, with a partner {wide_sw}; least lag odd for all: {wide_odd}',
      check('B wide', wide_odd))

# ---------------------------------------------------------------- C
print('\n== C. Sibling form of the reset switch')
okC, cntC, iffC = True, 0, True
for n in range(3, 200001, 2):
    r, t = run_t(n)
    if r < 1:
        continue
    m = (n - 1) // 2
    if m < 1:
        continue
    sib = Uk(n, r) == 4 * Uk(m, r) + 1
    big = reset_exp(n) >= 3
    iffC &= sib == big
    if big:
        cntC += 1
        okC &= U(Uk(n, r)) == U(Uk(m, r))
print(f'   odd n < 2*10^5, run >= 1: U^r(n) = 4 U^r(m) + 1 (m = (n-1)/2) iff reset(n) >= 3: {check("C iff", iffC)}; '
      f'merge at r+1 in all {cntC} reset>=3 cases: {check("C merge", okC)}')
print('   example: U^3(15) =', Uk(15, 3), '= 4*U^3(7) + 1 =', 4 * Uk(7, 3) + 1)
okA = True
for n in range(3, 200001, 2):
    r, t = run_t(n)
    m = (n - 1) // 2
    if r < 1 or m < 1:
        continue
    okA &= Uk(n, r) == 2 ** reset_exp(m) * Uk(m, r) + 1 if run_t(m)[0] == r - 1 else True
print('   Ahmed 2016 Thm 2.1, general form U^r(2m+1) = 2^(e(m)) U^r(m) + 1 (odd n < 2*10^5, run >= 1):', check('C ahmed', okA))

# ---------------------------------------------------------------- D
print('\n== D. {7, 21}: the two append maps')


def gen_group(p, gens):
    comp = lambda f, g: ((f[0] * g[0]) % p, (f[0] * g[1] + f[1]) % p)
    G, frontier = {(1, 0)}, [(1, 0)]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                c = comp(g, h)
                if c not in G:
                    G.add(c)
                    nxt.append(c)
        frontier = nxt
    return G


for k in (3, 5, 7):
    p = 2 ** k - 1
    G = gen_group(p, [(2, 1), (4, 1)])
    QR = {x * x % p for x in range(1, p)}
    autP = p * len(QR)
    print(f'   p = 2^{k} - 1 = {p}: |<A_1, R>| = {len(G)} (= kp = {k * p}); inside Aut(P_p): {all(a in QR for a, b in G)}; |Aut(P_p)| = {autP}; equal: {len(G) == autP}')
    check(f'D gen {k}', len(G) == k * p and all(a in QR for a, b in G) and ((len(G) == autP) == (k == 3)))
print('   k <= 300 with 2^(k-1) - 1 = k (<2> = QR_p):', [k for k in range(2, 301) if 2 ** (k - 1) - 1 == k])
print('   k <= 300 with (4^k-1)/3 = k(2^k-1), i.e. 2^k + 1 = 3k:', [k for k in range(1, 301) if 2 ** k + 1 == 3 * k],
      check('D 2k+1', [k for k in range(1, 301) if 2 ** k + 1 == 3 * k] == [1, 3]))
print('   A_1^3(n) - n = 7(n+1), R^3(n) - n = 21(3n+1):', check('D offs', all((8 * n + 7) - n == 7 * (n + 1) and (64 * n + 21) - n == 21 * (3 * n + 1) for n in range(-99, 99))))
print('   mod 63: A_1^6 = id, R^3 = x + 21:', check('D 63', all((64 * x + 63) % 63 == x and (64 * x + 21) % 63 == (x + 21) % 63 for x in range(63))))
print('   fixed points: A_1 fixes -1 (limit of 2^a - 1), R fixes -1/3 (limit of (4^k-1)/3); mod 7: A_1 fixes 6, R fixes 2')

# ---------------------------------------------------------------- E
print('\n== E. Equal-time classes of Mersenne numbers (sigma level sets)')
sig = {a: sigma(2 ** a - 1) for a in range(2, 2001)}
first = {}
for a in range(2, 2001):
    first.setdefault(sig[a], a)
roots = sorted(first.values())
print('   classes among exponents 2..A:', {A: sum(1 for r_ in roots if r_ <= A) for A in (100, 400, 1200, 2000)},
      ' (THM-4556 (vi) counts 23, 37, 58 for A = 100, 400, 1200 including a = 1)')
check('E counts', [sum(1 for r_ in roots if r_ <= A) for A in (100, 400, 1200)] == [22, 36, 57])
same = all((meet_index(orbM[a], orbM[b]) is not None) == (sig[a] == sig[b]) for a in range(3, 80) for b in range(2, a))
print('   equal-time meeting <=> equal odd-step time (a, b <= 80):', check('E iff', same))
tot = []
for a in range(3, 1201, 2):
    if first[sig[a]] == a:
        continue
    D = min(a - b for b in range(2, a) if sig[b] == sig[a])
    x, y, i = 2 ** a - 1, 2 ** (a - D) - 1, 0
    while x != y:
        x, y, i = U(x), U(y), i + 1
    ell = i - (a - 1)
    tot.append(sum(word(2 * 3 ** (a - 1) - 1, ell)))
tot.sort()
print(f'   switching odd a <= 1200: {len(tot)}; template total at the least-lag merge: median {tot[len(tot) // 2]}, '
      f'share <= 20: {sum(1 for v in tot if v <= 20) / 600:.3f} of all odd a (certified density 0.1199 at K = 20)')
# least-lag merge values and the collision (equal totals: A_n - A_m = D) for odd a <= 2000
minv, mina, coll_ok, n_sw2000 = None, None, True, 0
for a in range(3, 2001, 2):
    if first[sig[a]] == a:
        continue
    n_sw2000 += 1
    D = min(a - b for b in range(2, a) if sig[b] == sig[a])
    x, y, An, Am = 2 ** a - 1, 2 ** (a - D) - 1, 0, 0
    while x != y:
        yx, yy = 3 * x + 1, 3 * y + 1
        kx, ky = (yx & -yx).bit_length() - 1, (yy & -yy).bit_length() - 1
        x, y, An, Am = yx >> kx, yy >> ky, An + kx, Am + ky
    coll_ok &= An - Am == D
    if minv is None or x < minv:
        minv, mina = x, a
evens = {a: None for a in (4, 6, 8)}
for a in evens:
    x, y = 2 ** a - 1, 2 ** (a - 1) - 1
    while x != y:
        x, y = U(x), U(y)
    evens[a] = x
print(f'   odd a <= 2000 with a switch: {n_sw2000}; least-lag merge value >= {minv} (attained at a = {mina})',
      check('E 911', minv == 911 and mina == 107), f'; every least-lag switch is a collision (A_n - A_m = D): {check("E coll", coll_ok)}')
print('   even a merge with a - 1 (reset switch) at values:', evens)
sig4 = dict(sig)
for a in range(2001, 4001):
    sig4[a] = sigma(2 ** a - 1)
first4 = {}
for a in range(2, 4001):
    first4.setdefault(sig4[a], a)
roots4 = sorted(first4.values())
NA = {A: sum(1 for r_ in roots4 if r_ <= A) for A in (100, 200, 400, 800, 1000, 1600, 2000, 3000, 4000)}
print('   classes N(A) among exponents 2..A:', NA, ' N(2A)/N(A):', {A: round(NA[2 * A] / NA[A], 3) for A in (100, 200, 400, 800, 1000, 2000)})
xs = [math.log(A) for A in range(100, 4001, 50)]
ys = [math.log(sum(1 for r_ in roots4 if r_ <= A)) for A in range(100, 4001, 50)]
mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
slope = sum((u - mx) * (v - my) for u, v in zip(xs, ys)) / sum((u - mx) ** 2 for u in xs)
print(f'   least-squares log-log fit of N(A), A in [100, 4000]: exponent {slope:.3f} (HYP-9213: about 0.37)')
win = []
for lo in range(1, 4001, 500):
    w = [a for a in range(max(lo, 3), lo + 500, 2) if a <= 4000]
    win.append(round(sum(1 for a in w if first4[sig4[a]] != a) / len(w), 3))
print('   share of odd a switching, windows of 500 consecutive exponents (250 odd each), up to 4000:', win)

# ---------------------------------------------------------------- F
print('\n== F. CRT independence of the mirror clocks (forward switching vs backward-minimality, depth 30)')


def has_smaller(n, depth=30, budget=3 * 10 ** 6):
    N1 = n + 1
    stack, nodes = [(n, depth)], 0
    while stack:
        x, h = stack.pop()
        nodes += 1
        if nodes > budget:
            return None
        r_ = x % 3
        if r_ == 0 or h == 0:
            continue
        k = 2 if r_ == 1 else 1
        while True:
            xn = ((x << k) - 1) // 3
            if xn < n:
                return True
            if (xn + 1) * 2 ** (h - 1) >= 3 ** (h - 1) * N1:
                break
            stack.append((xn, h - 1))
            k += 2
    return False


hs = {a: has_smaller(2 ** a - 1) for a in range(3, 402, 2)}
print('   undecided (node budget exhausted):', [a for a, v in hs.items() if v is None], check('F budget', all(v is not None for v in hs.values())))
bmin = {a: (v is False) for a, v in hs.items()}
print('   backward-minimal to depth 30, odd 3 <= a <= 119:', sum(1 for a in range(3, 120, 2) if bmin[a]), 'of 59 (THM-4554 (vi): 35)',
      check('F 35', sum(1 for a in range(3, 120, 2) if bmin[a]) == 35))
for lo, hi in ((3, 121), (3, 401), (123, 401)):
    tab = Counter((first[sig[a]] != a, bmin[a]) for a in range(lo, hi + 1, 2))
    n_ = sum(tab.values())
    ps = (tab[(True, True)] + tab[(True, False)]) / n_
    pb = (tab[(True, True)] + tab[(False, True)]) / n_
    print(f'   odd a in [{lo},{hi}]: P(switch) = {ps:.3f}, P(bmin) = {pb:.3f}, P(both) = {tab[(True, True)] / n_:.3f}, product = {ps * pb:.3f}')
oroots = [r_ for r_ in roots if r_ % 2 and 3 <= r_ <= 401]
print(f'   class roots (odd, <= 401) that are backward-minimal to depth 30: {sum(1 for r_ in oroots if bmin[r_])} of {len(oroots)}')

# ---------------------------------------------------------------- G
print('\n== G. The word map at -1 (reduced words, first letter >= 2)')


def comps(total, first_min=2):
    def rec(rem):
        if rem == 0:
            yield ()
            return
        for c in range(1, rem + 1):
            for t_ in rec(rem - c):
                yield (c,) + t_
    for c0 in range(first_min, total + 1):
        for t_ in rec(total - c0):
            yield (c0,) + t_


def numer(u):
    X, e = -1, 0
    for c in u:
        X, e = 3 * X + (1 << e), e + c
    return X // 2


def canon(u):
    u = list(u)
    while len(u) >= 2 and u[0] == 2:
        u = [u[1] + 2] + u[2:]
    return tuple(u)


print('    A |       M |       V |  pairs | sporadic | V/M   | H2 - log2 M')
prev = None
for A in range(8, 19):
    cnt, cls = Counter(), defaultdict(set)
    for u in comps(A):
        N = numer(u)
        cnt[N] += 1
        cls[N].add(canon(u))
    M = sum(cnt.values())
    V = len(cnt)
    pairs = sum(m * (m - 1) // 2 for m in cnt.values())
    spor = sum(len(s) * (len(s) - 1) // 2 for s in cls.values())
    H2 = -math.log2(sum((m / M) ** 2 for m in cnt.values()))
    print(f'   {A:2d} | {M:7d} | {V:7d} | {pairs:6d} | {spor:8d} | {V / M:.3f} | {H2 - math.log2(M):+.3f}')
check('G M', M == 2 ** (18 - 2))

# ---------------------------------------------------------------- H
print('\n== H. Exact certified 2-adic switching density and the a = 95 (mod 128) family')


def word_mod(x, M):
    w, cur, bits = [], 0, M
    x %= (1 << M)
    while bits > 0:
        if x & 1:
            if cur:
                w.append(cur)
            cur = 1
            x = (3 * x + 1) >> 1
        else:
            cur += 1
            x >>= 1
        bits -= 1
        x &= (1 << bits) - 1 if bits > 0 else 0
    return w


def prefix_vals(w):
    out, X, e = {}, -1, 0
    for i, c in enumerate(w):
        X, e = 3 * X + (1 << e), e + c
        out[(e, i + 1)] = X
    return out


dens = []
for K in range(9, 19):
    Mb = K + 1
    mod = 1 << Mb
    inv3 = pow(3, -1, mod)
    good = total = 0
    for Xr in range(1, mod, 8):
        total += 1
        pv = prefix_vals(word_mod(2 * Xr - 1, Mb))
        hit = False
        for D in range(1, K, 2):
            qv = prefix_vals(word_mod(2 * Xr * pow(inv3, D, mod) - 1, Mb))
            if any(qv.get((e, L + D)) == v for (e, L), v in pv.items()):
                hit = True
                break
        good += hit
    dens.append(round(good / total, 4))
print('   certified share of odd exponents, K = 9..18:', dens, check('H dens', dens[:3] == [0.0156, 0.0312, 0.0469]))
u95, up95 = word(2 * 3 ** 94 - 1, 3), word(2 * 3 ** 93 - 1, 4)
print('   a = 95: post-run words', u95, 'and (partner 2^94 - 1)', up95, '; f(-1) =', fval(u95), fval(up95), '; normal forms', canon(u95), canon(up95))
check('H words', u95 == [2, 6, 1] and up95 == [4, 1, 1, 3] and fval(u95) == fval(up95) == Fr(125, 256))
okH = True
for a in (95 + 128 * t for t in range(0, 8)):
    x = Uk(2 * 3 ** (a - 1) - 1, 3)
    y = Uk(2 * 3 ** (a - 2) - 1, 4)
    okH &= x == y
print('   a = 95 + 128t, t = 0..7 (a up to 991): 2^a - 1 and 2^(a-1) - 1 meet at time a + 2:', check('H family', okH))

# ---------------------------------------------------------------- I
print('\n== I. openai/math spot checks')
BARK = {7: '+++--+-', 11: '+++---+--+-', 13: '+++++--++-+-+'}
for n_, sq in BARK.items():
    xs_ = [1 if c == '+' else -1 for c in sq]
    acf = [sum(xs_[i] * xs_[i + k] for i in range(n_ - k)) for k in range(1, n_)]
    pl = {i for i, c in enumerate(sq) if c == '+'}
    mi = {i for i, c in enumerate(sq) if c == '-'}
    QRn = {i * i % n_ for i in range(1, n_)}
    trans = lambda S, T: any({(t_ + c) % n_ for t_ in T} == S for c in range(n_))
    small = mi if len(mi) < len(pl) else pl
    dif = Counter((x_ - y_) % n_ for x_ in small for y_ in small if x_ != y_)
    print(f'   Barker {n_}: |C_k| <= 1: {all(abs(c) <= 1 for c in acf)}; minus set translate of QR: {trans(mi, QRn)}; '
          f'plus set translate of QR: {trans(pl, QRn)}; smaller set a ({n_},{len(small)},{set(dif.values())}) difference set')
check('I barker7', trans(mi, {1, 2, 4}) if False else any({(t_ + c) % 7 for t_ in (1, 2, 4)} == {3, 4, 6} for c in range(7)))
Avals = [0.3, 1.7, 2.2, 0.05, 1.1]
lhs = math.exp(-sum((a_ / sum(Avals)) * math.log(a_ / sum(Avals)) for a_ in Avals)) * math.prod(a_ ** (a_ / sum(Avals)) for a_ in Avals)
print(f'   Gibbs/Kraft identity max_q e^H(q) prod A_i^q_i = sum A_i at q_i = A_i/sum A: {lhs:.12f} vs {sum(Avals):.12f}', check('I gibbs', abs(lhs - sum(Avals)) < 1e-12))

print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')
