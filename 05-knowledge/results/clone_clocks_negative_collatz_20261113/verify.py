"""Receipt for the clone-clock / negative-Collatz packet (klein, 2026-11-13).

Run from the repository root or from this directory:
    python3 -B 05-knowledge/results/clone_clocks_negative_collatz_20261113/verify.py
    python3 -O -B 05-knowledge/results/clone_clocks_negative_collatz_20261113/verify.py
Add --full to rescreen the 402 order-10 digit-1 clocks for ready vertices (about a minute).
Checks raise explicit exceptions (no bare asserts), so -O does not weaken them.
"""
import gzip, itertools, os, random, re, sys
from collections import Counter
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'scripts'))
from funnel import ready_count, build_nested, append_ready, hp_stats, H_only  # noqa: E402
from twomod import visit_table  # noqa: E402
from certify_twoclock import build, complete_table, pform, substitute, TT, C3, P_poly  # noqa: E402
from digits import digits  # noqa: E402
from tern_check import parse, blowup  # noqa: E402


def need(cond, msg):
    if not cond:
        raise RuntimeError('CHECK FAILED: ' + msg)


def T(x):  # accelerated 3x+1 map on Z (floor division handles negatives)
    return (3 * x + 1) // 2 if x % 2 else x // 2


def Tn(x, k):
    for _ in range(k):
        x = T(x)
    return x


out = []

# ---------------------------------------------------------------- 1. Collatz dictionary
need(all(T(2 * h + 1) == 3 * h + 2 == 2 * T(h) + 1 for h in range(-2001, 2002, 2)), 'C2 = T o D = D o T on odd H')
for H in range(-2001, 2002, 2):
    r = (H + 1) // 2
    for sig in (0, 1, 2):
        need((3 * H + 2 * sig + 1) // 2 == 3 * r + sig - 1, 'C2 digit sig is r -> 3r + sig - 1')
# unique C2 parent of every r >= 2 with digit sig in {0,1,2}; nonzero digits exactly at r not divisible by 3
for r in range(2, 5000):
    par = [(p, s) for s in (0, 1, 2) for p in [(r - s + 1) / 3] if p == int(p) and p >= 1]
    need(len(par) == 1, 'unique C2 parent')
    need((par[0][1] != 1) == (r % 3 != 0), 'nonzero digit iff r not divisible by 3')
need(((7 + 1) // 2, (21 + 1) // 2) == (4, 11) and 4 == 3 * 1 + 1 and 11 == 3 * 4 - 1, 'r(7)=4=3*1+1, r(21)=11=3*4-1')
cycles = {}
for start in (-1, -5, -17, 1):
    x, cyc = start, [start]
    while True:
        x = T(x)
        if x == start:
            break
        cyc.append(x)
    cycles[start] = (len(cyc), sum(1 for c in cyc if c % 2))
need(cycles == {-1: (1, 1), -5: (3, 2), -17: (11, 7), 1: (2, 1)}, 'cycle lengths / odd steps')
for d, (L, A) in cycles.items():
    for c in (1, 5, 60, 462, 1386):
        for t in range(0, 4):
            need(Tn(c * 2 ** (L * t) + d, L * t) == c * 3 ** (A * t) + d, 'shadowing exchange 2^L -> 3^A at anchor %d' % d)
need([Tn(462 * 8 ** t - 5, 3 * t) for t in range(1, 6)] == [462 * 9 ** t - 5 for t in range(1, 6)], '462 ray: 8 -> 9')
need([3 * (-5) + 2, 3 * (-13) + 2] == [-13, -37], 'anchors -5 -> -13 -> -37 under C2')
need(-37 in [-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34], '-37 lies on the -17 cycle')
need(Tn(-13, 4) == -7, '-13 is preperiodic into the -5 cycle')
# ascending runs: n = -1 mod 2^s  <=> first s-1 steps of T odd ;  n = +1 mod 2^s <=> same for (3x-1)/2
def Tm(x):
    return (3 * x - 1) // 2 if x % 2 else x // 2
for s in range(2, 6):
    for n in range(1, 4000, 2):
        runp, runm, x, y = 0, 0, n, n
        while x % 2 and runp < s:
            x = T(x); runp += 1
        while y % 2 and runm < s:
            y = Tm(y); runm += 1
        need((n % 2 ** s == 2 ** s - 1) == (runp >= s), '3x+1 ascents vs n = -1 mod 2^s')
        need((n % 2 ** s == 1) == (runm >= s), '3x-1 ascents vs n = +1 mod 2^s')
out.append('1. Collatz dictionary: C2 = T.D = D.T; r-digits 3r-1/3r/3r+1; shadowing at -1,-5,-17,+1; ascent classes mod 2^s: OK')

# ---------------------------------------------------------------- 2. the 462 chain (inherited F3) in the funnel compiler
for k in range(1, 7):
    ws = ['1', '0' + '1' * k + '1011']
    R, _ = ready_count(ws)
    need(R == 462 * 2 ** k - 5, '462 ray via the six-state compiler')
    if k <= 4:
        adj, n, root = build_nested(ws)
        A2, n2, z = append_ready(adj, n, root)
        need(hp_stats(A2, n2)[0] == R, '462 ray literal count')
out.append('2. 462*2^k-5 = words (1, 0 1^k 1011) + ready z: compiler k<=6, literal k<=4: OK')

# ---------------------------------------------------------------- 3. the Pi skeleton (+1 anchor) and the one-arc 3x-1 move
order = ['a', 'v', 'c', 'u', 'd0', 'd1', 'x0', 'f0', 'f1']
back = {('d0', 'a'), ('d1', 'a'), ('f0', 'c'), ('f1', 'c')}
adj, idx = build(order, back)
n = len(order); u, v = idx['u'], idx['v']
tab = complete_table(adj, n, u, v)
need(tab == {(1, 1): 361, (1, 2): 180, (2, 1): 360, (2, 2): 180, (3, 1): 120, (3, 2): 60}, 'Pi visit table')
need(pform(tab) == {(1, 1): 1, (3, 2): 60}, 'Pi P-form = 60 P_M(3) P_N(2) + H(M) H(N)')
mods = [TT(1), TT(2), TT(3), C3, TT(4)]
for M in mods:
    for N in mods[:4]:
        A, tot = substitute(adj, n, u, v, M, N)
        if tot > 15:
            continue
        need(H_only(A, tot) == 60 * P_poly(M, 3) * P_poly(N, 2) + H_only(M, len(M)) * H_only(N, len(N)), 'Pi literal')
others = [x for x in range(n) if x not in (u, v)]
need([order[x] for x in others if (adj[u] >> x & 1) != (adj[v] >> x & 1)] == ['c'], 'u and v differ only at c')
need(adj[v] >> u & 1 == 1, 'v -> u')
for j in range(1, 4):
    for k in range(2, 5):
        A, tot = substitute(adj, n, u, v, TT(j), TT(k))
        H = H_only(A, tot)
        need(H == 1 + 60 * 2 ** k * 3 ** j, 'Pi transitive count')
        cpos = others.index(idx['c']); w = len(others) + j + k - 1
        need(A[w] >> cpos & 1 == 1, 'last V clone beats c')
        B = list(A); B[w] &= ~(1 << cpos); B[cpos] |= 1 << w
        H2 = H_only(B, tot)
        need(H2 == (3 * H - 1) // 2 == 1 + 60 * 3 ** (j + 1) * 2 ** (k - 1) and -H2 == T(-H), 'one-arc 3x-1 step')
out.append('3. Pi: H = 60 P_M(3) P_N(2) + H(M)H(N) (complete table + literal, 20 module pairs); one arc flip = (3H-1)/2: OK')

# ---------------------------------------------------------------- 3b. the digit-3 two-clock Psi (Session-197 skeleton + 2 vertices)
psi = [974, 1020, 568, 1008, 993, 961, 644, 772, 68, 256]
tpsi = complete_table(psi, 10, 3, 9)
need(pform(tpsi) == {(1, 1): -3, (3, 1): 40, (3, 2): 160}, 'Psi P-form = 160 P3 P2 + 40 P3 H(N) - 3 H(M)H(N)')
for M in (TT(1), TT(2), C3):
    for N in (TT(1), TT(2), TT(3)):
        A, tot = substitute(psi, 10, 3, 9, M, N)
        if tot > 14:
            continue
        need(H_only(A, tot) == 160 * P_poly(M, 3) * P_poly(N, 2) + 40 * P_poly(M, 3) * H_only(N, len(N))
             - 3 * H_only(M, len(M)) * H_only(N, len(N)), 'Psi literal')
out.append('3b. Psi: H = 3^j(160*2^k+40) - 3, exact digit-3 step H -> 3H+6 for every k; complete table + literal: OK')

# ---------------------------------------------------------------- 3c. the refutation witnesses of 'c11 never -1'
wit = {'genuine': ([1998, 2044, 1592, 2032, 2017, 1985, 1156, 1796, 1092, 320, 512], 3, 10, {(1, 1): -1, (1, 2): -2, (3, 1): 60, (3, 2): 440}),
       'product60': ([2042, 2044, 2041, 1776, 1760, 1472, 1408, 1792, 1560, 96, 512], 7, 0, {(1, 1): -1, (1, 2): -1, (3, 1): 60, (3, 2): 60}),
       'product100': ([1530, 2044, 2041, 2032, 1504, 192, 1920, 1792, 1568, 1073, 32], 8, 2, {(1, 1): -1, (1, 2): -1, (3, 1): 100, (3, 2): 100})}
for name, (rows, uu, vv, pf) in wit.items():
    need(pform(complete_table(rows, 11, uu, vv)) == pf, 'witness P-form ' + name)
    for j in (1, 2):
        for k in (1, 2):
            A, tot = substitute(rows, 11, uu, vv, TT(j), TT(k))
            need(H_only(A, tot) == sum(c * m ** j * mp ** k for (m, mp), c in pf.items()), 'witness literal ' + name)
need(all(60 * 3 ** j * 2 ** k + 60 * 3 ** j - 2 ** k - 1 == (60 * 3 ** j - 1) * (2 ** k + 1) for j in range(1, 5) for k in range(1, 5)),
     'product form (60*3^j-1)(2^k+1)')
out.append('3c. refutation of c11 != -1: products (60*3^j-1)(2^k+1), (100*3^j-1)(2^k+1) and the genuine '
           '3^j(440*2^k+60)-2^(k+1)-1, all certified: OK')

# ---------------------------------------------------------------- 4. Session 197 skeleton with v = b
order0 = ['a', 'v', 'c', 'u', 'd0', 'd1', 'f0', 'f1']
adj0, idx0 = build(order0, back)
tab0 = complete_table(adj0, len(order0), idx0['u'], idx0['v'])
need(pform(tab0) == {(1, 1): 1, (1, 2): -1, (3, 2): 30}, 'S60 with v=b P-form')
for M in mods[:4]:
    A, tot = substitute(adj0, len(order0), idx0['u'], idx0['v'], M, TT(1))
    need(H_only(A, tot) == 60 * P_poly(M, 3) - H_only(M, len(M)), 'Session 197: 60 P_M(3) - H(M)')
out.append('4. Session-197 clock = 60 P_M(3) - H(M) for M in TT1..TT3, C3; with v=b: 30P3P2 - H(M)P2 + H(M)H(N): OK')

# ---------------------------------------------------------------- 5. excursion factorials and the odd-digit theorem
rng = random.Random(2026)
for trial in range(300):
    m = rng.randint(3, 8)
    A = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < 0.5: A[i] |= 1 << j
            else: A[j] |= 1 << i
    uu = rng.randrange(m)
    t = visit_table(A + [0], m + 1, uu, m, CU=m + 1, CV=1)  # dummy isolated v never used: count words with v absent
    a = Counter()
    for (i, ip), c in t.items():
        if ip == 0 and i >= 1:
            a[i] += c
    for i, c in a.items():
        need(c % factorial(i - 1) == 0, '(i-1)! divides a_i')
out.append('5. (i-1)! | a_i on 300 random rooted tournaments: OK')

# ---------------------------------------------------------------- 6. census data
def lines(name):
    p = os.path.join(HERE, 'data', name)
    with gzip.open(p, 'rt') as f:
        return f.read().splitlines()

tern = {n: [l for l in lines('tern_all_orders.txt.gz') if l.startswith('TERN %d ' % n)] for n in range(5, 10)}
need([len(tern[n]) for n in range(5, 10)] == [0, 4, 30, 206, 1309], 'complete ternary census counts')
for n in range(5, 10):
    for l in tern[n]:
        A_ = int(re.search(r' A=(\d+)', l).group(1)); C_ = int(re.search(r' C=(-?\d+)', l).group(1))
        need(A_ % 2 == 0 and C_ % 2 != 0, 'A even, digit odd')
def cval(l):
    return int(re.search(r' C=(-?\d+)', l).group(1))
dig1 = Counter(int(re.search(r' A=(\d+)', l).group(1)) for n in range(5, 10) for l in tern[n] if cval(l) == -1)
need(dig1 == Counter({60: 6, 100: 2, 156: 2, 180: 4}), 'digit-1 amplitudes through order 9: %s' % dig1)
gen = [l for l in lines('gen_all_orders.txt.gz') if l.startswith('GEN')]
need(Counter(int(l.split()[1]) for l in gen) == Counter({8: 10, 9: 132}), 'complete two-clock census counts')
allgen = gen + [l for l in lines('backscan_n10_b5.txt.gz') + lines('backscan_n11_b4.txt.gz') + lines('extscan_e1.txt.gz')
                + lines('extscan_e2.txt.gz') + lines('extone_from_certified.txt.gz') if l.startswith('GEN')]
pairs = Counter(); cm1 = []
for l in allgen:
    s = int(re.search(r'sigma=(-?\d+)', l).group(1)); s2 = int(re.search(r'sig2=(-?\d+)', l).group(1))
    A_ = int(re.search(r' A=(\d+)', l).group(1)); al = int(re.search(r' al=(\d+)', l).group(1))
    pairs[(s, s2)] += 1
    need((s, s2) != (1, 0), 'no C2-line two-clock')
    need((s, s2) != (1, 0), 'no C2-line two-clock (sigma, sig2) = (1, 0)')
    if 2 * s2 - s == -1:
        cm1.append((A_, s, al, s2))
    for j in range(1, 8):
        lead, const = 3 ** j * al - s2, 3 ** j * (A_ - 2 * al) - (s - 2 * s2)
        for (Tt, Bb) in ((45, 7), (135, 23), (405, 71), (1215, 215), (462, -5), (1386, -13), (4158, -37)):
            for sh in range(0, 6):
                need(not (lead == Tt * 2 ** sh and const == Bb), 'no two-clock contains a 45/462 chain ray')
need(not any(l.startswith('TWO') for l in lines('backscan_n10_b5.txt.gz') + lines('backscan_n11_b4.txt.gz')), 'no TWO lines')
pure = Counter(s for (s, s2) in pairs.elements() if s2 == 0)
need(set(pure) == {-11, -9, -5, -3, -1, 3}, 'pure v-sector digits observed: %s' % dict(pure))
need(sorted(cm1) == [(180, 3, 60, 1), (300, 3, 100, 1), (940, 5, 440, 2)], 'the three c11 = -1 lines (refutation witnesses)')
out.append('6. census: ternary 0/4/30/206/1309 (n=5..9, all A even, all digits odd); two-clock lines %d, '
           '(sigma,sig2)=(1,0) never; c11=-1 only in the 3 refutation witnesses; pure v-sectors with sigma in %s; '
           'no 45/462 ray: OK' % (len(allgen), dict(sorted(pure.items()))))

# ---------------------------------------------------------------- 7. readiness screen
cnt = 0
for n in (8, 9):
    for l in tern[n]:
        if cval(l) != -1:
            continue
        code = l.split()[2]
        H_, sig, tau = digits(parse(code, n), n)
        cnt += 1
        need(1 not in sig and 1 not in tau, 'no digit-1 vertex in a digit-1 clock (n<=9)')
out.append('7. readiness: none of the %d digit-1 clocks of order <= 9 has a C2-digit-1 vertex: OK' % cnt)
if '--full' in sys.argv:
    seen = set(); c10 = 0
    for l in lines('backscan_n10_b5.txt.gz'):
        if not l.startswith('TERN') or ' sigma=1 ' not in l:
            continue
        uu = int(re.search(r' u=(\d+)', l).group(1))
        bk = tuple(sorted(tuple(map(int, p.split('>'))) for p in re.search(r'back=([\d>,]+)', l).group(1).strip(',').split(',')))
        if (uu, bk) in seen:
            continue
        seen.add((uu, bk))
        A = [0] * 10
        for i in range(10):
            for j in range(i + 1, 10):
                if (j, i) in bk: A[j] |= 1 << i
                else: A[i] |= 1 << j
        H_, sig, tau = digits(A, 10)
        need(1 not in sig and 1 not in tau, 'no digit-1 vertex (n=10)')
        c10 += 1
    out.append('7b. readiness: none of the %d near-transitive order-10 digit-1 clocks has a digit-1 vertex: OK' % c10)

# ---------------------------------------------------------------- 8. rigorous certificates of every two-clock tuple
import subprocess, tempfile
import certify_all
with tempfile.TemporaryDirectory() as d:
    exe = os.path.join(d, 'certify')
    subprocess.run(['cc', '-O2', '-o', exe, os.path.join(HERE, 'c', 'certify.c')], check=True)
    certified, unit, seeded = certify_all.main(exe)
need(certified == 868 and len(unit) == 217, 'certified tuple counts')
need(all(f['alpha'] * 3 ** (j + 1) - f['sig2'] == 3 * (f['alpha'] * 3 ** j - f['sig2']) + 2 * f['sig2'] for f in unit for j in range(1, 6)),
     'seed multipliers follow c -> 3c + 2 sig2')
need(any(f['alpha'] == 30 and f['sig2'] == 1 for f in unit), 'the C2-line multiplier family 30*3^j - 1')
out.append('8. all 868 two-clock tuples certified by complete tables (rigorous caps) + literal counts; 217 unit-offset '
           'families give seeds; multipliers follow c -> 3c + 2 sig2: OK')

print('\n'.join(out))
print('ALL CHECKS PASSED')
