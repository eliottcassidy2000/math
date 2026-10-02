#!/usr/bin/env python3
"""procgen_tcpc_20261001_run.py -- TCPC lane runner (collatz-procgen-20260922, 2026-10-01).

Re-verifies every FINITE-EXACT / VERIFIED claim of
05-knowledge/results/procgen_tcpc_20261001_tournament_clock_prime_collatz.md
and ends with 'ALL CHECKS PASSED'.  Prints to stdout only.

Flags:
  --n26       also run the N = 26 subset-DP checks (D_13 on T and on T^op, plus the interval control),
              about 13 minutes and 257 MB.
  --long      also run the r = 11 rotational cross-set census, the N = 13 structured searches, F(P_23),
              the m = 23 circulant scan and D_15 (N = 30) by the inclusion-exclusion engines.
  --n34       also run D_17 (N = 34) by karp3 on T and karp4 on T^op (about 25 minutes, 1 MB).
  --n38       also run D_19 (N = 38) by karp4 (several hours).
  --only-n38  run only the D_19 check.

Run from anywhere:  python3 04-computation/experiments/procgen_tcpc_20261001_run.py
"""
import itertools
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_tcpc_20261001_lib as L  # noqa: E402

T0 = time.time()
NCHECK = 0


def check(cond, msg):
    global NCHECK
    NCHECK += 1
    if not cond:
        raise RuntimeError('CHECK FAILED: ' + msg)


def section(title):
    print('\n== %s  [t=%.1fs]' % (title, time.time() - T0), flush=True)


LEAN_SRC = os.path.join(HERE, 'procgen_tcpc_20261001_par26.c')
LEAN_BIN = os.path.join(L.SCRATCH, 'procgen_tcpc_20261001_par26')


def ensure_lean():
    os.makedirs(L.SCRATCH, exist_ok=True)
    if (not os.path.exists(LEAN_BIN)) or os.path.getmtime(LEAN_BIN) < os.path.getmtime(LEAN_SRC):
        subprocess.run(['clang', '-O3', '-o', LEAN_BIN, LEAN_SRC], check=True)


def lean_parity(A, arcs):
    ensure_lean()
    n = len(A)
    inp = '%d\n' % n + '\n'.join(' '.join(map(str, r)) for r in A) + '\n%d\n' % len(arcs) + \
        '\n'.join('%d %d' % e for e in arcs) + '\n'
    r = subprocess.run([LEAN_BIN], input=inp, capture_output=True, text=True, check=True)
    H, out = None, {}
    for ln in r.stdout.split('\n'):
        p = ln.split()
        if not p:
            continue
        if p[0] == 'H':
            H = int(p[1])
        else:
            out[(int(p[1]), int(p[2]))] = int(p[3])
    return H, out


def rooted(A):
    L.ensure_hp_binary()
    n = len(A)
    inp = '%d\n' % n + '\n'.join(' '.join(map(str, r)) for r in A) + '\n'
    r = subprocess.run([L.HP_BIN, 'rooted'], input=inp, capture_output=True, text=True, check=True)
    return int(r.stdout.split()[1])


KARP3_SRC = os.path.join(HERE, 'procgen_tcpc_20261001_karp3.c')
KARP3_BIN = os.path.join(L.SCRATCH, 'procgen_tcpc_20261001_karp3')
KARP4_SRC = os.path.join(HERE, 'procgen_tcpc_20261001_karp4.c')
KARP4_BIN = os.path.join(L.SCRATCH, 'procgen_tcpc_20261001_karp4')


def ensure_karp():
    os.makedirs(L.SCRATCH, exist_ok=True)
    for src, binp, flags in [(KARP3_SRC, KARP3_BIN, ['-O3']), (KARP4_SRC, KARP4_BIN, ['-O3', '-march=armv8-a+crypto'])]:
        if (not os.path.exists(binp)) or os.path.getmtime(binp) < os.path.getmtime(src):
            subprocess.run(['clang'] + flags + ['-o', binp, src], check=True)


def karp_run(binp, A, reps):
    """inclusion-exclusion parity engine (karp3 or karp4): returns (H mod 2, {arc: c mod 2}, #subset reps)."""
    ensure_karp()
    N = len(A)
    inp = '%d\n' % N + '\n'.join(' '.join(map(str, row)) for row in A) + '\n%d\n' % len(reps) + \
        '\n'.join('%d %d' % e for e in reps) + '\n'
    out = subprocess.run([binp], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    H, C, R = None, {}, None
    for ln in out:
        p = ln.split()
        if not p:
            continue
        if p[0] == 'H':
            H = int(p[1])
        elif p[0] == 'C':
            C[(int(p[1]), int(p[2]))] = int(p[3])
        elif p[0] == 'REPS':
            R = int(p[1])
    return H, C, R


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def half_paths(G_elems, add, neg, zero, C):
    """F(D): paths 0 -> w1 -> ... -> wk (k=(n-1)/2) in Clk(G,C) picking one of each {x,-x}."""
    mult = {}
    for c in C:
        mult[c] = mult.get(c, 0) + 1
    k = (len(G_elems) - 1) // 2
    cnt = 0

    def dfs(cur, used, depth, w):
        nonlocal cnt
        if depth == k:
            cnt += w
            return
        for c, mu in mult.items():
            if c == zero:
                continue
            nx = add(cur, c)
            if nx == zero or nx in used or neg(nx) in used:
                continue
            used.add(nx)
            dfs(nx, used, depth + 1, w * mu)
            used.discard(nx)
    dfs(zero, set(), 0, 1)
    return cnt


def zm(m):
    return list(range(m)), (lambda x, y: (x + y) % m), (lambda x: (-x) % m), 0


def multiplier_group(m, C):
    mult = sorted(C)
    return [u for u in range(1, m) if math.gcd(u, m) == 1 and sorted((u * c) % m for c in C) == mult]


def T_sheet(r, S, X):
    """two-sheet clock: a_j = 2j, b_j = 2j+1; a_j->a_k iff k-j in S; b_j->b_k iff j-k in S;
    a_j -> b_k iff k-j in X else b_k -> a_j."""
    N = 2 * r
    A = [[0] * N for _ in range(N)]
    for j in range(r):
        for k in range(r):
            if j != k:
                if (k - j) % r in S:
                    A[2 * j][2 * k] = 1
                if (j - k) % r in S:
                    A[2 * j + 1][2 * k + 1] = 1
            if (k - j) % r in X:
                A[2 * j][2 * k + 1] = 1
            else:
                A[2 * k + 1][2 * j] = 1
    return A


def D_r(r):
    m = (r - 1) // 2
    S = set(range(1, m + 1))
    Y = set()
    for d in range(1, m // 2 + 1):
        Y |= {d % r, (-d) % r}
    if m % 2 == 1:
        Y |= {0}
    return T_sheet(r, S, Y)


def arc_orbit_reps(r, A):
    N = len(A)
    reps, seen = [], set()
    for u in range(N):
        for v in range(N):
            if u != v and A[u][v]:
                key = (u % 2, v % 2, (v // 2 - u // 2) % r)
                if key not in seen:
                    seen.add(key)
                    reps.append((u, v))
    return reps


def all_odd_full(A):
    res = L.hp_c(A, 'par')
    return all(v == 1 for v in res['arc'].values()), res


def paley_minus_0(q):
    QR = set(L.qr_set(q))
    V = list(range(1, q))
    return [[1 if (x != y and ((y - x) % q) in QR) else 0 for y in V] for x in V]


def paley_sheet_data(q):
    g = L.primitive_root(q)
    r = (q - 1) // 2
    QR = set(L.qr_set(q))
    chi = lambda x: 1 if x % q in QR else -1
    SA = {d for d in range(1, r) if chi(pow(g, 2 * d, q) - 1) == 1}
    X = {e for e in range(r) if chi(pow(g, 2 * e + 1, q) - 1) == 1}
    return SA, X


def H_memo(A):
    """independent engine: memoized recursion."""
    from functools import lru_cache
    n = len(A)
    out = [[j for j in range(n) if A[i][j]] for i in range(n)]

    @lru_cache(maxsize=None)
    def ext(rem, v):
        if rem == 0:
            return 1
        s = 0
        for w in out[v]:
            if rem >> w & 1:
                s += A[v][w] * ext(rem & ~(1 << w), w)
        return s
    full = (1 << n) - 1
    val = sum(ext(full & ~(1 << v), v) for v in range(n))
    ext.cache_clear()
    return val


def main():
    args = set(sys.argv[1:])
    random.seed(20261001)
    print('procgen_tcpc_20261001_run.py -- TCPC lane runner')
    L.ensure_hp_binary()
    ensure_lean()
    ensure_karp()
    if '--only-n38' in args:
        section('D_19 (N = 38) by the inclusion-exclusion engine karp4')
        A = D_r(19)
        check(L.is_tournament(A), 'D_19 tournament')
        reps = arc_orbit_reps(19, A)
        H4, C4, R4 = karp_run(KARP4_BIN, A, reps)
        check(H4 == 1 and len(reps) == 37 and all(C4.values()), 'D_19 karp4')
        print('  D_19 (N = 38): H odd, all 37 arc orbits odd (%d subset reps)' % R4)
        print('\n%d checks, %.1f s' % (NCHECK, time.time() - T0))
        print('ALL CHECKS PASSED')
        return

    # ------------------------------------------------------------------
    section('A. engines: C exact / C parity / Python DP / brute force / memo agree')
    for trial in range(60):
        n = random.randint(2, 7)
        A = [[(random.choice([0, 0, 1, 1, 2]) if i != j else random.choice([0, 1])) for j in range(n)] for i in range(n)]
        rc = L.hp_c(A, 'exact')
        rp = L.hp_py(A)
        rb = L.hp_brute(A)
        check(rc['H'] == rp['H'] == rb == H_memo(A), 'H engines trial %d' % trial)
        check(rc['HC'] == rp['HC'], 'HC engines')
        check(rc['arc'] == rp['arc'] and rc['start'] == rp['start'] and rc['end'] == rp['end'], 'arc engines')
        rq = L.hp_c(A, 'par')
        check(rq['H'] == rc['H'] % 2 and all(rq['arc'][e] == rc['arc'][e] % 2 for e in rc['arc']), 'parity engine')
        Hl, cl = lean_parity(A, list(rc['arc'].keys()))
        check(Hl == rc['H'] % 2 and all(cl[e] == rc['arc'][e] % 2 for e in rc['arc']), 'lean engine')
    print('  60 random multidigraphs (n <= 7, loops and multiplicities): all five engines agree')
    A = L.clock_matrix(7, [1, 2, 4])
    check(L.hp_c(A)['H'] == 189 and L.hp_c(A)['HC'] == 24, 'QR7')
    check(19 * rooted(L.clock_matrix(19, L.qr_set(19))) == 1172695746915, 'H(QR19) repo value')
    print('  H(QR_7) = 189, HC = 24; H(QR_19) = 19 x rooted = 1172695746915 (matches HYP-9028 record)')

    # ------------------------------------------------------------------
    section('B. L1 reflection law on all clocks of odd order')
    total = 0
    for m in [3, 5, 7, 9, 11, 13]:
        G = zm(m)
        cnt = 0
        for rr in range(0, m):
            for C in itertools.combinations(range(1, m), rr):
                A = L.clock_matrix(m, C)
                res = L.hp_c(A)
                F = half_paths(*G, list(C))
                Mg = multiplier_group(m, C)
                check(all(v % 2 == 0 for v in res['arc'].values()), 'L1(1) m=%d C=%s' % (m, C))
                check((res['H'] - m * F) % (2 * m * len(Mg)) == 0, 'L1(2) m=%d C=%s' % (m, C))
                check(F % len(Mg) == 0, 'L1 |M| divides F')
                if L.is_tournament(A):
                    check((F // len(Mg)) % 2 == 1, 'L1(3)')
                cnt += 1
        total += cnt
        print('  Z/%d: %d connection sets, L1 (1)(2)(3) hold' % (m, cnt))
    check(total == 5460, 'count 5460')
    # random multisets with loops
    for t in range(300):
        m = random.choice([5, 7, 9, 11])
        C = [random.randrange(m) for _ in range(random.randint(1, 6))]
        A = L.clock_matrix(m, C)
        res = L.hp_c(A)
        F = half_paths(*zm(m), C)
        check(all(v % 2 == 0 for v in res['arc'].values()), 'L1 multiset arcs')
        check((res['H'] - m * F) % (2 * m) == 0, 'L1 multiset H')
    print('  300 random multisets with loops: L1 holds')
    # non-cyclic group Z/3 x Z/3: all 256 subsets
    elems, add = L.zm_product(3, 3)
    neg = lambda x: ((-x[0]) % 3, (-x[1]) % 3)
    nz = [e for e in elems if e != (0, 0)]
    for rr in range(0, 9):
        for C in itertools.combinations(nz, rr):
            A = L.cayley_matrix(elems, add, C)
            res = L.hp_c(A)
            F = half_paths(elems, add, neg, (0, 0), list(C))
            check(all(v % 2 == 0 for v in res['arc'].values()), 'L1 Z3xZ3 arcs')
            check((res['H'] - 9 * F) % 18 == 0, 'L1 Z3xZ3 H')
    print('  Z/3 x Z/3: all 256 connection sets, L1 holds')
    if '--long' in args:
        F23 = half_paths(*zm(23), L.qr_set(23))
        check(F23 == 63871533 and (F23 // 11) % 2 == 1 and F23 % 11 == 0, 'F(P_23)')
        print('  [--long] Paley p=23: F = 63871533 = 11 x 5806503 (odd quotient)')
    Fexp = {3: 1, 7: 9, 11: 185, 19: 573057}
    for p, Fv in Fexp.items():
        QR = L.qr_set(p)
        F = half_paths(*zm(p), QR)
        check(F == Fv, 'F(P_%d)' % p)
        H = p * rooted(L.clock_matrix(p, QR)) if p > 13 else L.hp_c(L.clock_matrix(p, QR))['H']
        Mg = (p - 1) // 2
        check((F // Mg) % 2 == 1 and (H - p * F) % (2 * p * Mg) == 0, 'Paley L1 p=%d' % p)
        print('  Paley p=%2d: F = %d, F/|M| = %d (odd), H = %d, H - pF divisible by 2p|M|' % (p, F, F // Mg, H))

    # even order: c(u->v) = #(rho_(u+v)-symmetric HPs with middle arc u->v) mod 2
    def mid_sym(m, C, u, v):
        g = (u + v) % m
        mult = {}
        for c in C:
            mult[c % m] = mult.get(c % m, 0) + 1
        k = m // 2 - 1
        cnt = 0

        def dfs(cur, used, depth, w):
            nonlocal cnt
            if depth == k:
                cnt += w
                return
            for c, mu in mult.items():
                y = (cur - c) % m
                if y in used or (g - y) % m in used or (g - y) % m == y:
                    continue
                used.add(y)
                dfs(y, used, depth + 1, w * mu)
                used.discard(y)
        if (g - u) % m != v:
            return 0
        dfs(u, {u, v}, 0, mult.get((v - u) % m, 0))
        return cnt
    narcs = 0
    for m in [4, 6, 8, 10]:
        for rr in range(1, m):
            for C in itertools.combinations(range(1, m), rr):
                res = L.hp_c(L.clock_matrix(m, C))
                for (u, v), c in res['arc'].items():
                    if u == 0:
                        check((c - mid_sym(m, C, u, v)) % 2 == 0, 'even-order L1 m=%d C=%s' % (m, C))
                        narcs += 1
    check(narcs == 2844, 'even-order arc count')
    print('  even order (Z/4, Z/6, Z/8, Z/10, all connection sets, %d arcs out of 0): c(u->v) = number of' % narcs)
    print('    rho_(u+v)-symmetric HPs with middle arc u->v (mod 2)')

    # ------------------------------------------------------------------
    section('C. L2 lexicographic law and the mod-9 clock table')

    def walks_formula(T, Ss):
        n = len(T)
        pcs = [L.path_cover_numbers(S) for S in Ss]
        sizes = [len(S) for S in Ss]
        total = 0
        # enumerate block walks by DFS with visit counts bounded by block sizes
        def dfs(last, counts, w):
            nonlocal total
            if all(counts[v] >= 1 for v in range(n)):
                prod = w
                for v in range(n):
                    k = counts[v]
                    prod *= math.factorial(k) * pcs[v].get(k, 0)
                total += prod
            for v in range(n):
                if counts[v] < sizes[v] and (last is None or (v != last and T[last][v])):
                    counts[v] += 1
                    dfs(v, counts, w * (1 if last is None else T[last][v]))
                    counts[v] -= 1
        dfs(None, [0] * n, 1)
        return total

    for t in range(40):
        n = random.randint(2, 3)
        T = [[(random.choice([0, 1, 1, 2]) if i != j else 0) for j in range(n)] for i in range(n)]
        Ss = []
        for v in range(n):
            s = random.randint(1, 3)
            Ss.append([[(random.choice([0, 1, 1]) if i != j else 0) for j in range(s)] for i in range(s)])
        A = L.lex_product_multi(T, Ss)
        check(walks_formula(T, Ss) == L.hp_c(A)['H'], 'L2b formula trial %d' % t)
    print('  L2b composition formula = direct count on 40 random compositions (missing/doubled pairs, multiplicities)')
    # arc parity corollary on random tournaments (inner arcs; corrected cross-arc formula)
    def rt(k):
        M = [[0] * k for _ in range(k)]
        for i in range(k):
            for j in range(i + 1, k):
                if random.random() < 0.5:
                    M[i][j] = 1
                else:
                    M[j][i] = 1
        return M

    def pc2_end(S, x, start=False):
        n_ = len(S)
        seen, cnt = set(), 0
        for perm in itertools.permutations(range(n_)):
            for cut in range(1, n_):
                P1, P2 = perm[:cut], perm[cut:]
                if all(S[P1[i]][P1[i + 1]] for i in range(len(P1) - 1)) and \
                        all(S[P2[i]][P2[i + 1]] for i in range(len(P2) - 1)):
                    key = frozenset([P1, P2])
                    if key in seen:
                        continue
                    seen.add(key)
                    if (start and x in (P1[0], P2[0])) or ((not start) and x in (P1[-1], P2[-1])):
                        cnt += 1
        return cnt

    def walk_counts(T, u, v):
        n_ = len(T)
        Nw = {}

        def dfs(walk, counts):
            if all(counts[w] >= 1 for w in range(n_)) and all(counts[w] == 1 for w in range(n_) if w not in (u, v)):
                marks = sum(1 for i in range(len(walk) - 1) if walk[i] == u and walk[i + 1] == v)
                if marks:
                    key = (counts[u], counts[v])
                    Nw[key] = Nw.get(key, 0) + marks
            for w in range(n_):
                lim = 2 if w in (u, v) else 1
                if counts[w] < lim and (not walk or (w != walk[-1] and T[walk[-1]][w])):
                    counts[w] += 1
                    walk.append(w)
                    dfs(walk, counts)
                    walk.pop()
                    counts[w] -= 1
        dfs([], [0] * n_)
        return Nw
    naive_fail = 0
    ncross = 0
    for t in range(40):
        n, s = random.randint(2, 4), random.randint(2, 4)
        T, S = rt(n), rt(s)
        A = L.lex_product(T, S)
        rA, rT, rS = L.hp_c(A), L.hp_c(T), L.hp_c(S)
        check(rA['H'] % 2 == 1, 'Redei lex')
        for (x, y), c in rA['arc'].items():
            u, a = divmod(x, s)
            v, b = divmod(y, s)
            if u == v:
                check(c % 2 == rS['arc'][(a, b)] % 2, 'L2 inner arc parity')
            else:
                Nw = walk_counts(T, u, v)
                pe = {1: rS['end'][a], 2: pc2_end(S, a)}
                ps = {1: rS['start'][b], 2: pc2_end(S, b, start=True)}
                pred = sum(Nw.get((i, j), 0) * pe[i] * ps[j] for i in (1, 2) for j in (1, 2)) % 2
                check(c % 2 == pred, 'L2 cross arc parity (corrected)')
                ncross += 1
                if c % 2 != (rT['arc'][(u, v)] * rS['end'][a] * rS['start'][b]) % 2:
                    naive_fail += 1
    check(naive_fail > 0, 'naive (1,1)-only cross formula is refuted')
    print('  L2 inner-arc parity and the corrected cross-arc parity (terms k_u, k_v in {1,2}) hold on 40 random')
    print('    lexicographic products (%d cross arcs); the naive (1,1)-only formula fails on %d of them (REFUTED)'
          % (ncross, naive_fail))
    # lexicographic products are not all-odd in the tested ranges
    def parse_tourn(st, n_):
        M = [[0] * n_ for _ in range(n_)]
        k = 0
        for i in range(n_):
            for j in range(i + 1, n_):
                if st[k] == '1':
                    M[i][j] = 1
                else:
                    M[j][i] = 1
                k += 1
        return M
    TT2 = [[0, 1], [0, 0]]
    QR7v = [[1 if (x != y and ((y - x) % 7) in (1, 2, 4)) else 0 for y in range(1, 7)] for x in range(1, 7)]
    nprod = 0
    for n_ in range(2, 8):
        outs = subprocess.run(['gentourng', '-q', str(n_)], capture_output=True, text=True, check=True).stdout.split()
        for st in outs:
            T = parse_tourn(st, n_)
            prods = [L.lex_product(T, TT2), L.lex_product(TT2, T)]
            if n_ <= 4:
                prods.append(L.lex_product(T, QR7v))
            for A in prods:
                ok, _ = all_odd_full(A)
                check(not ok, 'lex product all-odd?')
                nprod += 1
    print('  no all-odd lexicographic product among T[TT_2], TT_2[T] (|T| <= 7), T[QR_7 - v] (|T| <= 4): %d products' % nprod)
    table = [([1, 4, 7], 648, 72), ([2, 5, 8], 648, 72), ([1, 3, 4, 7], 3159, 207), ([1, 3, 4, 6, 7], 14256, 1152),
             ([1, 8], 18, 2), ([1, 2, 4, 5, 7, 8], 37584, 3168), ([1, 4, 7, 8], 2268, 154), ([1, 2, 3, 4], 3267, 222)]
    for C, H, HC in table:
        r = L.hp_c(L.clock_matrix(9, C))
        check(r['H'] == H and r['HC'] == HC, 'mod-9 table %s' % C)
        print('  Clk(9,%-18s) H = %6d  HC = %5d' % (C, H, HC))
    # C3[S] closed form
    def c3_formula(S):
        pc = L.path_cover_numbers(S)
        P = lambda k: math.factorial(k) * pc.get(k, 0)
        tot = 0
        for q in range(1, len(S) + 1):
            for rho in range(3):
                tot += P(q + 1) ** rho * P(q) ** (3 - rho)
        return 3 * tot
    K3 = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
    C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    E3 = [[0] * 3 for _ in range(3)]
    check(c3_formula(E3) == 648 and c3_formula(C3) == 3159 and c3_formula(K3) == 14256, 'C3[S] formula')
    smirnov = sum(1 for w in set(itertools.permutations('aaabbbccc')) if all(w[i] != w[i + 1] for i in range(8)))
    check(smirnov == 174 and 216 * smirnov == 37584, 'Smirnov 174')
    # L2a: squares clock = C3[3K1] etc. (nauty isomorphism)
    lex_sq = L.lex_product(C3, E3)
    check(L.canon_digraph(L.clock_matrix(9, [1, 4, 7]))[0] == L.canon_digraph(lex_sq)[0], 'Sq9 = C3[3K1]')
    check(L.canon_digraph(L.clock_matrix(9, [1, 3, 4, 7]))[0] == L.canon_digraph(L.lex_product(C3, C3))[0], 'C3[C3]')
    print('  C3[S] closed form gives 648, 3159, 14256; Smirnov count 174, 216 x 174 = 37584; L2a isomorphisms (nauty)')
    # lex all-even closure (Cayley): C3[C3] arcs even
    check(all(v % 2 == 0 for v in L.hp_c(L.clock_matrix(9, [1, 3, 4, 7]))['arc'].values()), 'C3[C3] all-even')

    # ------------------------------------------------------------------
    section('D. L3 power-residue clocks mod 3^k')
    for k in range(2, 7):
        M = 3 ** k
        U = [x for x in range(M) if x % 3]
        sq = {pow(x, 2, M) for x in U}
        cu = {pow(x, 3, M) for x in U}
        check(sq == {x for x in U if x % 3 == 1}, 'squares = 1 mod 3, k=%d' % k)
        check(cu == {x for x in U if x % 9 in (1, 8)}, 'cubes = +-1 mod 9, k=%d' % k)
        prod = {(a * b) % M for a in sq for b in cu}
        check(prod == set(U), 'Sq Cu = U')
        check(len(sq & cu) == 3 ** (k - 2), 'intersection order')
        for a in range(0, 3):
            for j in range(0, k):
                d = 2 ** a * 3 ** j
                Pd = {pow(x, d, M) for x in U}
                if a >= 1:
                    pred = {x for x in U if x % (3 ** (j + 1)) == 1}
                else:
                    pred = {x for x in U if x % (3 ** (j + 1)) in (1, 3 ** (j + 1) - 1)}
                check(Pd == pred, 'd-th powers k=%d d=%d' % (k, d))
    print('  k = 2..6: squares = 1 mod 3, cubes = +-1 mod 9, all d = 2^a 3^j power sets as predicted;')
    print('            Sq Cu = U(3^k), |Sq cap Cu| = 3^(k-2) (direct product only for k = 2)')
    for k in [2, 3]:
        M = 3 ** k
        cubes = sorted({pow(x, 3, M) for x in range(M) if x % 3})
        blow = L.lex_product([[0, 1, 0, 0, 0, 0, 0, 0, 1], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0],
                              [0, 0, 1, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 1, 0, 0],
                              [0, 0, 0, 0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [1, 0, 0, 0, 0, 0, 0, 1, 0]],
                             [[0] * (3 ** (k - 2)) for _ in range(3 ** (k - 2))])
        check(L.canon_digraph(L.clock_matrix(M, cubes))[0] == L.canon_digraph(blow)[0], 'cube clock blow-up k=%d' % k)
    sqk = sorted({pow(x, 2, 27) for x in range(27) if x % 3})
    check(L.canon_digraph(L.clock_matrix(27, sqk))[0] ==
          L.canon_digraph(L.lex_product(C3, [[0] * 9 for _ in range(9)]))[0], 'square clock mod 27')
    print('  cube clock mod 9, 27 = undirected C_9 blown up; square clock mod 27 = C_3[9K_1] (nauty)')
    check(all(x == (pow(x, 3, 9) * pow(x, 4, 9)) % 9 for x in range(9) if x % 3), 'x = x^3 x^4 mod 9')
    check(sorted((x * x) % 9 for x in range(9)) == [0, 0, 0, 1, 1, 4, 4, 7, 7], 'owner squares multiset')
    check([(x * x) % 9 for x in range(10)] == [0, 1, 4, 0, 7, 7, 0, 4, 1, 0], 'owner square sequence')
    check([(x ** 3) % 9 for x in range(4)] == [0, 1, 8, 0], 'owner cube sequence')
    print('  x = x^3 x^4 (mod 9); squares mod 9 = 0,1,4,0,7,7,0,4,1,0; cubes = 0,1,8,0')

    # ------------------------------------------------------------------
    section('E. L4 multiplicative law and twisted spectral sum')
    import cmath
    def arcs(m, C):
        return sorted((x, (x + c) % m) for x in range(m) for c in C)
    sq9, cu9 = [1, 4, 7], [1, 8]
    AB = [(a * b) % 9 for a in sq9 for b in cu9]
    union = sorted(sum(([((a * x) % 9, (a * ((x + b) % 9)) % 9) for x in range(9) for b in cu9] for a in sq9), []))
    check(arcs(9, AB) == union, 'L4 union')
    union2 = sorted(sum(([((a * x) % 9, (a * ((x + b) % 9)) % 9) for x in range(9) for b in sq9] for a in cu9), []))
    check(arcs(9, AB) == union2, 'L4 union other order')
    for (m, A_, B_) in [(9, sq9, cu9), (13, [1, 3, 9], [1, 12]), (21, [1, 4, 16], [1, 2, 8, 11])]:
        prodAB = [(a * b) % m for a in A_ for b in B_]
        for j in range(m):
            lam = lambda C, jj: sum(cmath.exp(2j * cmath.pi * jj * c / m) for c in C)
            check(abs(lam(prodAB, j) - sum(lam(B_, (a * j) % m) for a in A_)) < 1e-9, 'twisted sum')
    print('  Clk(9, Sq Cu) = union of a Clk(9,Cu) over a in Sq = union of c Clk(9,Sq) over c in Cu; twisted spectral sums')

    # ------------------------------------------------------------------
    section('F. L5 complement duality H(D) = sum_k (-1)^(n-k) k! pc_k(complement)')
    for t in range(60):
        n = random.randint(2, 7)
        A = [[(1 if (i != j and random.random() < 0.5) else 0) for j in range(n)] for i in range(n)]
        Cm = L.complement(A)
        pc = L.path_cover_numbers(Cm)
        val = sum((-1) ** (n - k) * math.factorial(k) * pc.get(k, 0) for k in range(1, n + 1))
        H = L.hp_c(A)['H']
        check(val == H, 'L5 identity')
        check(H % 2 == L.hp_c(Cm)['H'] % 2, 'L5 mod 2')
    print('  60 random digraphs: identity exact; H(D) = H(complement) mod 2')

    # ------------------------------------------------------------------
    section('G. two-sheet clocks and Redei-rigid (all-odd) tournaments')
    for q in [7, 11, 19, 23]:
        SA, X = paley_sheet_data(q)
        check(L.canon_digraph(paley_minus_0(q))[0] == L.canon_digraph(T_sheet((q - 1) // 2, SA, X))[0],
              'Prop 1.1 q=%d' % q)
    print('  Prop 1.1: QR_q - 0 = T(S_A, X) in discrete-log coordinates for q = 7, 11, 19, 23 (nauty)')
    # T1 vertical arcs
    nv = 0
    for r in [3, 5, 7, 9]:
        pairs = [(d, r - d) for d in range(1, (r + 1) // 2)]
        for bitsS in itertools.product([0, 1], repeat=len(pairs)):
            S = {p[b] for p, b in zip(pairs, bitsS)}
            for z in (0, 1):
                for bits in itertools.product([0, 1], repeat=len(pairs)):
                    X = set([0] if z else [])
                    for p, b in zip(pairs, bits):
                        if b:
                            X |= {p[0], p[1]}
                    A = T_sheet(r, S, X)
                    check(L.is_tournament(A), 'two-sheet tournament')
                    res = L.hp_c(A, 'par')
                    for j in range(r):
                        e = (2 * j, 2 * j + 1) if A[2 * j][2 * j + 1] else (2 * j + 1, 2 * j)
                        check(res['arc'][e] == 1, 'T1 vertical arc odd')
                        nv += 1
    check(nv == 5688, 'T1 count')
    print('  T1: all %d vertical arcs of all symmetric reversed-sheet two-sheet clocks, r <= 9, are odd' % nv)

    # exhaustive census of general two-sheet clocks on Z/2r
    def family(N):
        evens = list(range(2, N, 2))
        pairs, seen = [], set()
        for d in evens:
            if d in seen:
                continue
            pairs.append((d, (-d) % N))
            seen |= {d, (-d) % N}
        odds = list(range(1, N, 2))
        for e0 in itertools.product([0, 1], repeat=len(pairs)):
            for e1 in itertools.product([0, 1], repeat=len(pairs)):
                for o in itertools.product([0, 1], repeat=len(odds)):
                    C0 = [p[b] for p, b in zip(pairs, e0)] + [d for d, b in zip(odds, o) if b]
                    C1 = [p[b] for p, b in zip(pairs, e1)] + [(-d) % N for d, b in zip(odds, o) if not b]
                    yield C0, C1
    expect = {6: (32, 12, 45, '3'), 10: (512, 40, 15745, '5'), 14: (8192, 84, 24540117, '7')}
    for N, (tot_e, odd_e, H_e, aut_e) in expect.items():
        tot = 0
        classes = {}
        nodd = 0
        for C0, C1 in family(N):
            A = L.twisted_clock(N, C0, C1)
            tot += 1
            ok, _ = all_odd_full(A)
            if ok:
                nodd += 1
                c, g = L.canon_digraph(A)
                if c not in classes:
                    classes[c] = (g, L.hp_c(A)['H'])
        check(tot == tot_e and nodd == odd_e and len(classes) == 1, 'census N=%d' % N)
        (g, H), = classes.values()
        check(H == H_e and g == aut_e, 'census class N=%d' % N)
        r = N // 2
        check(nodd == (r - 1) * N, 'labelled count (r-1)N')
        print('  N = %2d: %5d two-sheet clocks, %3d all-odd labelled = (r-1)N, ONE class: H = %d, |Aut| = %s'
              % (N, tot, nodd, H, g))
    # D_7: the N = 14 witness, three engines
    A = L.twisted_clock(14, [2, 4, 6, 7, 9, 11], [1, 8, 9, 10, 11, 12, 13])
    check(L.is_tournament(A), 'D7 tournament')
    rc, rp = L.hp_c(A, 'exact'), L.hp_py(A)
    check(rc['H'] == rp['H'] == 24540117 and rc['HC'] == 1001369, 'D7 H')
    check(rc['arc'] == rp['arc'] and len(rc['arc']) == 91 and all(v % 2 == 1 for v in rc['arc'].values()), 'D7 arcs')
    check(len(set(rc['arc'].values())) == 7 and min(rc['arc'].values()) == 3085307 and max(rc['arc'].values()) == 4087295,
          'D7 arc values')
    H0 = H_memo(A)
    check(H0 == 24540117, 'D7 memo H')
    for (u, v), c in rc['arc'].items():
        B = [row[:] for row in A]
        B[u][v] = 0
        check(H0 - H_memo(B) == c, 'D7 memo arc')
    check(L.canon_digraph(A)[0] == L.canon_digraph(D_r(7))[0], 'witness = D_7')
    print('  N = 14 witness C0 = {2,4,6,7,9,11}, C1 = {1,8,...,13} = D_7: H = 24540117, HC = 1001369,')
    print('    all 91 arcs odd (7 values, 3085307..4087295) by C 128-bit, Python DP and memoized H(T)-H(T-e)')
    # D_r family r = 3..11 (two engines), exact H for r <= 9
    Hexp = {3: 45, 5: 15745, 7: 24540117, 9: 116670839805}
    for r in [3, 5, 7, 9, 11]:
        A = D_r(r)
        check(L.is_tournament(A), 'D_r tournament')
        ok, res = all_odd_full(A)
        reps = arc_orbit_reps(r, A)
        Hl, cl = lean_parity(A, reps)
        check(ok and Hl == 1 and all(cl.values()), 'D_%d all-odd' % r)
        line = '  D_%-2d (N = %2d): all-odd by the full parity engine and the lean engine (%d orbit reps)' % (r, 2 * r, len(reps))
        if r in Hexp:
            re = L.hp_c(A, 'exact')
            check(re['H'] == Hexp[r] and all(v % 2 == 1 for v in re['arc'].values()), 'D_r exact')
            line += '; exact H = %d' % re['H']
        print(line)
        # interval control
        m = (r - 1) // 2
        okI, resI = all_odd_full(T_sheet(r, set(range(1, m + 1)), set(range(m))))
        check(okI == (r % 4 == 3), 'interval control r=%d' % r)
    print('  control: interval cross set {0..m-1} is all-odd exactly for r = 3 mod 4 (r = 3,5,7,9,11)')
    for r in [3, 5]:
        A = D_r(r)
        rp = L.hp_py(A)
        check(rp['H'] == Hexp[r] and all(v % 2 == 1 for v in rp['arc'].values()), 'D_r python DP')
    print('  D_3, D_5 also by the pure-Python DP')
    # arc parities are not affine in X over GF(2) (r = 3, 5)
    for r in [3, 5]:
        S = set(range(1, (r - 1) // 2 + 1))
        N = 2 * r
        prs = [(u, v) for u in range(N) for v in range(u + 1, N)]

        def pv(X):
            A = T_sheet(r, S, X)
            res = L.hp_c(A, 'par')
            return [res['arc'][(u, v)] if A[u][v] else res['arc'][(v, u)] for (u, v) in prs]
        base = pv(set())
        lin = [[a ^ b for a, b in zip(pv({i}), base)] for i in range(r)]
        nonaff = 0
        for bits in itertools.product([0, 1], repeat=r):
            pred = base[:]
            for i in range(r):
                if bits[i]:
                    pred = [a ^ b for a, b in zip(pred, lin[i])]
            if pred != pv({i for i in range(r) if bits[i]}):
                nonaff += 1
        check(nonaff > 0, 'not affine r=%d' % r)
    print('  arc parities of T(S, X) are not affine in X over GF(2) (r = 3, 5)')
    # (Q all-even, Q - v all-odd) over all 7-vertex tournaments
    outs = subprocess.run(['gentourng', '-q', '7'], capture_output=True, text=True, check=True).stdout.split()
    comb = {}
    for st in outs:
        Q = parse_tourn(st, 7)
        qe = all(v == 0 for v in L.hp_c(Q, 'par')['arc'].values())
        for v in range(7):
            keep = [i for i in range(7) if i != v]
            ok, _ = all_odd_full([[Q[i][j] for j in keep] for i in keep])
            comb[(qe, ok)] = comb.get((qe, ok), 0) + 1
    check(comb == {(True, True): 7, (True, False): 189, (False, True): 23, (False, False): 2973}, '(Q,v) census')
    print('  all 456 seven-vertex tournaments x 7 vertices: (Q all-even, Q-v all-odd) = both 7, even only 189,')
    print('    odd only 23, neither 2973')
    # quadrant form of D_r: a->a iff Im z > 0, b->b iff Im z < 0, a_j->b_k (k != j) iff Re z > 0,
    # vertical a_j -> b_j iff r = 3 mod 4, where z = exp(2 pi i (k-j)/r)
    import cmath as _cm
    for r in range(3, 40, 2):
        N = 2 * r
        Q = [[0] * N for _ in range(N)]
        for j in range(r):
            for k in range(r):
                z = _cm.exp(2j * _cm.pi * (k - j) / r)
                if k != j:
                    if z.imag > 0:
                        Q[2 * j][2 * k] = 1
                    if z.imag < 0:
                        Q[2 * j + 1][2 * k + 1] = 1
                    if z.real > 0:
                        Q[2 * j][2 * k + 1] = 1
                    else:
                        Q[2 * k + 1][2 * j] = 1
                else:
                    if r % 4 == 3:
                        Q[2 * j][2 * j + 1] = 1
                    else:
                        Q[2 * j + 1][2 * j] = 1
        check(Q == D_r(r), 'quadrant form r=%d' % r)
    print('  quadrant form of D_r (signs of Im z, -Im z, Re z; vertical arcs by r mod 4) checked for r = 3..39')
    # identities with Paley, apex completion
    check(L.canon_digraph(D_r(3))[0] == L.canon_digraph(paley_minus_0(7))[0], 'D3 = QR7-v')
    check(L.canon_digraph(D_r(5))[0] == L.canon_digraph(paley_minus_0(11))[0], 'D5 = QR11-v')
    check(L.canon_digraph(D_r(9))[0] != L.canon_digraph(paley_minus_0(19))[0], 'D9 != QR19-v')
    check(L.canon_digraph(D_r(11))[0] != L.canon_digraph(paley_minus_0(23))[0], 'D11 != QR23-v')
    for r in [3, 5, 7, 9, 11]:
        A = D_r(r)
        N = len(A)
        Q = [row[:] + [0] for row in A] + [[0] * (N + 1)]
        for v in range(N):
            if v % 2 == 0:
                Q[v][N] = 1
            else:
                Q[N][v] = 1
        check(L.is_tournament(Q) and len(set(sum(row) for row in Q)) == 1, 'Q_r regular')
        g = L.canon_digraph(Q)[1]
        gD = L.canon_digraph(A)[1]
        check(gD == str(r), '|Aut D_r| = r')
        if r == 3:
            check(L.canon_digraph(Q)[0] == L.canon_digraph(L.clock_matrix(7, [1, 2, 4]))[0], 'Q3 = QR7')
        elif r == 5:
            check(L.canon_digraph(Q)[0] == L.canon_digraph(L.clock_matrix(11, L.qr_set(11)))[0], 'Q5 = QR11')
        else:
            check(g == str(r), '|Aut Q_r| = r')
    print('  D_3 = QR_7 - v, D_5 = QR_11 - v, D_9 != QR_19 - v, D_11 != QR_23 - v; |Aut D_r| = r;')
    print('  apex completion Q_r regular, Q_3 = QR_7, Q_5 = QR_11, |Aut Q_r| = r for r = 7, 9, 11')
    # r = 9 reversed-sheet census over all sheets and all X (orbit reps)
    r = 9
    def tsets(r):
        pairs = [(d, r - d) for d in range(1, (r + 1) // 2)]
        for bits in itertools.product([0, 1], repeat=len(pairs)):
            yield frozenset(p[b] for p, b in zip(pairs, bits))
    units = [u for u in range(1, r) if math.gcd(u, r) == 1]
    seen, reps = set(), []
    for SA in tsets(r):
        for bits in itertools.product([0, 1], repeat=r):
            X = [i for i in range(r) if bits[i]]
            best = None
            for u in units:
                S2 = tuple(sorted((u * s) % r for s in SA))
                X1 = [(u * x) % r for x in X]
                for c in range(r):
                    key = (S2, tuple(sorted((x + c) % r for x in X1)))
                    if best is None or key < best:
                        best = key
            if best not in seen:
                seen.add(best)
                reps.append(best)
    check(len(reps) == 176, 'r=9 reps')
    cls = set()
    for SA, X in reps:
        A = T_sheet(r, set(SA), set(X))
        ok, _ = all_odd_full(A)
        if ok:
            cls.add(L.canon_digraph(A)[0])
    check(cls == {L.canon_digraph(D_r(9))[0], L.canon_digraph(paley_minus_0(19))[0]}, 'r=9 census classes')
    print('  r = 9: all sheets x all cross sets (176 orbit reps): exactly two all-odd classes, D_9 and QR_19 - v')
    if '--long' in args:
        r = 11
        S = set(range(1, 6))
        seenX, hits = set(), []
        for bits in itertools.product([0, 1], repeat=r):
            X = [i for i in range(r) if bits[i]]
            key = min(tuple(sorted((x + c) % r for x in X)) for c in range(r))
            if key in seenX:
                continue
            seenX.add(key)
            ok, _ = all_odd_full(T_sheet(r, S, set(key)))
            if ok:
                hits.append(key)
        check(len(seenX) == 188 and len(hits) == 6, 'r=11 census')
        isoc = {}
        for X in hits:
            isoc.setdefault(L.canon_digraph(T_sheet(r, S, set(X)))[0], []).append(X)
        check(len(isoc) == 3, 'r=11 iso classes')
        check(all(c != L.canon_digraph(paley_minus_0(23))[0] for c in isoc), 'r=11 not Paley')
        print('  [--long] r = 11 rotational sheets: 188 translation classes of X, 6 all-odd, 3 isomorphism classes:',
              list(isoc.values()))
    if '--n26' in args:
        r = 13
        A = D_r(r)
        reps = arc_orbit_reps(r, A)
        Hl, cl = lean_parity(A, reps)
        check(Hl == 1 and len(reps) == 25 and all(cl.values()), 'D_13 lean')
        AT = L.transpose(A)
        Hl2, cl2 = lean_parity(AT, [(v, u) for (u, v) in reps])
        check(Hl2 == 1 and all(cl2.values()), 'D_13 transposed')
        AI = T_sheet(r, set(range(1, 7)), set(range(6)))
        repsI = arc_orbit_reps(r, AI)
        okI = lean_parity(AI, repsI)[1]
        check(len(repsI) == 25 and sum(okI.values()) == 23, 'D_13 interval control')
        print('  [--n26] D_13 (N = 26): all 25 arc orbits odd on T and on T^op; interval control 23/25 odd')

    # one-sheet clocks: which circulant tournaments become all-odd after deleting a vertex
    def circ_scan(m):
        pairs = [(d, m - d) for d in range(1, (m + 1) // 2)]
        units = [u for u in range(1, m) if math.gcd(u, m) == 1]
        seen, reps = set(), []
        for bits in itertools.product([0, 1], repeat=len(pairs)):
            C = [p_[b] for p_, b in zip(pairs, bits)]
            key = min(tuple(sorted((u * c) % m for c in C)) for u in units)
            if key not in seen:
                seen.add(key)
                reps.append(key)
        hits = []
        for C in reps:
            A = L.clock_matrix(m, C)
            ok, _ = all_odd_full([[A[i][j] for j in range(1, m)] for i in range(1, m)])
            if ok:
                hits.append(C)
        return len(reps), hits
    circ_ms = [7, 11, 15, 19] + ([23] if '--long' in args else [])
    for m in circ_ms:
        nrep, hits = circ_scan(m)
        if L.is_prime(m):
            check(hits == [tuple(sorted(L.qr_set(m)))], 'circulant minus vertex m=%d' % m)
        else:
            check(hits == [], 'circulant minus vertex m=%d' % m)
        print('  circulants on Z/%d (%d classes up to multipliers): minus a vertex all-odd only for %s'
              % (m, nrep, 'Paley' if hits else 'none'))
    # Paley-type products chi_p(a) chi_q(b) on Z/p x Z/q (p = 3 mod 4, q = 1 mod 4): not doubly regular
    def chi(a, p_):
        a %= p_
        return 0 if a == 0 else (1 if pow(a, (p_ - 1) // 2, p_) == 1 else -1)
    for (p_, q_, Sq, H_e, odd_e) in [(3, 5, {1, 2}, 197728485, 73), (3, 5, {1, 3}, 197094945, 41)]:
        elems = [(a, b) for a in range(p_) for b in range(q_)]
        C = [(a, b) for (a, b) in elems if (a and b and chi(a, p_) * chi(b, q_) == 1) or
             (b == 0 and a and chi(a, p_) == 1) or (a == 0 and b in Sq)]
        A = L.cayley_matrix(elems, lambda x, y: ((x[0] + y[0]) % p_, (x[1] + y[1]) % q_), C)
        check(L.is_tournament(A), 'product tournament')
        n = len(A)
        common = {sum(1 for w in range(n) if A[u][w] and A[v][w]) for u in range(n) for v in range(u + 1, n)}
        check(len(common) > 1, 'not doubly regular')
        check(L.hp_c(A)['H'] == H_e, 'product H')
        B = [[A[i][j] for j in range(1, n)] for i in range(1, n)]
        check(sum(L.hp_c(B, 'par')['arc'].values()) == odd_e, 'product minus vertex')
    print('  Paley-type products on Z/3 x Z/5 (chi_3(a) chi_5(b), two choices on 0 x Z/5): not doubly regular,')
    print('    H = 197728485 / 197094945, minus a vertex 73/91 and 41/91 arcs odd')
    # apex completions: Q_r all-even (r <= 9); Q_7 - v all-odd for every v; Q_9 - a_0 not
    def Q_r(r):
        A = D_r(r)
        N = len(A)
        Q = [row[:] + [0] for row in A] + [[0] * (N + 1)]
        for v in range(N):
            if v % 2 == 0:
                Q[v][N] = 1
            else:
                Q[N][v] = 1
        return Q
    for r in [3, 5, 7, 9]:
        Q = Q_r(r)
        rq = L.hp_c(Q, 'par')
        check(all(v == 0 for v in rq['arc'].values()), 'Q_r all-even r=%d' % r)
        n = len(Q)
        for u, name in [(n - 1, 'inf'), (0, 'a0'), (1, 'b0')]:
            keep = [i for i in range(n) if i != u]
            B = [[Q[i][j] for j in keep] for i in keep]
            nodd = sum(L.hp_c(B, 'par')['arc'].values())
            if r <= 7 or name == 'inf':
                check(nodd == len(B) * (len(B) - 1) // 2, 'Q_r - v all-odd')
            else:
                check(nodd == 97, 'Q_9 - a0 not all-odd')
    print('  Q_r all-even for r = 3,5,7,9; Q_7 - v all-odd for v = inf, a_0, b_0; Q_9 - a_0 has 97/153 odd arcs')
    # sextic cyclotomic tournaments on Z/19
    p = 19
    g = L.primitive_root(p)
    cls = [sorted({pow(g, 6 * k + i, p) for k in range(3)}) for i in range(6)]
    found = {}
    for bits in itertools.product([0, 1], repeat=3):
        I = [i + 3 * b for i, b in zip(range(3), bits)]
        C = sorted(set().union(*[cls[i] for i in I]))
        A = L.clock_matrix(p, C)
        check(L.is_tournament(A), 'cyclotomic tournament')
        cert, grp = L.canon_digraph(A)
        H = p * rooted(A)
        B = [[A[i][j] for j in range(1, p)] for i in range(1, p)]
        nodd = sum(L.hp_c(B, 'par')['arc'].values())
        found.setdefault(cert, set()).add((H, grp, nodd))
    vals = sorted(list(v)[0] for v in found.values())
    check(len(found) == 2 and vals == [(1167595581285, '57', 117), (1172695746915, '171', 153)], 'sextic table')
    print('  Z/19 sextic-union tournaments: 2 classes; Paley H = 1172695746915 (|Aut| 171, minus a vertex all-odd)')
    print('    versus mixed unions H = 1167595581285 (|Aut| 57, minus a vertex 117/153 odd)')
    # mod 7 versus mod 9: the same Z/6 = squares x cubes on a field and on a ring
    U7 = list(range(1, 7))
    sq7 = {pow(x, 2, 7) for x in U7}
    cu7 = {pow(x, 3, 7) for x in U7}
    check(sq7 == {1, 2, 4} and cu7 == {1, 6} and {(a * b) % 7 for a in sq7 for b in cu7} == set(U7), 'U(7)')
    check(L.is_tournament(L.clock_matrix(7, sorted(sq7))) and not L.is_tournament(L.clock_matrix(9, [1, 4, 7])), 'field vs ring')
    check({pow(2, k, 7) for k in range(3)} == sq7 and {pow(2, k, 9) for k in range(6)} == {1, 2, 4, 5, 7, 8}, '<2>')
    print('  mod 7: squares {1,2,4} = <2> (Paley, trivial-cycle code), cubes {1,6}; mod 9: squares {1,4,7} = <4>,')
    print('    cubes {1,8}, <2> = U(9): the same Z/3 x Z/2 clock group, a tournament on the field, a blow-up on the ring')
    if '--long' in args:
        def run_batch(mats):
            inp = [str(len(mats))]
            for A in mats:
                inp.append(str(len(A)))
                inp.extend(' '.join(map(str, r)) for r in A)
            out = subprocess.run([L.HP_BIN, 'batchpar'], input='\n'.join(inp) + '\n', capture_output=True,
                                 text=True, check=True).stdout.split('\n')
            return [(int(x.split()[1]), int(x.split()[2]), int(x.split()[3])) for x in out if x.strip()]

        def tsets_(r):
            pairs = [(d, r - d) for d in range(1, (r + 1) // 2)]
            for bits in itertools.product([0, 1], repeat=len(pairs)):
                yield set(p_[b] for p_, b in zip(pairs, bits))

        def fam5():
            for S1 in tsets_(5):
                for S2 in tsets_(5):
                    for Z in range(32):
                        Zs = {d for d in range(5) if Z >> d & 1}
                        for att in range(64):
                            for ff in range(8):
                                A = [[0] * 13 for _ in range(13)]
                                for i in range(5):
                                    for j in range(5):
                                        if (j - i) % 5 in S1:
                                            A[i][j] = 1
                                        if (j - i) % 5 in S2:
                                            A[5 + i][5 + j] = 1
                                        if (j - i) % 5 in Zs:
                                            A[i][5 + j] = 1
                                        else:
                                            A[5 + j][i] = 1
                                for k in range(3):
                                    f = 10 + k
                                    for i in range(5):
                                        if att >> (2 * k) & 1:
                                            A[f][i] = 1
                                        else:
                                            A[i][f] = 1
                                        if att >> (2 * k + 1) & 1:
                                            A[f][5 + i] = 1
                                        else:
                                            A[5 + i][f] = 1
                                for bit, (u, v) in enumerate([(10, 11), (10, 12), (11, 12)]):
                                    if ff >> bit & 1:
                                        A[u][v] = 1
                                    else:
                                        A[v][u] = 1
                                yield A

        def fam11():
            for S in tsets_(11):
                for att in range(4):
                    for ff in range(2):
                        A = [[0] * 13 for _ in range(13)]
                        for i in range(11):
                            for j in range(11):
                                if (j - i) % 11 in S:
                                    A[i][j] = 1
                        for k in range(2):
                            for i in range(11):
                                if att >> k & 1:
                                    A[11 + k][i] = 1
                                else:
                                    A[i][11 + k] = 1
                        if ff:
                            A[11][12] = 1
                        else:
                            A[12][11] = 1
                        yield A

        def fam7():
            six = [parse_tourn(st, 6) for st in subprocess.run(['gentourng', '-q', '6'], capture_output=True,
                                                               text=True, check=True).stdout.split()]
            for S in tsets_(7):
                for T6 in six:
                    for att in range(64):
                        A = [[0] * 13 for _ in range(13)]
                        for i in range(7):
                            for j in range(7):
                                if (j - i) % 7 in S:
                                    A[i][j] = 1
                        for a in range(6):
                            for b in range(6):
                                A[7 + a][7 + b] = T6[a][b]
                            for i in range(7):
                                if att >> a & 1:
                                    A[7 + a][i] = 1
                                else:
                                    A[i][7 + a] = 1
                        yield A
        for name, gen, tot_e, best_e in [('Z/11 + 2', fam11(), 256, 22), ('Z/7 + 6', fam7(), 28672, 54),
                                         ('Z/5 x 2 + 3', fam5(), 262144, 72)]:
            cnt = best = found_ = 0
            batch = []
            for A in gen:
                batch.append(A)
                cnt += 1
                if len(batch) == 4096:
                    for h, odd, tot in run_batch(batch):
                        check(h == 1 and tot == 78, 'N=13 batch')
                        best = max(best, odd)
                        found_ += (odd == tot)
                    batch = []
            for h, odd, tot in run_batch(batch):
                best = max(best, odd)
                found_ += (odd == tot)
            check(cnt == tot_e and found_ == 0 and best == best_e, 'N=13 family %s' % name)
            print('  [--long] N = 13, %-12s: %6d tournaments, none all-odd, best %d/78 odd arcs' % (name, cnt, best))

    # ------------------------------------------------------------------
    section('G2. inclusion-exclusion engines (Karp/Bax over GF(2)[e], O(1) memory) and D_r beyond N = 26')
    rng2 = random.Random(31)
    for r in [3, 5, 7, 9, 11, 13]:
        m = (r - 1) // 2
        tests = [('D_r', D_r(r)), ('interval', T_sheet(r, set(range(1, m + 1)), set(range(m))))]
        for k in range(3):
            pairs = [(d, r - d) for d in range(1, (r + 1) // 2)]
            S = {p_[rng2.randint(0, 1)] for p_ in pairs}
            X = {x for x in range(r) if rng2.random() < 0.5}
            tests.append(('random%d' % k, T_sheet(r, S, X)))
        line = []
        for name, A in tests:
            reps = arc_orbit_reps(r, A)
            H3, C3, R3 = karp_run(KARP3_BIN, A, reps)
            H4, C4, R4 = karp_run(KARP4_BIN, A, reps)
            check((H3, C3, R3) == (H4, C4, R4) and H3 == 1, 'karp3 = karp4 r=%d %s' % (r, name))
            if r <= 11:
                full = L.hp_c(A, 'par')
                check(all(full['arc'][e] == C3[e] for e in reps), 'karp = subset DP r=%d %s' % (r, name))
            else:
                Hl, cl = lean_parity(A, reps) if '--n26' in args else (None, None)
                if cl is not None:
                    check(all(cl[e] == C3[e] for e in reps), 'karp = lean r=13 %s' % name)
            if name == 'D_r':
                check(all(C3.values()), 'D_%d all-odd (karp)' % r)
            line.append('%s %d/%d' % (name, sum(C3.values()), len(reps)))
        print('  r = %2d: karp3 = karp4%s; odd orbits: %s' % (r, ' = subset DP' if r <= 11 else (' = lean' if '--n26' in args else ''), ', '.join(line)))
    if '--long' in args:
        A = D_r(15)
        reps = arc_orbit_reps(15, A)
        H3, C3, R3 = karp_run(KARP3_BIN, A, reps)
        H4, C4, R4 = karp_run(KARP4_BIN, A, reps)
        check(H3 == H4 == 1 and C3 == C4 and all(C3.values()) and len(reps) == 29, 'D_15')
        print('  [--long] D_15 (N = 30): all 29 arc orbits odd by karp3 and karp4 (%d subset reps)' % R3)
    if '--n34' in args:
        A = D_r(17)
        reps = arc_orbit_reps(17, A)
        H3, C3, R3 = karp_run(KARP3_BIN, A, reps)
        check(H3 == 1 and len(reps) == 33 and all(C3.values()), 'D_17 karp3')
        AT = L.transpose(A)
        H4, C4, R4 = karp_run(KARP4_BIN, AT, [(v, u) for (u, v) in reps])
        check(H4 == 1 and all(C4.values()), 'D_17^op karp4')
        print('  [--n34] D_17 (N = 34): all 33 arc orbits odd by karp3 on T and karp4 on T^op (%d subset reps)' % R3)
    if '--n38' in args:
        A = D_r(19)
        reps = arc_orbit_reps(19, A)
        H4, C4, R4 = karp_run(KARP4_BIN, A, reps)
        check(H4 == 1 and len(reps) == 37 and all(C4.values()), 'D_19 karp4')
        print('  [--n38] D_19 (N = 38): all 37 arc orbits odd by karp4 (%d subset reps)' % R4)

    # ------------------------------------------------------------------
    section('H. the Syracuse clock')

    def syrac_numerators(K):
        laws = [({0: 1}, 1)]
        num, den = {0: 1}, 1
        for n in range(K):
            M, Mn, Lp = 3 ** (n + 1), 3 ** n, 2 * 3 ** n
            new = {}
            for x in range(M):
                if x % 3 == 0:
                    continue
                s = 0
                for a in range(1, Lp + 1):
                    t = (pow(2, a, M) * x) % M
                    if t % 3 != 1:
                        continue
                    c = num.get(((t - 1) // 3) % Mn, 0)
                    if c:
                        s += c << (Lp - a)
                if s:
                    new[x] = s
            num, den = new, den * ((1 << Lp) - 1)
            laws.append((num, den))
        return laws
    K = 8
    laws = syrac_numerators(K)
    num, den = laws[2]
    tao = [Fr(num.get(x, 0), den) for x in range(9)]
    check(tao == [Fr(0), Fr(8, 63), Fr(16, 63), Fr(0), Fr(11, 63), Fr(4, 63), Fr(0), Fr(2, 63), Fr(22, 63)], 'Tao mod 9')
    print('  Syrac(Z/9) = (0, 8, 16, 0, 11, 4, 0, 2, 22)/63  (Tao 2019, after Lemma 1.12) reproduced')
    for k in range(1, K + 1):
        num, den = laws[k]
        M = 3 ** k
        check(sum(num.values()) == den and all(x % 3 for x in num), 'mass / class 0 null')
        check(Fr(sum(c for x, c in num.items() if x % 3 == 1), den) == Fr(1, 3), 'class-1 mass')
        check(all(num.get((2 * x) % M, 0) == 2 * num.get(x, 0) for x in range(1, M, 3)), 'G2 doubling k=%d' % k)
        if k >= 2:
            pnum, pden = laws[k - 1]
            Mp = 3 ** (k - 1)
            proj = {}
            for x, c in num.items():
                proj[x % Mp] = proj.get(x % Mp, 0) + c
            check(all(Fr(proj.get(y, 0), den) == Fr(pnum.get(y, 0), pden) for y in range(Mp)), 'consistency (1.23)')
            inv4 = pow(4, -1, M)
            for x in range(1, M, 3):
                lhs = 4 * Fr(num.get((x * inv4) % M, 0), den)
                y3 = ((x - 1) // 3) % Mp
                y6 = (((x - 1) // 3) * pow(2, -1, Mp)) % Mp
                rhs = Fr(num.get(x, 0), den)
                rhs += Fr(pnum.get(y3, 0), pden) if y3 % 3 == 1 else 0
                rhs += 2 * Fr(pnum.get(y6, 0), pden) if y6 % 3 == 1 else 0
                check(lhs == rhs, 'G3 one-tick k=%d' % k)
    print('  k = 1..%d: class 0 null; class 1 has mass 1/3; G2 P(2x) = 2P(x) for x = 1 mod 3; G3 one-tick IFS identity;' % K)
    print('            projection consistency (Tao (1.23))')
    # 3-state clock at level 9
    a = Fr(8, 21); b = Fr(11, 21); c = Fr(2, 21)
    check(a == b / 4 + Fr(1, 4) and b == c / 4 + Fr(1, 2) and c == a / 4, '3-state square clock')
    print('  level 9: mu1(1) = mu1(4)/4 + 1/4, mu1(4) = mu1(7)/4 + 1/2, mu1(7) = mu1(1)/4  ->  (8, 11, 2)/21')
    # coupling one tick
    law = {x: tao[x] for x in range(9)}
    c0 = {s: law[pow(4, s, 9)] * 3 for s in range(3)}
    c1 = {s: law[(-pow(4, s, 9)) % 9] * Fr(3, 2) for s in range(3)}
    check(all(c1[s] == c0[(s + 1) % 3] for s in range(3)), 'one-tick coupling')
    check(law[1] != (law[1] + law[4] + law[7]) * (law[1] + law[8]), 'not a product law')
    print('  coupling: P(s | c=1) = P(s+1 | c=0) = (11, 2, 8)/21; the law is not a product of its marginals')
    # P^2 = 1 pi
    units = [1, 2, 4, 5, 7, 8]
    P = {x: {y: Fr(0) for y in units} for x in units}
    for x in units:
        for aa in range(1, 7):
            y = ((3 * x + 1) * pow(pow(2, aa, 9), -1, 9)) % 9
            P[x][y] += Fr(1, 2 ** aa) / (1 - Fr(1, 64))
    check(len(set(tuple(P[x][y] for y in units) for x in units)) == 2, 'rank 2 rows')
    P2 = {x: tuple(sum(P[x][z] * P[z][y] for z in units) for y in units) for x in units}
    check(all(P2[x] == tuple(law[y] for y in units) for x in units), 'P^2 = 1 pi')
    print('  mod-9 residue chain: two distinct rows (out-twins over the sign bit), P^2 = 1 pi exactly')
    # step law on integers
    coord = {}
    for cc in (0, 1):
        for s in range(3):
            coord[((-1) ** cc * pow(4, s, 9)) % 9] = (cc, s)
    bad = 0
    for A_ in range(1, 2 * 10 ** 6, 2):
        t = 3 * A_ + 1
        v = (t & -t).bit_length() - 1
        B_ = t >> v
        exp = (v % 2, v % 3) if A_ % 3 == 0 else (v % 2, (1 + coord[A_ % 9][0] + v) % 3)
        if coord[B_ % 9] != exp:
            bad += 1
    check(bad == 0, 'step law integers')
    print('  step law (c,s) -> (v mod 2, 1+c+v mod 3) holds for all 10^6 odd A < 2*10^6')
    # general level: x = (-1)^c (-2)^l mod 81, l in Z/27; S(A): c' = v mod 2, l' = lambda(A) - v, lambda = -A mod 3
    M81 = 81
    logm2 = {}
    xx = 1
    for l in range(27):
        logm2[xx] = l
        xx = (xx * (-2)) % M81
    check(len(logm2) == 27, '-2 generates 1-units mod 81')
    def sl(x):
        c = 0 if x % 3 == 1 else 1
        return c, logm2[(x * (-1) ** c) % M81]
    for A_ in range(1, 200001, 2):
        t = 3 * A_ + 1
        v = (t & -t).bit_length() - 1
        B_ = t >> v
        lam = logm2[t % M81]
        check(lam % 3 == (-A_) % 3, 'lambda = -A mod 3')
        check(sl(B_) == (v % 2, (lam - v) % 27), 'general step law mod 81')
    print('  mod 81: S(A) = (-1)^v (-2)^(lambda(A) - v), lambda(A) = log_(-2)(3A+1) = -A (mod 3): all odd A < 2*10^5')

    bad = cnt = 0
    for A_ in range(1, 4 * 10 ** 5, 2):
        a_ = A_
        vs, xs = [], [a_]
        for _ in range(12):
            t = 3 * a_ + 1
            v = (t & -t).bit_length() - 1
            a_ = t >> v
            vs.append(v)
            xs.append(a_)
        for i in range(1, 12):
            pred = ((-1) ** (vs[i] % 2) * pow(4, (1 + vs[i - 1] % 2 + vs[i]) % 3, 9)) % 9
            cnt += 1
            if xs[i + 1] % 9 != pred:
                bad += 1
    check(bad == 0 and cnt == 2200000, 'two-hand formula')
    print('  A_(n+1) mod 9 = (-1)^(v_n) 4^(1 + v_(n-1) + v_n): all 2,200,000 orbit transitions')
    rng = random.Random(7)
    Ns = 200000
    cntl = {}
    for _ in range(Ns):
        a_ = rng.getrandbits(80) | 1
        for _ in range(6):
            t = 3 * a_ + 1
            a_ = t >> ((t & -t).bit_length() - 1)
        cntl[a_ % 9] = cntl.get(a_ % 9, 0) + 1
    tv = sum(abs(cntl.get(x, 0) / Ns - float(tao[x])) for x in range(9))
    check(tv < 0.01, 'empirical law')
    print('  empirical law of S^6(A) mod 9, 200000 random 80-bit A: TV distance %.4f to Tao' % tv)

    # ------------------------------------------------------------------
    section('I. triplet instances')
    # Pythagorean
    npyth = 0
    for m_ in range(2, 101):
        for n_ in range(1, m_):
            if (m_ - n_) % 2 == 1 and math.gcd(m_, n_) == 1:
                a_, b_, c_ = m_ * m_ - n_ * n_, 2 * m_ * n_, m_ * m_ + n_ * n_
                issq = lambda z: z >= 0 and math.isqrt(z) ** 2 == z
                check(c_ % 2 == 1 and a_ % 2 == 1 and b_ % 4 == 0, 'parities')
                check(issq(c_ + b_) and issq(c_ - b_), 'c +- b squares')
                check((c_ + a_) % 2 == 0 and issq((c_ + a_) // 2) and issq((c_ - a_) // 2), 'c +- a twice squares')
                check(not issq(c_ + a_) and not issq(c_ - a_) or a_ == 0, 'c +- a not squares')
                npyth += 1
    print('  %d primitive triples (m < 101): c +- (even leg) squares, c +- (odd leg) twice squares, never squares' % npyth)
    # squares mod 9 under squaring
    sqmap = {x: (x * x) % 9 for x in [1, 4, 7]}
    check(sqmap == {1: 1, 4: 7, 7: 4}, 'squaring involution')
    check(all((2 * x) % 3 == (-x) % 3 for x in range(3)), '2 = -1 mod 3')
    check([p for p in range(3, 200) if L.is_prime(p) and (2 + 1) % p == 0] == [3], 'only p=3 has 2=-1')
    print('  squaring on {1,4,7}: 1 fixed, 4 <-> 7; doubling on Z/3 = negation; p = 3 is the only prime with 2 = -1')
    # F = S + U palindrome, S = 2^U - 1, solutions p, p^3, p^2qr
    LIM = 100000
    spf = list(range(LIM + 1))
    for i in range(2, int(LIM ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, LIM + 1, i):
                if spf[j] == j:
                    spf[j] = i
    shapes = {}
    for N in range(2, LIM + 1):
        x, ex = N, []
        while x > 1:
            p_ = spf[x]
            e = 0
            while x % p_ == 0:
                x //= p_
                e += 1
            ex.append(e)
        r_ = len(ex)
        tau = 1
        for e in ex:
            tau *= e + 1
        F_ = tau - 2
        sqf = all(e == 1 for e in ex)
        S_ = 2 ** r_ - 1 - (1 if sqf else 0)
        U_ = r_ - (1 if (r_ == 1 and ex[0] == 1) else 0)
        sq_cont = F_ - S_
        pal = (U_ == sq_cont)
        check((F_ == S_ + U_) == pal, 'palindrome')
        if not sqf and not (r_ == 1 and ex[0] == 1):
            check(S_ == 2 ** U_ - 1, 'S = 2^U - 1')
        if F_ == S_ + U_:
            shapes[tuple(sorted(ex, reverse=True))] = shapes.get(tuple(sorted(ex, reverse=True)), 0) + 1
    check(set(shapes) == {(1,), (3,), (2, 1, 1)}, 'F=S+U shapes')
    N = 2 * 2 * 3 * 5
    divs = [d for d in range(2, N) if N % d == 0]
    check(len(divs) == 10, '10 divisors')
    primes_ = [d for d in divs if L.is_prime(d)]
    sqfc = [d for d in divs if not L.is_prime(d) and all(d % (q * q) for q in (2, 3, 5))]
    sqc = [d for d in divs if any(d % (q * q) == 0 for q in (2, 3, 5))]
    check((len(primes_), len(sqfc), len(sqc)) == (3, 4, 3), '3-4-3')
    Om = lambda d: sum(1 for q in (2, 3, 5) for k in range(1, 3) if d % (q ** k) == 0)
    ranks = [sum(1 for d in divs if Om(d) == k) for k in (1, 2, 3)]
    check(ranks == [3, 4, 3], 'rank 3-4-3')
    check(sorted(sqc) != sorted(d for d in divs if Om(d) == 3), 'different partitions')
    check(sorted(d for d in divs if Om(d) == 2) == [4, 6, 10, 15] and sorted(sqfc) == [6, 10, 15, 30], 'p^2 <-> pqr')
    print('  N <= 10^5: F = S + U  <=>  (U, S-U, F-S) palindromic; S = 2^U - 1 off squarefree;')
    print('             solutions exactly p, p^3, p^2qr; 60 = 2^2*3*5: split 3-4-3 and rank split 3-4-3 differ by p^2 <-> pqr')
    for p_ in [q for q in range(5, 200) if L.is_prime(q)]:
        check(pow(p_, 3, 9) == (1 if p_ % 3 == 1 else 8) and pow(p_, 2, 9) in (1, 4, 7), 'p^3 = (p/3) mod 9')
    print('  primes 5 <= p < 200: p^3 = (p/3) (mod 9), p^2 in {1,4,7}')

    print('\n%d checks, %.1f s' % (NCHECK, time.time() - T0))
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
