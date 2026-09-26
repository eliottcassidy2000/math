"""procgen_smallgraph_20260926_t1.py -- T1: the square-sum problem at 15, forced endpoints,
the Leaf Lemma, the Zigzag Threshold Theorem, and target-set variants.

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).  Called by the runner
procgen_smallgraph_20260926_run.py; every printed claim is a check(...).
"""
from math import isqrt

from procgen_smallgraph_20260926_lib import (
    check, family, squares, sum_graph, leaves, components, verify_path, degree_formula,
    ham_count, ham_count_capped, ham_heur, forcing_ham, cpsat_ham, ham_exist, is_square, path_uniqueness)


def compress(lst):
    out, i = [], 0
    while i < len(lst):
        j = i
        while j + 1 < len(lst) and lst[j + 1] == lst[j] + 1:
            j += 1
        out.append(f'{lst[i]}' if i == j else f'{lst[i]}-{lst[j]}')
        i = j + 1
    return ','.join(out) if out else '(none)'


def second_target(x, tlist, exclude_double=True):
    """(sigma1, sigma2): the two smallest targets > x, excluding 2x."""
    out = []
    for t in tlist:
        if t > x and not (exclude_double and t == 2 * x):
            out.append(t)
            if len(out) == 2:
                break
    while len(out) < 2:
        out.append(10 ** 18)
    return out[0], out[1]


# ---------------------------------------------------------------------------------------
# decision helpers: every NONE beyond a hand-checkable local obstruction is confirmed by
# two independent exact methods (C exhaustive DFS and Python forcing), plus CP-SAT.
# ---------------------------------------------------------------------------------------
def local_obstruction(n, E):
    L, deg = leaves(n, E)
    iso = [v for v in range(1, n + 1) if deg[v] == 0]
    if n >= 2 and iso:
        return f'isolated vertex {iso[0]}'
    if len(L) >= 3:
        return f'three vertices of degree 1: {L[:3]}'
    if len(components(n, E)) > 1:
        return 'disconnected'
    return None


def decide_path(n, E, T, dfs_limit=20_000_000, want_methods=True):
    """Returns (status, seq, how).  status PATH or NONE; raises if undecided."""
    if n == 1:
        return 'PATH', [1], 'trivial'
    ob = local_obstruction(n, E)
    if ob:
        return 'NONE', None, ob
    st, seq, _ = ham_heur(n, E, seed=n, restarts=20, budget=20000)
    if st == 'PATH' and verify_path(n, seq, T):
        return 'PATH', seq, 'heuristic witness'
    st2, seq2, nodes2 = forcing_ham(n, E)
    if st2 == 'PATH':
        assert verify_path(n, seq2, T)
        return 'PATH', seq2, 'forcing witness'
    if st2 != 'NONE':
        raise RuntimeError(f'undecided n={n}')
    methods = ['forcing']
    st3, _ = cpsat_ham(n, E)
    if st3 != 'NONE':
        raise RuntimeError(f'CP-SAT disagrees at n={n}: {st3}')
    methods.append('cpsat')
    if want_methods:
        st1, _, info1 = ham_exist(n, E, limit=dfs_limit)
        if st1 == 'NONE':
            methods.append('dfs')
        elif st1 == 'PATH':
            raise RuntimeError(f'DFS disagrees at n={n}')
    return 'NONE', None, 'exact:' + '+'.join(methods)


def decide_cycle(n, E, T, dfs_limit=20_000_000):
    if n < 3:
        return 'NONE', None, 'n<3'
    L, deg = leaves(n, E)
    if L:
        return 'NONE', None, f'vertex {L[0]} of degree <= 1'
    st, seq, _ = ham_heur(n, E, cycle=True, seed=n, restarts=20, budget=20000)
    if st == 'CYCLE' and verify_path(n, seq, T, cycle=True):
        return 'CYCLE', seq, 'heuristic witness'
    st2, seq2, _ = forcing_ham(n, E, cycle=True)
    if st2 == 'CYCLE':
        assert verify_path(n, seq2, T, cycle=True)
        return 'CYCLE', seq2, 'forcing witness'
    if st2 != 'NONE':
        raise RuntimeError(f'undecided cycle n={n}')
    st3, _ = cpsat_ham(n, E, cycle=True)
    if st3 != 'NONE':
        raise RuntimeError(f'CP-SAT disagrees (cycle) n={n}')
    methods = ['forcing', 'cpsat']
    st1, _, _ = ham_exist(n, E, cycle=True, limit=dfs_limit)
    if st1 == 'NONE':
        methods.append('dfs')
    elif st1 == 'CYCLE':
        raise RuntimeError('DFS disagrees')
    return 'NONE', None, 'exact:' + '+'.join(methods)


# ---------------------------------------------------------------------------------------
# T0  cross-validation of the three exact methods
# ---------------------------------------------------------------------------------------
def t0():
    print('== T0 solver cross-validation ==')
    fams = ['squares', 'triangular', 'perfect_powers', 'sq_c2', 'pow2or3', 'sq_c-5']
    agree, tot = 0, 0
    for fam in fams:
        f = family(fam)
        for n in range(2, 46):
            T = set(f(2 * n))
            E = sum_graph(n, T)
            st1, _, _ = ham_exist(n, E)
            st2, seq2, _ = forcing_ham(n, E)
            c1, _, _ = ham_exist(n, E, cycle=True)
            c2, cs2, _ = forcing_ham(n, E, cycle=True)
            ok = ((st1 == 'PATH') == (st2 == 'PATH')) and ((c1 == 'CYCLE') == (c2 == 'CYCLE'))
            ok = ok and (seq2 is None or verify_path(n, seq2, T)) and (cs2 is None or verify_path(n, cs2, T, cycle=True))
            agree += ok
            tot += 1
    check('T0.1 exhaustive C DFS and forcing solver agree (paths and cycles)', agree == tot,
          f'{tot} graphs: 6 families, 2 <= n <= 45')
    bad, tot = 0, 0
    for fam in ['squares', 'sq_c2', 'sq_c-5', 'triangular', 'perfect_powers', 'sq_c7', 'sq_c-3', 'pow2or3', 'sq_plus_1']:
        f = family(fam)
        for n in range(3, 27):
            T = set(f(2 * n))
            E = sum_graph(n, T)
            cnt, comp, _ = ham_count(n, E, limit=50_000_000)
            if not comp:
                continue
            st, _, _ = path_uniqueness(n, E)
            tot += 1
            bad += st != {0: 'NONE', 1: 'UNIQUE'}.get(cnt, 'MULTIPLE')
    check('T0.2 uniqueness test (#paths >= 2 iff G-e is traceable for an edge e of one path) = exhaustive counts',
          bad == 0 and tot > 200, f'{tot} graphs')


# ---------------------------------------------------------------------------------------
# T1.A  degree formula and the complete leaf table of Q_n
# ---------------------------------------------------------------------------------------
LEAF_TABLE = {1: 7, 2: 13, 3: 5, 4: 11, 5: 10, 6: 9, 7: 8, 8: 16, 9: 15, 10: 14,
              11: 13, 12: 12, 16: 19, 17: 18, 18: 30}
ISOLATED_TABLE = {1: 2, 2: 6, 4: 4}


def t1a():
    print('== T1.A degree formula and leaves of Q_n ==')
    NMAX = 300
    sq = squares(4 * NMAX)
    sqs = set(sq)
    bad = 0
    for n in range(1, NMAX + 1):
        E = sum_graph(n, sqs)
        L, deg = leaves(n, E)
        for x in range(1, n + 1):
            f = isqrt(x + n) - isqrt(x) - (1 if is_square(2 * x) else 0)
            if f != deg[x]:
                bad += 1
    check('T1.A1 Lemma 1 (degree formula) for squares',
          bad == 0, f'deg_n(x) = floor(sqrt(x+n)) - floor(sqrt x) - [2x square] for all 1<=x<=n<={NMAX}')
    # hand inequalities (exact integer forms)
    check('T1.A2 exact forms of the counting bounds', 68 ** 2 > 2 * 48 ** 2 and 150 ** 2 > 2 * 106 ** 2,
          '68^2=4624>4608=2*48^2 gives (sqrt2-1)sqrt(x)>=2 for x>=24, so >=2 squares in (x,2x]; '
          '150^2=22500>22472=2*106^2 gives >=3 squares for x>=53')
    cnt2 = [x for x in range(24, 200001) if isqrt(2 * x) - isqrt(x) < 2]
    cnt3 = [x for x in range(53, 200001) if isqrt(2 * x) - isqrt(x) < 3]
    check('T1.A2b numeric sanity of the two counting bounds', not cnt2 and not cnt3,
          'no x in [24,2e5] with <2 squares in (x,2x], none in [53,2e5] with <3')
    # leaf table from sigma2, for x <= 60 (x >= 24 impossible by A2; 32, 50 checked explicitly)
    table = {}
    iso = {}
    for x in range(1, 61):
        s1, s2 = second_target(x, sq)
        top = s2 - x - 1
        if top >= x:
            table[x] = top
        topi = s1 - x - 1
        if topi >= x:
            iso[x] = topi
    check('T1.A3 complete leaf table of Q_n: x is a leaf exactly for x <= n <= N(x)', table == LEAF_TABLE,
          ' '.join(f'{x}:[{x},{N}]' for x, N in sorted(table.items())))
    check('T1.A4 isolated vertices of Q_n', iso == ISOLATED_TABLE,
          'x=1 for n<=2, x=2 for n<=6, x=4 at n=4 only')
    ok = True
    for n in range(1, NMAX + 1):
        E = sum_graph(n, sqs)
        L, deg = leaves(n, E)
        pred = sorted(x for x, N in LEAF_TABLE.items() if x <= n <= N)
        if L != pred:
            ok = False
    check('T1.A5 predicted leaf sets = computed leaf sets', ok,
          f'all n <= {NMAX}; in particular Q_n has minimum degree >= 2 for every n >= 31 (hand: A2 + table)')


# ---------------------------------------------------------------------------------------
# T1.B  why 15, and the anatomy of the unique chain
# ---------------------------------------------------------------------------------------
CHAIN15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]


def t1b():
    print('== T1.B why n = 15: obstruction below, anatomy at 15, 16, 17 ==')
    sqs = set(squares(200))
    reasons = {}
    for n in range(2, 15):
        E = sum_graph(n, sqs)
        L, deg = leaves(n, E)
        if n <= 6:
            ok = deg[2] == 0
            reasons[n] = 'vertex 2 isolated'
        else:
            pred = sorted(x for x, N in LEAF_TABLE.items() if x <= n <= N)
            ok = len(pred) >= 3 and L == pred
            reasons[n] = f'leaves {pred}'
        check(f'T1.B1 n={n}: no Hamiltonian path (local)', ok, reasons[n])
    n = 15
    E = sum_graph(n, sqs)
    L, deg = leaves(n, E)
    degs = {x: deg[x] for x in range(1, 16)}
    check('T1.B2 Q_15 degrees', L == [8, 9] and degs[1] == 3 and degs[3] == 3 and
          all(degs[x] == 2 for x in range(1, 16) if x not in (1, 3, 8, 9)),
          'deg 1: {8,9}; deg 3: {1,3}; all other 11 vertices deg 2')
    check('T1.B3 Q_15 has exactly 15 edges', len(E) == 15, 'sum of degrees 30')
    pedges = sorted(tuple(sorted(p)) for p in zip(CHAIN15, CHAIN15[1:]))
    rest = sorted(set(E) - {(1, 3)})
    check('T1.B4 uniqueness: Q_15 minus edge {1,3} IS the chain', verify_path(15, CHAIN15, sqs) and rest == pedges,
          'a Hamiltonian path omits exactly one of 15 edges; it must lower deg(1)=deg(3)=3, so it is {1,3}')
    cnt, comp, eps = ham_count(15, E)
    check('T1.B5 exhaustive count at n=15', cnt == 1 and comp and eps == {(8, 9): 1}, 'one chain up to reversal, ends 8 and 9')
    sums = [a + b for a, b in zip(CHAIN15, CHAIN15[1:])]
    check('T1.B6 the chain is the 3-square zigzag', sums == [9, 16, 25, 16] * 3 + [9, 16],
          f'sums {sums}')
    by = {}
    for (x, y) in E:
        by.setdefault(x + y, []).append((x, y))
    low, high = set(range(1, 9)), set(range(10, 16))
    ok = (all(x in low and y in low for x, y in by[9]) and all(x in high and y in high for x, y in by[25])
          and all(x <= 7 and y >= 9 for x, y in by[16]) and by[4] == [(1, 3)] and sorted(by) == [4, 9, 16, 25])
    check('T1.B7 edge classes of Q_15', ok,
          '9-edges inside L={1..8}; 25-edges inside H={10..15}; 16-edges join {1..7} to {9..15}; 4-edge {1,3} only')
    check('T1.B8 why the endpoints are 8 and 9', [p for p in E if 8 in p] == [(1, 8)] and [p for p in E if 9 in p] == [(7, 9)],
          '8: 16-partner is itself (2*8=4^2), 25-partner 17>15, only 9-partner 1; 9: 9-partner 0, 25-partner 16>15, only 16-partner 7')
    for n, ends, path in [(16, (8, 16), CHAIN15 + [16]), (17, (16, 17), [17] + CHAIN15 + [16])]:
        E = sum_graph(n, sqs)
        cnt, comp, eps = ham_count(n, E)
        check(f'T1.B9 n={n}: unique chain = zigzag extended', cnt == 1 and comp and list(eps) == [ends]
              and verify_path(n, path, sqs), f'{path}; forced ends {ends}')
    E = sum_graph(18, sqs)
    L, _ = leaves(18, E)
    check('T1.B10 n=18: three leaves', L == [16, 17, 18], 'window {15,16,17} closes')


# ---------------------------------------------------------------------------------------
# T1.C  Leaf Lemma for arbitrary target sets
# ---------------------------------------------------------------------------------------
VARIANTS = ['squares', 'cubes', 'triangular', 'pentagonal', 'pow2', 'pow3', 'pow2or3',
            'perfect_powers', 'sq_plus_1', 'sq_minus_1', 'sq_c2', 'sq_c-5', 'sq_c7']


def t1c():
    print('== T1.C Leaf Lemma (all target sets) ==')
    viol = 0
    viol2 = 0
    tested = 0
    for fam in VARIANTS:
        f = family(fam)
        for n in range(2, 161):
            T = sorted(set(f(4 * n + 4)))
            Ts = set(T)
            E = sum_graph(n, Ts)
            L, deg = leaves(n, E)
            Lset = set(L)
            for x in range(1, n + 1):
                if deg[x] != degree_formula(n, x, Ts):
                    viol += 1
            for x in L:
                if 2 * x in Ts:
                    continue
                below = [t for t in T if t <= x]
                si = below[-1] if below else 1
                for y in range(max(1, si), x + 1):
                    if y not in Lset:
                        viol += 1
            if len(L) <= 2:
                for x in L:
                    ok = (x in Ts) or (2 * x in Ts) or (x - 1 in Ts and x - 1 in Lset) or x <= 2
                    if not ok:
                        viol2 += 1
            tested += 1
    check('T1.C1 Lemma 1 + Leaf Lemma on 13 target families, 2<=n<=160', viol == 0,
          f'{tested} graphs: degree formula exact; every non-half-target leaf x has all of [s_i, x] as leaves')
    check('T1.C2 Corollary: with <= 2 leaves, each leaf is a target, a half-target, a target+1 leaf, or <= 2',
          viol2 == 0, 'no exception')


# ---------------------------------------------------------------------------------------
# T1.D  Zigzag Lemma and Zigzag Threshold Theorem
# ---------------------------------------------------------------------------------------
def zigzag(k):
    """Z_k for h = 2k: positions 4j+1: h-2j, 4j+2: 2j+1, 4j+3: 2h-1-2j, 4j+4: h+2+2j."""
    h = 2 * k
    n = 2 * h - 1
    seq = []
    j = 0
    while len(seq) < n:
        for v in (h - 2 * j, 2 * j + 1, 2 * h - 1 - 2 * j, h + 2 + 2 * j):
            if len(seq) < n:
                seq.append(v)
        j += 1
    return seq


def t1d():
    print('== T1.D Zigzag Lemma and Zigzag Threshold Theorem ==')
    ok = True
    for k in range(1, 401):
        z = zigzag(k)
        T = {2 * k + 1, 4 * k, 6 * k + 1}
        if not (verify_path(4 * k - 1, z, T) and z[0] == 2 * k and z[-1] == 2 * k + 1):
            ok = False
    check('T1.D1 Zigzag Lemma: {2k+1,4k,6k+1} in S => Z_k is a Hamiltonian path of G_S(4k-1) from 2k to 2k+1',
          ok, 'explicit construction verified for k = 1..400 with S = exactly these three targets')
    check('T1.D1b Z_4 is the square chain at 15', zigzag(4) == CHAIN15, f'Z_4 = {zigzag(4)}')
    # targets of S_c for c = k(4-k): 2ki + (i-2)^2
    ok = True
    for k in range(2, 200):
        c = k * (4 - k)
        S = set(t for t in (j * j + c for j in range(1, 3 * k + 10)) if t >= 3)
        pred = set(t for t in (2 * k * i + (i - 2) ** 2 for i in range(3 - k, 3 * k)) if t >= 3)
        pred = set(t for t in pred if t <= 8 * k + 4)
        if set(t for t in S if t <= 8 * k + 4) != pred or not {2 * k + 1, 4 * k, 6 * k + 1} <= S:
            ok = False
        if 3 in S:
            ok = False
    check('T1.D2 S_c, c=k(4-k): relevant targets <= 8k+4 are {4 (k>=3), 2k+1, 4k, 6k+1, 8k+4}; 3 not in S_c',
          ok, 'k = 2..199')
    # theorem, verified exhaustively for k = 2..40 (hand proof covers all k)
    rows = []
    allok = True
    for k in range(2, 41):
        c = k * (4 - k)
        f = family(f'sq_c{c}')
        T = set(f(20 * k))
        first = None
        for n in range(2, 4 * k):
            E = sum_graph(n, T)
            st, seq, how = decide_path(n, E, T, want_methods=(n > 4 * k - 8))
            if st == 'PATH':
                first = n
                break
        n = 4 * k - 1
        E = sum_graph(n, T)
        cnt, comp, eps = ham_count(n, E)
        good = first == 4 * k - 1 and cnt == 1 and comp and list(eps) == [(2 * k, 2 * k + 1)]
        allok = allok and good
        if k <= 9 or k in (20, 40):
            rows.append(f'k={k} c={c} first n={first} paths={cnt} ends={list(eps)}')
    check('T1.D3 Zigzag Threshold Theorem: for c=k(4-k) the first n>=2 with a path is 4k-1, the path is unique (=Z_k)',
          allok, 'k=2..40 exhaustive; ' + '; '.join(rows))
    # hand-proof ingredients at general k: degree structure at n = 4k-1 and the leaves below
    ok = True
    for k in range(3, 120):
        c = k * (4 - k)
        T = set(family(f'sq_c{c}')(20 * k))
        n = 4 * k - 1
        E = sum_graph(n, T)
        L, deg = leaves(n, E)
        if L != [2 * k, 2 * k + 1] or deg[1] != 3 or deg[3] != 3 or len(E) != n:
            ok = False
        if any(deg[x] != 2 for x in range(1, n + 1) if x not in (1, 3, 2 * k, 2 * k + 1)):
            ok = False
        for m in range(2 * k + 2, 4 * k - 1):
            L2, _ = leaves(m, sum_graph(m, T))
            if not {2 * k, 2 * k + 1, 2 * k + 2} <= set(L2):
                ok = False
        L3, _ = leaves(2 * k + 1, sum_graph(2 * k + 1, T))
        if not {2, 4, 2 * k} <= set(L3):
            ok = False
        for m in (2 * k - 1, 2 * k):
            L4, _ = leaves(m, sum_graph(m, T))
            if not {2, 4, 5} <= set(L4):
                ok = False
        for m in range(2, 2 * k - 1):
            _, d5 = leaves(m, sum_graph(m, T))
            if d5[2] != 0:
                ok = False
    check('T1.D4 hand-proof ingredients hold for k = 3..119', ok,
          'n<=2k-2: vertex 2 isolated; n=2k-1,2k: leaves 2,4,5; n=2k+1: leaves 2,4,2k; 2k+2<=n<=4k-2: leaves 2k,2k+1,2k+2; '
          'n=4k-1: leaves {2k,2k+1}, deg(1)=deg(3)=3, others 2, exactly n edges')
    # converse over a range of c
    sig = []
    nonsig = 0
    for c in range(-300, 101):
        f = family(f'sq_c{c}')
        T = set(f(4000))
        first = None
        for n in range(3, 400):
            E = sum_graph(n, T)
            st, seq, how = decide_path(n, E, T, want_methods=False)
            if st == 'PATH':
                first = n
                break
        if first is None:
            continue
        E = sum_graph(first, T)
        status, Pu, Qu = path_uniqueness(first, E)
        zig = status == 'UNIQUE'
        if zig:
            u, v = sorted((Pu[0], Pu[-1]))
            zig = (v == u + 1 and 2 * u in T and v in T)
        if zig:
            sig.append(c)
        else:
            nonsig += 1
    pred = [k * (4 - k) for k in range(2, 40) if -300 <= k * (4 - k) <= 100]
    check('T1.D5 converse (-300<=c<=100): the zigzag signature at the first n>=3 occurs iff c = k(4-k)',
          sorted(sig) == sorted(set(pred)), f'signature set {sorted(sig)}; {nonsig} other c without it')


# ---------------------------------------------------------------------------------------
# T1.E  Pell / Catalan
# ---------------------------------------------------------------------------------------
def t1e():
    print('== T1.E Pell and Catalan readings of (8, 9) ==')
    # zigzag triples inside the squares: h = 2s^2, h+1 = t^2, 3h+1 = u^2  (simultaneous Pell)
    t, s = 3, 2
    hits = []
    for idx in range(1, 401):
        if is_square(6 * s * s + 1):
            hits.append((s, t))
        t, s = 3 * t + 4 * s, 2 * t + 3 * s
    check('T1.E1 square zigzag triples {h+1,2h,3h+1} (h=2s^2 Pell) for the first 400 Pell solutions', hits == [(2, 3)],
          'only s=2 (h=8: 9,16,25); the system t^2-2s^2=1, u^2-6s^2=1 checked to s ~ 10^306')
    sols = [(tt, m) for m in range(1, 300) for tt in [isqrt(2 ** m + 1)] if tt * tt == 2 ** m + 1]
    check('T1.E2 t^2 - 2^m = 1 only for (3,3)', sols == [(3, 3)],
          '(t-1)(t+1)=2^m forces t-1=2, t+1=4; checked m<300')
    gers = [(p, a) for p in range(1, 200) for a in range(1, 130) if abs(2 ** p - 3 ** a) == 1]
    check('T1.E3 Gersonides |2^p-3^a|=1 (THM-4484)', gers == [(1, 1), (2, 1), (3, 2)], f'{gers}')
    # the leaf 8 comes from 2*8 = 16 = 4^2, not from 8 = 2^3: in Q_n replace the cube structure
    check('T1.E4 the only doubled square that is a power of two and one below a square is 8',
          [x for x in range(1, 10 ** 6) if (x & (x - 1)) == 0 and is_square(2 * x) and is_square(x + 1)] == [8],
          'x=2s^2=2^m with x+1 square -> x=8 (x < 10^6); by T1.E2 for all x')


# ---------------------------------------------------------------------------------------
# T1.F  square-sum thresholds and forced endpoints up to NMAX
# ---------------------------------------------------------------------------------------
A071983 = {15: 1, 16: 1, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 3, 24: 0, 25: 10, 26: 12, 27: 35,
           28: 52, 29: 19, 30: 20, 31: 349, 32: 392, 33: 669, 34: 4041, 35: 17175, 36: 12960}
A071984 = {32: 1, 33: 1, 34: 11, 35: 57, 36: 31, 37: 20, 38: 25, 39: 50, 40: 64}


def t1f(NMAX=300):
    print(f'== T1.F square-sum thresholds, counts and forced endpoints (n <= {NMAX}) ==')
    sqs = set(squares(4 * NMAX))
    P, C, NP, NC = [], [], [], []
    hows = {}
    for n in range(1, NMAX + 1):
        E = sum_graph(n, sqs)
        st, seq, how = decide_path(n, E, sqs)
        (P if st == 'PATH' else NP).append(n)
        if st == 'NONE':
            hows[n] = how
        stc, seqc, howc = decide_cycle(n, E, sqs)
        (C if stc == 'CYCLE' else NC).append(n)
        if stc == 'NONE' and n >= 3 and 'degree' not in howc:
            hows[-n] = howc
    check('T1.F1 Hamiltonian path in Q_n', P == [1, 15, 16, 17, 23] + list(range(25, NMAX + 1)),
          f'exactly n in {compress(P)}; nonlocal NONE certificates: ' +
          ', '.join(f'{n}:{hows[n]}' for n in sorted(k for k in hows if k > 0) if 'exact' in hows[n]))
    check('T1.F2 Hamiltonian cycle in Q_n', C == list(range(32, NMAX + 1)),
          f'exactly n in {compress(C)}; n<=30 has a leaf; n=31: {hows.get(-31)}')
    ok = True
    for n, v in A071983.items():
        E = sum_graph(n, sqs)
        cnt, comp, eps = ham_count(n, E)
        ok = ok and comp and cnt == v
    check('T1.F3 path counts n=15..36 equal A071983 (CITED)', ok, 'all 22 values')
    ok = True
    for n, v in A071984.items():
        E = sum_graph(n, sqs)
        cnt, comp, _ = ham_count(n, E, cycle=True)
        ok = ok and comp and cnt == v
    check('T1.F4 cycle counts n=32..40 equal A071984 (CITED)', ok, 'all 9 values')
    forced = {}
    endsets = {}
    for n in [15, 16, 17, 23] + list(range(25, 37)):
        E = sum_graph(n, sqs)
        cnt, comp, eps = ham_count(n, E)
        inter = None
        ends = set()
        for (u, v) in eps:
            inter = {u, v} if inter is None else inter & {u, v}
            ends |= {u, v}
        forced[n] = sorted(inter)
        endsets[n] = len(ends)
    check('T1.F5 forced endpoints (in every chain), exhaustive for n<=36',
          forced == {15: [8, 9], 16: [8, 16], 17: [16, 17], 23: [18], 25: [18], 26: [18], 27: [18], 28: [18],
                     29: [18], 30: [18], 31: [], 32: [], 33: [], 34: [], 35: [], 36: []},
          '15:{8,9} 16:{8,16} 17:{16,17} 23,25..30:{18} 31..36: none; #possible endpoints: ' +
          ' '.join(f'{n}:{endsets[n]}' for n in sorted(endsets)))
    ok = True
    for n in range(37, NMAX + 1):
        E = sum_graph(n, sqs)
        st, seq, _ = ham_heur(n, E, seed=7 * n, restarts=50, budget=50000)
        assert st == 'PATH' and verify_path(n, seq, sqs)
        a, b = seq[0], seq[-1]
        st2, seq2, _ = ham_heur(n, E, seed=11 * n, restarts=200, budget=50000, interior=(a, b))
        if not (st2 == 'PATH' and verify_path(n, seq2, sqs) and {seq2[0], seq2[-1]}.isdisjoint({a, b})):
            ok = False
    check('T1.F6 no forced endpoint for 37 <= n <= NMAX', ok,
          f'two chains with disjoint endpoint pairs found and verified for every n in [37,{NMAX}]')


# ---------------------------------------------------------------------------------------
# T1.G  variants: first possible n, endpoints, Catalan-type signature
# ---------------------------------------------------------------------------------------
def t1g():
    print('== T1.G variants: first possible n, endpoints, signature ==')
    specs = [('squares', 60), ('triangular', 60), ('pentagonal', 80), ('perfect_powers', 60),
             ('sq_plus_1', 80), ('sq_minus_1', 80), ('pow2', 200), ('pow3', 200), ('pow2or3', 100)]
    summary = {}
    for fam, N in specs:
        f = family(fam)
        T = set(f(8 * N))
        P, C = [], []
        for n in range(1, N + 1):
            E = sum_graph(n, T)
            st, seq, how = decide_path(n, E, T, want_methods=(n <= 40))
            if st == 'PATH':
                P.append(n)
            stc, _, _ = decide_cycle(n, E, T)
            if stc == 'CYCLE':
                C.append(n)
        nontriv = [n for n in P if n >= 3]
        first = nontriv[0] if nontriv else None
        info = ''
        sigok = False
        if first:
            E = sum_graph(first, T)
            cnt, comp, eps = ham_count(first, E)
            L, _ = leaves(first, E)
            info = f'first n>=3: {first}, chains={cnt}, end pairs={sorted(eps)}, leaves={L}'
            if cnt == 1 and len(eps) == 1:
                (u, v), = eps
                sigok = (v == u + 1 and 2 * u in T and v in T)
        summary[fam] = (compress(P), compress(C), first, sigok)
        check(f'T1.G {fam} (n <= {N})', first is None or comp,
              f'paths: {compress(P)}; cycles: {compress(C)}; {info}; zigzag signature: {sigok}')
    check('T1.G* zigzag signature at the first n>=3 among the variants',
          [f for f in summary if summary[f][3]] == ['squares', 'pow2or3'],
          'squares (n=15, Z_4 on 9,16,25) and powers of 2 or 3 (n=3, Z_1 on the Gersonides pair 3,4); no other variant')
    # cubes: threshold 305 (A304120) and cycles from 473
    f = family('cubes')
    T = set(f(4000))
    P, NP_nonlocal = [], []
    for n in range(1, 400):
        E = sum_graph(n, T)
        st, seq, how = decide_path(n, E, T, want_methods=False)
        if st == 'PATH':
            P.append(n)
        elif 'exact' in how:
            NP_nonlocal.append(n)
    check('T1.G cubes: Hamiltonian path for n < 400 exactly at', P == [1, 305, 333] + list(range(385, 400)),
          f'{compress(P)}; local obstruction for all other n<295; exact search (forcing+CP-SAT) NONE at {compress(NP_nonlocal)}')
    E = sum_graph(305, T)
    cnt, comp, eps = ham_count(305, E)
    L, _ = leaves(305, E)
    check('T1.G cubes at 305 (A304120(3)=305, CITED)', comp and cnt == 3 and L == [256] and all(256 in p for p in eps),
          f'3 chains, leaf 256 = 512/2 (half-target of 8^3) is a forced end; end pairs {sorted(eps)}')
    C = []
    for n in range(3, 601):
        E = sum_graph(n, T)
        st, seq, how = decide_cycle(n, E, T)
        if st == 'CYCLE':
            C.append(n)
    check('T1.G cubes: Hamiltonian cycle for n <= 600 exactly at', C == list(range(473, 601)),
          f'{compress(C)} (n <= 472: a vertex of degree <= 1); the cycle threshold equals the min-degree-2 threshold')
