"""procgen_smallgraph_20260926_t3.py -- T3: small graphs that encode arithmetic.

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).  Called by the runner.
New instance: the Collatz-alphabet sum graph C_n (x ~ y iff x + y is a power of 2 or a power
of 3), whose Hamiltonicity windows are cut out by powers of 3 and by the gaps |2^p - 3^a|.
"""
import bisect
import random
from math import isqrt, log

from procgen_smallgraph_20260926_lib import (
    check, family, sum_graph, leaves, components, verify_path, forcing_ham, cpsat_ham,
    ham_exist, ham_count, ham_heur, squares, powers, pow23, perfect_powers)
from procgen_smallgraph_20260926_t1 import compress, decide_path, decide_cycle, local_obstruction


# ---------------------------------------------------------------------------------------
# T3.A  the Collatz parity graph: degree carries no threshold information
# ---------------------------------------------------------------------------------------
def parity_graph_indegrees(q, k, sigma):
    M = 1 << k
    half = M >> 1
    indeg = [0] * M
    for s in range(M):
        if s % 2 == 0:
            t = (s // 2) % half
        else:
            t = ((q * s + sigma[s]) // 2) % half
        indeg[t] += 1
        indeg[t + half] += 1
    return indeg


def t3a():
    print('== T3.A parity graph of q n +- 1 sign strategies: degrees ==')
    ok = True
    for q in (3, 5, 7):
        for k in range(2, 13):
            for sg in (1, -1):
                sigma = {s: sg for s in range(1, 1 << k, 2)}
                if set(parity_graph_indegrees(q, k, sigma)) != {2}:
                    ok = False
    check('T3.A1 constant signs (Collatz q=3, and 3n-1, 5n+-1, 7n+-1): every node has in-degree 2 and out-degree 2',
          ok, 'levels k=2..12: the parity graph is 2-regular (de Bruijn type); degree cannot carry a threshold')
    rng = random.Random(20260926)
    rng_ok = True
    seen = set()
    for trial in range(300):
        q = rng.choice((3, 5, 7, 9))
        k = rng.randint(3, 10)
        sigma = {s: rng.choice((1, -1)) for s in range(1, 1 << k, 2)}
        ind = parity_graph_indegrees(q, k, sigma)
        seen |= set(ind)
        if not (set(ind) <= {1, 2, 3} and sum(ind) == 2 << k and ind.count(3) == ind.count(1)):
            rng_ok = False
    check('T3.A2 general sign strategies: in-degrees in {1,2,3}, #3 = #1 (merge pairs, THM-4481)', rng_ok and seen == {1, 2, 3},
          '300 random strategies; flips merge target pairs, which is the only degree information there is')


# ---------------------------------------------------------------------------------------
# T3.B  lacunary targets give forests;  T3.C  powers of p split by valuation
# ---------------------------------------------------------------------------------------
def vp(x, p):
    c = 0
    while x % p == 0:
        x //= p
        c += 1
    return c


def t3bc():
    print('== T3.B lacunary target sets and T3.C valuation splitting ==')
    rng = random.Random(7)
    sets = {'pow2': powers(2, 4000), 'pow3': powers(3, 4000), 'pow5': powers(5, 4000)}
    for trial in range(5):
        v, S = rng.randint(3, 9), []
        while v < 4000:
            S.append(v)
            v = int(v * rng.uniform(2.0, 3.5)) + 1
        sets[f'random_lacunary_{trial}'] = S
    ok = True
    for name, S in sets.items():
        Ss = sorted(S)
        lac = all(Ss[i + 1] >= 2 * Ss[i] for i in range(len(Ss) - 1))
        for n in range(2, 1001):
            E = sum_graph(n, set(S))
            if len(E) != n - len(components(n, E)):
                ok = False
        ok = ok and lac
    check('T3.B1 |S cap (M,2M)| <= 1 for all M  =>  G_S(n) is a forest', ok,
          f'{len(sets)} target sets (powers of 2,3,5 and 5 random lacunary sets), all n <= 1000: |E| = n - #components')
    sq = set(squares(400))
    first_cycle = None
    for n in range(2, 60):
        E = sum_graph(n, sq)
        if len(E) > n - len(components(n, E)):
            first_cycle = n
            break
    E15 = sum_graph(15, sq)
    cyc = [1, 3, 6, 10, 15]
    unicyclic = (len(E15) == 15 and len(components(15, E15)) == 1 and
                 all((min(a, b), max(a, b)) in set(E15) for a, b in zip(cyc, cyc[1:] + cyc[:1])))
    check('T3.B2 squares are not lacunary: Q_n is a forest for n <= 14, and Q_15 is unicyclic', first_cycle == 15 and unicyclic,
          'the unique cycle of Q_15 is (1,3,6,10,15) with sums 4,9,16,25,16: first cycle and first chain appear at the same n')
    ok = True
    for p in (2, 3, 5):
        S = set(powers(p, 10 ** 6))
        for n in range(2, 600):
            E = sum_graph(n, S)
            if any(vp(x, p) != vp(y, p) for x, y in E):
                ok = False
    check('T3.C1 S = powers of p: every edge joins x, y with v_p(x) = v_p(y)', ok, 'p = 2, 3, 5, n < 600 (x+y = p^j, x,y < p^j)')
    ok = True
    S = set(powers(2, 10 ** 6))
    for n in range(2, 600):
        comps = components(n, sum_graph(n, S))
        classes = {}
        for x in range(1, n + 1):
            classes.setdefault(vp(x, 2), []).append(x)
        if sorted(comps) != sorted(classes.values()):
            ok = False
    check('T3.C2 powers of 2: the components of G_S(n) are exactly the 2-adic valuation classes (each a tree)', ok,
          'n < 600; hence no Hamiltonian path for any n >= 2 (1 and 2 are separated)')
    S3 = set(powers(3, 10 ** 6))
    paths3 = [n for n in range(1, 200) if decide_path(n, sum_graph(n, S3), S3, want_methods=False)[0] == 'PATH']
    check('T3.C3 powers of 3: Hamiltonian path only for n in {1, 2}', paths3 == [1, 2],
          'for n >= 3 the multiples of 3 form a closed class (valuation)')


# ---------------------------------------------------------------------------------------
# T3.D  consecutive targets: the Catalan / Gersonides zigzag
# ---------------------------------------------------------------------------------------
def consecutive_zigzag(m):
    """path on [1..m] using sums m and m+1: m, 1, m-1, 2, m-2, ..."""
    seq, lo, hi = [], 1, m
    turn = 1
    seq.append(m)
    hi = m - 1
    while lo <= hi:
        if turn:
            seq.append(lo)
            lo += 1
        else:
            seq.append(hi)
            hi -= 1
        turn ^= 1
    return seq


def t3d():
    print('== T3.D consecutive targets: Catalan and Gersonides zigzags ==')
    ok = True
    for m in range(2, 601):
        z = consecutive_zigzag(m)
        T = {m, m + 1}
        if not (verify_path(m, z, T) and verify_path(m - 1, z[1:], T)):
            ok = False
    check('T3.D1 m, m+1 in S  =>  explicit zigzag Hamiltonian paths of G_S(m) and G_S(m-1)', ok,
          'm = 2..600; e.g. m=8: ' + str(consecutive_zigzag(8)))
    pp = perfect_powers(10 ** 7)
    cons = [a for a, b in zip(pp, pp[1:]) if b == a + 1 and a >= 2]
    check('T3.D2 consecutive perfect powers up to 1e7', cons == [8], '(8,9) only; for all sizes this is Mihailescu (CITED)')
    PP = set(perfect_powers(10 ** 4))
    rows = []
    good = True
    for n in (7, 8):
        E = sum_graph(n, PP)
        cnt, comp, eps = ham_count(n, E)
        rows.append(f'n={n}: {cnt} chains, ends {sorted(eps)}')
        good = good and comp
    E8 = sum_graph(8, PP)
    cnt8, _, _ = ham_count(8, E8)
    firstpp = min(n for n in range(3, 40) if decide_path(n, sum_graph(n, PP), PP, want_methods=False)[0] == 'PATH')
    check('T3.D3 perfect-power sums: the first chains (n=7,8) come from the Catalan pair 8,9',
          good and firstpp == 7 and cnt8 == 1 and verify_path(8, consecutive_zigzag(8)[::-1], PP),
          '; '.join(rows) + '; the n=8 chain is the (8,9) zigzag [4,5,3,6,2,7,1,8] (sums 9,8,9,...)')
    C = set(pow23(10 ** 7))
    cons23 = [(a, b) for a, b in zip(sorted(C), sorted(C)[1:]) if b == a + 1]
    check('T3.D4 consecutive targets among powers of 2 and 3', cons23 == [(1, 2), (2, 3), (3, 4), (8, 9)],
          'Gersonides (THM-4484): these give zigzag chains of C_n at n = 1, 2, 3, 7, 8')


# ---------------------------------------------------------------------------------------
# T3.E  the Collatz-alphabet sum graph C_n
# ---------------------------------------------------------------------------------------
def leaf_profile(tlist, N):
    """L[n] = #vertices of degree <= 1 of G_S(n), I[n] = #isolated, via sigma1/sigma2."""
    diffL = [0] * (N + 3)
    diffI = [0] * (N + 3)
    for x in range(1, N + 1):
        i = bisect.bisect_right(tlist, x)
        c = []
        while len(c) < 2 and i < len(tlist):
            if tlist[i] != 2 * x:
                c.append(tlist[i])
            i += 1
        g2 = (c[1] - x) if len(c) > 1 else 10 ** 18
        g1 = (c[0] - x) if c else 10 ** 18
        hi = min(N, g2 - 1)
        if hi >= x:
            diffL[x] += 1
            diffL[hi + 1] -= 1
        hi = min(N, g1 - 1)
        if hi >= x:
            diffI[x] += 1
            diffI[hi + 1] -= 1
    L, I, a, b = [0] * (N + 1), [0] * (N + 1), 0, 0
    for n in range(1, N + 1):
        a += diffL[n]
        b += diffI[n]
        L[n], I[n] = a, b
    return L, I


def predicted_windows(amax):
    wins = [(2, 3)]
    for a in range(2, amax + 1):
        B = 3 ** (a - 1)
        P = [2 ** p for p in range(0, 4 * a + 4) if B < 2 ** p < 3 * B]
        if len(P) == 2:
            lo = max(3 * B - P[0], P[1] - B)
            hi = 3 * B - 1
        else:
            lo = 2 * B
            hi = 3 * B
        wins.append((lo, hi))
    return wins


def windows_of(flags, N):
    out, n = [], 2
    while n <= N:
        if flags[n]:
            m = n
            while m + 1 <= N and flags[m + 1]:
                m += 1
            out.append((n, m))
            n = m + 1
        else:
            n += 1
    return out


def t3e(NHAM=2200, AMAX=13):
    print('== T3.E the Collatz-alphabet sum graph C_n (sums = powers of 2 or 3) ==')
    N = 3 ** AMAX
    tl = sorted(set(pow23(4 * N + 10)))
    ts = set(tl)
    # E1: the largest power of 2 <= n has degree <= 1
    ok = True
    for n in range(2, 20001):
        x = 1 << (n.bit_length() - 1)
        d = sum(1 for t in tl[bisect.bisect_right(tl, x):bisect.bisect_right(tl, x + n)] if t != 2 * x)
        if d > 1:
            ok = False
    check('T3.E1 the largest power of two 2^p <= n has degree <= 1 in C_n', ok,
          'hand: targets in (2^p, 2^p+n] lie in (2^p, 3*2^p): 2^(p+1)=2x is excluded, 2^(p+2) is too big, at most one power of 3; '
          'checked n <= 20000.  Hence C_n has NO Hamiltonian cycle (n >= 3) and 2^p ends every Hamiltonian path')
    L, I = leaf_profile(tl, N)
    ok2 = True
    for a in range(2, AMAX + 1):
        B = 3 ** (a - 1)
        for n in range(2 * B, 3 * B):
            if not (L[n] <= 2 and I[n] == 0):
                ok2 = False
        onepow = len([p for p in range(0, 4 * a + 4) if B < 2 ** p < 3 * B]) == 1
        if (L[3 * B] <= 2 and I[3 * B] == 0) != onepow:
            ok2 = False
        if 3 * B + 1 <= N and L[3 * B + 1] < 3:
            ok2 = False
    check('T3.E2 for 2*3^(a-1) <= n < 3^a: at most two leaves (powers of 2 in (n/3,n]) and no isolated vertex; '
          'n = 3^a qualifies iff (3^(a-1),3^a) holds one power of 2; n = 3^a + 1 has >= 3 leaves',
          ok2, f'hand proof in the note; verified for a = 2..{AMAX} (n <= {N})')
    flags = [False] * (N + 1)
    for n in range(2, N + 1):
        flags[n] = (L[n] <= 2 and I[n] == 0)
    win = windows_of(flags, N)
    pred = predicted_windows(AMAX)
    check('T3.E3 the candidate windows (<= 2 leaves, no isolated vertex) of C_n', win == pred,
          f'n <= 3^{AMAX}: ' + ' '.join(f'[{a},{b}]' for a, b in win) +
          '; left ends = max(3^a - 2^p, 2^(p+1) - 3^(a-1)) when (3^(a-1),3^a) holds 2^p<2^(p+1), else 2*3^(a-1)')
    gaps = []
    for (lo, hi), a in zip(pred[1:], range(2, AMAX + 1)):
        B = 3 ** (a - 1)
        reps = []
        if lo == 2 * B:
            reps.append(f'2*3^{a - 1}')
        for p in range(0, 4 * a + 4):
            if lo == 3 ** a - 2 ** p:
                reps.append(f'3^{a}-2^{p}')
            if lo == 2 ** p - 3 ** (a - 1):
                reps.append(f'2^{p}-3^{a - 1}')
        if reps:
            gaps.append(f'{lo}=' + '='.join(reps))
    check('T3.E3b the left ends as Collatz gaps', len(gaps) == AMAX - 1, ' '.join(gaps))
    # E4: exact Hamiltonian-path set for n <= NHAM.  Every non-local "no" is decided by the exact forcing
    # solver; CP-SAT confirms all of them for n <= NCP and a deterministic sample (every 8th) above.
    NCP = 1000
    P, NONE_local, NONE_exact = [], 0, []
    cnt_cp = 0
    for n in range(1, NHAM + 1):
        T = set(t for t in tl if t <= 2 * n)
        E = sum_graph(n, T)
        if n == 1:
            P.append(1)
            continue
        ob = local_obstruction(n, E)
        if ob:
            NONE_local += 1
            continue
        st, seq, nodes = forcing_ham(n, E)
        if st == 'PATH':
            assert verify_path(n, seq, T)
            x = 1 << (n.bit_length() - 1)
            assert x in (seq[0], seq[-1])
            P.append(n)
            continue
        assert st == 'NONE'
        NONE_exact.append(n)
        if n <= NCP or len(NONE_exact) % 8 == 0:
            st3, _ = cpsat_ham(n, E)
            assert st3 == 'NONE', (n, st3)
            cnt_cp += 1
    small = sum(1 for n in NONE_exact if n <= NCP)
    check(f'T3.E4 Hamiltonian paths of C_n, n <= {NHAM} (FINITE-EXACT)',
          compress(P) == '1-3,5-8,18-26,49-63,65-66,68-80,179-194,224-243,473-575,665-728,1319-1418,1536-1616,1620-1663,1703-1713,1792-1802,1920-2114,2120-2186'
          if NHAM == 2200 else True,
          f'{compress(P)}; {NONE_local} n excluded by a local obstruction, {len(NONE_exact)} by the exact forcing solver; '
          f'CP-SAT agrees on all {small} of them with n <= {NCP} and on a sample of {cnt_cp - small} above; '
          f'the largest power of 2 <= n is an end of every witness')
    inside = all(any(lo <= n <= hi for lo, hi in pred) for n in P if n >= 2)
    check('T3.E5 every Hamiltonian n lies in a candidate window; windows are only partly Hamiltonian', inside,
          'e.g. [162,243] -> 179-194, 224-243; [473,728] -> 473-575, 665-728')
    small = [n for n in P if n <= 8]
    check('T3.E6 the smallest Hamiltonian n are the Gersonides zigzags', small == [1, 2, 3, 5, 6, 7, 8],
          'n = 1,2 (targets 2,3), 2,3 (3,4), 7,8 (8,9) by T3.D1; n = 5, 6 are the only others below 9')


# ---------------------------------------------------------------------------------------
# T3.F  hitting times: degree thresholds versus Hamiltonicity thresholds
# ---------------------------------------------------------------------------------------
def t3f():
    print('== T3.F hitting times (min-degree vs Hamiltonicity) ==')
    sq = squares(10 ** 6)
    # delta(Q_n) >= k hitting times via the bound delta > (sqrt2-1)sqrt(n) - 2
    hits = {}
    for k in range(1, 9):
        nb = int(((k + 1) / (2 ** 0.5 - 1)) ** 2) + 2   # beyond nb the bound gives delta >= k
        last_bad = 0
        for n in range(1, nb + 1):
            dmin = min(isqrt(x + n) - isqrt(x) - (1 if isqrt(2 * x) ** 2 == 2 * x else 0) for x in range(1, n + 1))
            if dmin < k:
                last_bad = n
        hits[k] = last_bad + 1
    check('T3.F1 square-sum graph: delta(Q_n) >= k for all n >= h_k', hits[2] == 31,
          ' '.join(f'k={k}:h={h}' for k, h in hits.items()) + ' (hand: delta(Q_n) > (sqrt2-1)sqrt(n) - 2, rest finite)')
    rows = []
    got = {}
    fams = [('squares', 120), ('triangular', 120), ('pentagonal', 160), ('perfect_powers', 120),
            ('sq_plus_1', 120), ('sq_minus_1', 120), ('sq_c-5', 120)]
    for fam, N in fams:
        f = family(fam)
        T = set(f(8 * N))
        last_leafy = last_nocyc = last_local = last_nopath = 0
        for n in range(3, N + 1):
            E = sum_graph(n, T)
            Lv, deg = leaves(n, E)
            if Lv:
                last_leafy = n
            if local_obstruction(n, E):
                last_local = n
            if decide_path(n, E, T, want_methods=False)[0] != 'PATH':
                last_nopath = n
            if decide_cycle(n, E, T)[0] != 'CYCLE':
                last_nocyc = n
        got[fam] = (last_leafy, last_nocyc, last_local, last_nopath)
        rows.append(f'{fam}(n<={N}): deg<=1 until {last_leafy}, no cycle until {last_nocyc}; '
                    f'local path obstruction until {last_local}, no path until {last_nopath}')
    expected = {'squares': (30, 31, 18, 24), 'triangular': (11, 14, 8, 8), 'pentagonal': (56, 56, 27, 44),
                'perfect_powers': (16, 16, 14, 14), 'sq_plus_1': (23, 31, 19, 19), 'sq_minus_1': (38, 38, 19, 30),
                'sq_c-5': (36, 36, 22, 30)}
    check('T3.F2 hitting-time table: cycles appear at or soon after min degree 2; paths lag the leaf count more',
          got == expected, ' | '.join(rows))


# ---------------------------------------------------------------------------------------
# T3.G  hybrids requested in the brief
# ---------------------------------------------------------------------------------------
def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def t3g():
    print('== T3.G hybrids ==')
    sq = set(squares(10 ** 5))
    P = []
    for n in range(1, 151):
        E = set(sum_graph(n, sq))
        for x in range(1, n + 1):
            y = T(x)
            if y <= n and y != x:
                E.add((min(x, y), max(x, y)))
        E = sorted(E)
        if n == 1:
            P.append(n)
            continue
        if local_obstruction(n, E):
            continue
        st, seq, _ = forcing_ham(n, E)
        if st == 'PATH':
            Es = set(E)
            assert all((min(a, b), max(a, b)) in Es for a, b in zip(seq, seq[1:]))
            P.append(n)
    check('T3.G1 square-sum + Collatz edges {x, T(x)}: Hamiltonian path for n <= 150 at', compress(P) == '1-6,8-11,13-150',
          f'{compress(P)} (square-sum alone: 1,15-17,23,25-...); the Collatz edges fill the low-degree holes')
    ok = []
    for k in range(2, 10):
        M = 1 << k
        sqres = {(x * x) % M for x in range(M)}
        E = set()
        for x in range(M):
            for y in range(x + 1, M):
                if (x + y) % M in sqres or (3 * x + 1) % M == y or (3 * y + 1) % M == x:
                    E.add((x + 1, y + 1))
        st, seq, _ = forcing_ham(M, sorted(E), cycle=True)
        ok.append((k, st))
    check('T3.G2 Z/2^k with edges x->3x+1 and x+y = square residue: Hamiltonian cycle',
          ok == [(2, 'NONE')] + [(k, 'CYCLE') for k in range(3, 10)],
          ' '.join(f'k={k}:{st}' for k, st in ok) + ' (dense: about 2^k/6 square residues per vertex)')
