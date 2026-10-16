"""Referee check (C): vertices of out-degree 1 (HYP-9280 'Post-audit continuation'; results file section 11.3).

w has the unique out-neighbour o, U = V - {w, o} (all in-neighbours of w), U+ = {u in U : o -> u}, U- = {u in U : u -> o}.
Methods: explicit DFS over Hamiltonian paths and brute force over all n! permutations (pa_core), never the session's
subset DP.  Modes (each a short run; the .out file is the concatenation of the runs, in the order of __main__ usage):
  identity   the PROVED identity (U- empty): C(u,w) - C(w,u) = H(T[U]) for every u in U.  Random T[U], random labels;
             C by brute force over all n! permutations (n <= 10) or by DFS (n = 11, 12); H(T[U]) by a separate DFS of
             the sub-tournament.  Also the decomposition used in the proof, H(T) = H(T[U]) + sum_X H(T[X]) H(T[U-X]).
  reduction  U- nonempty: the exact decomposition behind the remark that 'the swap involution reduces it to a comparison
             of paths ending in U- with paths starting in U-':
               H = N1 + N2,   C(u,w) - C(w,u) = N1 + sum over type-N pairs (A, B') of ([u in A] - [u in B']),
             N1 = sum_{X subset U, X != U} E(X) S+(U-X)   (paths Q1 o Q2 w; E(X) = # paths of T[X] ending in U-, E(empty) = 1;
                                                          S+(Y) = # paths of T[Y] starting in U+),
             type-N pairs: A a path of T[X] starting in U-, B' a path of T[U-X] that is empty or starts in U+.
  random     the empirical claims at every out-degree-1 vertex of random strong tournaments, n = 6..12, both with a
             planted out-degree-1 vertex and unplanted (random strong tournaments that happen to have one).
  near n k   the same over all strong 'transitive (higher beats lower) + k reversed arcs lo -> hi, hi - lo >= 2'.
  exh8       EXHAUSTIVE at n = 8: w = 7 with out-neighbour o = 0, T - w any of the 2^21 labelled tournaments on {0..6};
             every (T, w) with outdeg(w) = 1 on 8 vertices is isomorphic to one of these.
Claims checked: (a) a witness (p(u before w) > 1/2) lies in U+;  (b) if U- is nonempty, the average of 2p - 1 over U-
is > 0 (hence U- has a witness);  (c) the average over U+ can be negative (counted separately: < 0 and exactly = 0).
Also: U- empty <=> indeg(o) = 1;  U+ is never empty in a strong tournament.
"""
import sys, time, itertools
import numpy as np
import numba as nb
from pa_core import (brute_force_perms, enumerate_paths_dfs, count_paths_dfs, strong_landau, strong_warshall,
                     random_tournament)

NR = 14
# res fields
R_PAIRS, R_UMEMPTY, R_INDO1, R_MISM, R_NOWPLUS, R_UPEMPTY, R_UMNONEMPTY, R_AVGUM_LE0, R_NOWMINUS, R_AVGUP_LT0, \
    R_AVGUP_EQ0, R_IDFAIL, R_STRONG, R_SPARE = range(NR)


@nb.njit(cache=False)
def sub_count_nb(n, A, mask):
    m = 0
    for v in range(n):
        if (mask >> v) & 1:
            m += 1
    if m == 0:
        return 1
    idx = np.zeros(m, dtype=np.int64)
    k = 0
    for v in range(n):
        if (mask >> v) & 1:
            idx[k] = v
            k += 1
    B = np.zeros((m, m), dtype=np.int8)
    for i in range(m):
        for j in range(m):
            B[i, j] = A[idx[i], idx[j]]
    return count_paths_dfs(m, B)


@nb.njit(cache=False)
def od1_eval(n, A, C, H, res, mins, wonly):
    """mins[0] = min over cases with U- nonempty of (avg over U- of 2p-1);  mins[1] = min over cases of the best
    margin 2p-1 inside U+;  mins[2] = min over cases of the avg over U+ of 2p-1."""
    for w in range(n):
        if wonly >= 0 and w != wonly:
            continue
        od = 0
        o = -1
        for x in range(n):
            if A[w, x]:
                od += 1
                o = x
        if od != 1:
            continue
        res[R_PAIRS] += 1
        um = 0; up = 0; sp = 0; sm = 0; wp = False; wm = False; bestp = -(1 << 62)
        umask = 0
        for u in range(n):
            if u == w or u == o:
                continue
            umask |= 1 << u
            mg = 2 * C[u, w] - H
            if A[u, o]:
                um += 1
                sm += mg
                if mg > 0:
                    wm = True
            else:
                up += 1
                sp += mg
                if mg > 0:
                    wp = True
                if mg > bestp:
                    bestp = mg
        indo = 0
        for x in range(n):
            if A[x, o]:
                indo += 1
        if um == 0:
            res[R_UMEMPTY] += 1
            HU = sub_count_nb(n, A, umask)
            for u in range(n):
                if (umask >> u) & 1:
                    if 2 * C[u, w] - H != HU:
                        res[R_IDFAIL] += 1
                        break
        if indo == 1:
            res[R_INDO1] += 1
        if (um == 0) != (indo == 1):
            res[R_MISM] += 1
        if not wp:
            res[R_NOWPLUS] += 1
        if up == 0:
            res[R_UPEMPTY] += 1
        else:
            v = bestp / H
            if v < mins[1]:
                mins[1] = v
            v2 = sp / (up * H)
            if v2 < mins[2]:
                mins[2] = v2
        if um > 0:
            res[R_UMNONEMPTY] += 1
            if sm <= 0:
                res[R_AVGUM_LE0] += 1
            if not wm:
                res[R_NOWMINUS] += 1
            v = sm / (um * H)
            if v < mins[0]:
                mins[0] = v
        if up > 0 and sp < 0:
            res[R_AVGUP_LT0] += 1
        if up > 0 and sp == 0:
            res[R_AVGUP_EQ0] += 1


def show(tag, r, mins, t0):
    print(f'  {tag}: out-degree-1 (T, w) pairs {int(r[R_PAIRS])}; U- empty {int(r[R_UMEMPTY])} (identity failures {int(r[R_IDFAIL])}; '
          f'indeg(o) = 1: {int(r[R_INDO1])}, mismatches {int(r[R_MISM])}); NO witness in U+: {int(r[R_NOWPLUS])} '
          f'(U+ empty: {int(r[R_UPEMPTY])}); U- nonempty {int(r[R_UMNONEMPTY])}: avg over U- <= 0 in {int(r[R_AVGUM_LE0])}, '
          f'no witness in U- in {int(r[R_NOWMINUS])}; avg over U+ < 0 in {int(r[R_AVGUP_LT0])}, = 0 exactly in {int(r[R_AVGUP_EQ0])}; '
          f'min avg(2p-1) over U- = {mins[0]:.6f}, min best margin in U+ = {mins[1]:.6f}, min avg over U+ = {mins[2]:.6f}  '
          f'[{time.time() - t0:.1f} s]', flush=True)


# ----------------------------------------------------------------------------------------------------------------
def shape_tournament(n, rng, planted_um_empty=True):
    labels = rng.permutation(n)
    w, o = int(labels[0]), int(labels[1])
    U = [int(x) for x in labels[2:]]
    A = np.zeros((n, n), dtype=np.int8)
    A[w, o] = 1
    for u in U:
        A[u, w] = 1
        A[o, u] = 1
    for i in range(len(U)):
        for j in range(i + 1, len(U)):
            a, b = U[i], U[j]
            if rng.random() < 0.5:
                A[a, b] = 1
            else:
                A[b, a] = 1
    return A, w, o, U


def sub_count(A, S):
    S = list(S)
    if not S:
        return 1
    B = A[np.ix_(S, S)].copy()
    return int(count_paths_dfs(len(S), B))


def mode_identity(rng):
    print('MODE identity: C(u,w) - C(w,u) = H(T[U]) when o beats all of U (U- empty)', flush=True)
    tot = 0; fails = 0; decomp_fail = 0; decomp_tested = 0; strong_fail = 0; witness_fail = 0
    t0 = time.time()
    plan = [(n, 60) for n in range(4, 10)] + [(10, 25), (11, 12), (12, 6)]
    for n, trials in plan:
        nf = 0
        t1 = time.time()
        for _ in range(trials):
            A, w, o, U = shape_tournament(n, rng)
            if not (strong_landau(n, A) and strong_warshall(n, A)):
                strong_fail += 1
            if n <= 10:
                H, C, Nd, L = brute_force_perms(n, A)
                meth = 'brute force n!'
            else:
                C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
                H = enumerate_paths_dfs(n, A, C, Nd, L)
                meth = 'DFS'
            HU = sub_count(A, U)
            tot += 1
            if any(int(C[u, w]) - int(C[w, u]) != HU for u in U):
                fails += 1; nf += 1
            if any(2 * int(C[u, w]) <= int(H) for u in U):
                witness_fail += 1
            if n <= 10:
                decomp_tested += 1
                s = HU
                for r in range(len(U) + 1):
                    for X in itertools.combinations(U, r):
                        Y = [u for u in U if u not in X]
                        s += sub_count(A, X) * sub_count(A, Y)
                if s != int(H):
                    decomp_fail += 1
        print(f'  n = {n}: {trials} random shapes ({meth}), identity failures {nf}  [{time.time() - t1:.1f} s]', flush=True)
    print(f'  TOTAL {tot} tournaments: identity failures {fails}; some u in U not a witness: {witness_fail}; decomposition '
          f'H = H(T[U]) + sum_X H(T[X])H(T[U-X]) checked {decomp_tested}, failures {decomp_fail}; shape tournaments not strong '
          f'(Landau and Warshall): {strong_fail}  [{time.time() - t0:.1f} s]', flush=True)
    # T_n special case: w = 0, o = n-1 in T_n; the value (2^(n-3)+1)/(2^(n-2)+1)
    from fractions import Fraction as Fr
    ok = True
    for n in range(4, 14):
        A = np.zeros((n, n), dtype=np.int8)
        for a in range(n):
            for b in range(a):
                A[a, b] = 1
        A[n - 1, 0] = 0; A[0, n - 1] = 1
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        H = int(enumerate_paths_dfs(n, A, C, Nd, L))
        vals = {Fr(int(C[u, 0]), H) for u in range(1, n - 1)}
        if vals != {Fr(2 ** (n - 3) + 1, 2 ** (n - 2) + 1)} or H != 2 ** (n - 2) + 1:
            ok = False
    print(f'  T_n (n = 4..13), w = 0, o = n-1: every in-neighbour has p = (2^(n-3)+1)/(2^(n-2)+1) and H = 2^(n-2)+1: {ok}', flush=True)


@nb.njit(cache=False)
def first_last_counts(n, A, mask, FL):
    """FL[a, b] = # Hamiltonian paths of T[mask] from a to b (explicit DFS)."""
    for a in range(n):
        for b in range(n):
            FL[a, b] = 0
    m = 0
    for v in range(n):
        if (mask >> v) & 1:
            m += 1
    if m == 0:
        return
    outm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j] and (mask >> i) & 1 and (mask >> j) & 1:
                outm[i] |= 1 << j
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    for s in range(n):
        if not (mask >> s) & 1:
            continue
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == m - 1:
                FL[s, path[depth]] += 1
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            c = cand[depth]
            if c == 0:
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            lb = c & (-c)
            cand[depth] = c & ~lb
            v = 0
            while (lb >> v) != 1:
                v += 1
            depth += 1
            path[depth] = v
            used |= lb
            cand[depth] = outm[v] & ~used


def mode_reduction(rng):
    print("MODE reduction: U- nonempty; H = N1 + N2 and C(u,w) - C(w,u) = N1 + sum over type-N pairs ([u in A] - [u in B2]), B2 = B'", flush=True)
    t0 = time.time()
    tested = 0; bad = 0; badH = 0; neg_typeN = 0
    for n, trials in [(5, 200), (6, 200), (7, 150), (8, 100), (9, 60), (10, 30)]:
        t1 = time.time(); got = 0
        while got < trials:
            A = random_tournament(n, rng)
            w = int(rng.integers(n)); o = int(rng.integers(n - 1)); o = o if o < w else o + 1
            for x in range(n):
                if x != w:
                    A[w, x] = 0; A[x, w] = 1
            A[w, o] = 1; A[o, w] = 0
            U = [u for u in range(n) if u not in (w, o)]
            Um = [u for u in U if A[u, o]]
            Up = [u for u in U if A[o, u]]
            if not Um or not strong_landau(n, A):
                continue
            got += 1
            C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
            H = int(enumerate_paths_dfs(n, A, C, Nd, L))
            FL = np.zeros((n, n), dtype=np.int64)
            E = {}; Sm = {}; Sp = {}; Hs = {}
            for r in range(len(U) + 1):
                for X in itertools.combinations(U, r):
                    mask = sum(1 << x for x in X)
                    first_last_counts(n, A, mask, FL)
                    if r == 0:
                        E[mask] = 1; Sm[mask] = 0; Sp[mask] = 0; Hs[mask] = 1
                    else:
                        E[mask] = int(sum(FL[a, b] for a in X for b in X if b in Um))
                        Sm[mask] = int(sum(FL[a, b] for a in X for b in X if a in Um))
                        Sp[mask] = int(sum(FL[a, b] for a in X for b in X if a in Up))
                        Hs[mask] = int(FL.sum())
            full = sum(1 << u for u in U)
            N1 = sum(E[X] * Sp[full ^ X] for X in E if X != full)
            N2 = sum(Hs[X] * (1 if (full ^ X) == 0 else Sp[full ^ X]) for X in E)
            if N1 + N2 != H:
                badH += 1
            tN = {u: 0 for u in U}
            for X in E:
                if X == 0:
                    continue
                cnt = Sm[X] * (1 if (full ^ X) == 0 else Sp[full ^ X])
                for u in U:
                    tN[u] += cnt if (X >> u) & 1 else -cnt
            for u in U:
                if 2 * int(C[u, w]) - H != N1 + tN[u]:
                    bad += 1
                    break
            if sum(tN[u] for u in Um) < 0:
                neg_typeN += 1
            tested += 1
        print(f'  n = {n}: {trials} tournaments with U- nonempty  [{time.time() - t1:.1f} s]', flush=True)
    print(f'  TOTAL {tested}: H = N1 + N2 failures {badH}; per-u formula failures {bad}; cases where the type-N term summed over U- '
          f'is negative (so the case-1 term N1 must compensate): {neg_typeN}  [{time.time() - t0:.1f} s]', flush=True)


def mode_random(rng):
    print('MODE random: empirical out-degree-1 claims on random strong tournaments', flush=True)
    for planted in (True, False):
        plan = ([(6, 3000), (7, 3000), (8, 3000), (9, 2000), (10, 800), (11, 300), (12, 80)] if planted else
                [(6, 3000), (7, 3000), (8, 2000), (9, 1000), (10, 400)])
        print(f' {"planted out-degree-1 vertex (random w, o, rest random)" if planted else "unplanted (random strong tournaments having an out-degree-1 vertex)"}', flush=True)
        for n, trials in plan:
            t0 = time.time()
            r = np.zeros(NR, dtype=np.int64); mins = np.full(3, 9.0)
            got = 0
            while got < trials:
                A = random_tournament(n, rng)
                if planted:
                    w = int(rng.integers(n)); o = int(rng.integers(n - 1)); o = o if o < w else o + 1
                    for x in range(n):
                        if x != w:
                            A[w, x] = 0; A[x, w] = 1
                    A[w, o] = 1; A[o, w] = 0
                else:
                    if not any(int(A[v].sum()) == 1 for v in range(n)):
                        continue
                if not strong_landau(n, A):
                    continue
                got += 1
                C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
                H = enumerate_paths_dfs(n, A, C, Nd, L)
                od1_eval(n, A, C, H, r, mins, -1)
            show(f'n = {n}, {trials} tournaments', r, mins, t0)


@nb.njit(parallel=True, cache=False)
def near_sweep(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]
    out = np.zeros((nthreads, NR), dtype=np.int64)
    mins = np.full((nthreads, 3), 9.0)
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        A = np.zeros((n, n), dtype=np.int8)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        for idx in range(t * chunk, min(m, (t + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 1 if b < a else 0
            for e in range(k):
                lo = arcsets[idx, e, 0]; hi = arcsets[idx, e, 1]
                A[hi, lo] = 0; A[lo, hi] = 1
            if not strong_landau(n, A):
                continue
            out[t, R_STRONG] += 1
            has = False
            for w in range(n):
                s = 0
                for x in range(n):
                    s += A[w, x]
                if s == 1:
                    has = True
            if not has:
                continue
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            od1_eval(n, A, C, H, out[t], mins[t], -1)
    res = np.zeros(NR, dtype=np.int64)
    mm = np.full(3, 9.0)
    for t in range(nthreads):
        for i in range(NR):
            res[i] += out[t, i]
        for i in range(3):
            if mins[t, i] < mm[i]:
                mm[i] = mins[t, i]
    return res, mm


def mode_near(n, k):
    t0 = time.time()
    cand = [(lo, hi) for lo in range(n) for hi in range(lo + 2, n)]
    combos = np.array(list(itertools.combinations(cand, k)), dtype=np.int64).reshape(-1, k, 2)
    r, mins = near_sweep(n, combos, 12)
    show(f'near-transitive n = {n}, k = {k} ({len(combos)} arc sets, {int(r[R_STRONG])} strong)', r, mins, t0)


@nb.njit(parallel=True, cache=False)
def exh8_sweep(nchunk):
    n = 8
    E = 21
    pu = np.zeros(E, dtype=np.int64); pv = np.zeros(E, dtype=np.int64)
    e = 0
    for i in range(7):
        for j in range(i + 1, 7):
            pu[e] = i; pv[e] = j
            e += 1
    total = 1 << E
    out = np.zeros((nchunk, NR), dtype=np.int64)
    mins = np.full((nchunk, 3), 9.0)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        A = np.zeros((n, n), dtype=np.int8)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 0
            for q in range(E):
                if (code >> q) & 1:
                    A[pu[q], pv[q]] = 1
                else:
                    A[pv[q], pu[q]] = 1
            A[7, 0] = 1
            for u in range(1, 7):
                A[u, 7] = 1
            if not strong_landau(n, A):
                continue
            out[c, R_STRONG] += 1
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            od1_eval(n, A, C, H, out[c], mins[c], 7)
    res = np.zeros(NR, dtype=np.int64)
    mm = np.full(3, 9.0)
    for t in range(nchunk):
        for i in range(NR):
            res[i] += out[t, i]
        for i in range(3):
            if mins[t, i] < mm[i]:
                mm[i] = mins[t, i]
    return res, mm


@nb.njit(cache=False)
def typeN_direct(n, A, w, o, out3):
    """Direct per-path classification at an out-degree-1 vertex w (out-neighbour o), explicit DFS:
    out3[0] = N1 = # paths ending at w;  out3[1] = S_N = sum over type-N paths (w not last, first vertex in U-) of
    (# U- before w) - (# U- after w);  out3[2] = sum over type-S paths (w not last, w first or first vertex in U+) of the
    same quantity (must be 0 by the swap involution);  out3[3] = sum_{u in U-} (C(u,w) - C(w,u)) computed per path."""
    umask = 0
    for u in range(n):
        if u != w and u != o and A[u, o]:
            umask |= 1 << u
    outm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                outm[i] |= 1 << j
    for i in range(4):
        out3[i] = 0
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    for s in range(n):
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == n - 1:
                before = 0; after = 0; seenw = False
                for t in range(n):
                    y = path[t]
                    if y == w:
                        seenw = True
                    elif (umask >> y) & 1:
                        if seenw:
                            after += 1
                        else:
                            before += 1
                q = before - after
                out3[3] += q
                if path[n - 1] == w:
                    out3[0] += 1
                elif (umask >> path[0]) & 1:
                    out3[1] += q
                else:
                    out3[2] += q
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            c = cand[depth]
            if c == 0:
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            lb = c & (-c)
            cand[depth] = c & ~lb
            v = 0
            while (lb >> v) != 1:
                v += 1
            depth += 1
            path[depth] = v
            used |= lb
            cand[depth] = outm[v] & ~used


@nb.njit(parallel=True, cache=False)
def typeN_exh8(nchunk):
    """res: [cases with U- nonempty, S_N < 0, type-S sum != 0, identity sum_{U-} != |U-| N1 + S_N, N1 == 0];
    mn: min of S_N / (|U-| N1)."""
    n = 8
    E = 21
    pu = np.zeros(E, dtype=np.int64); pv = np.zeros(E, dtype=np.int64)
    e = 0
    for i in range(7):
        for j in range(i + 1, 7):
            pu[e] = i; pv[e] = j
            e += 1
    total = 1 << E
    out = np.zeros((nchunk, 5), dtype=np.int64)
    mn = np.full(nchunk, 9.0)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        A = np.zeros((n, n), dtype=np.int8)
        o3 = np.zeros(4, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 0
            for q in range(E):
                if (code >> q) & 1:
                    A[pu[q], pv[q]] = 1
                else:
                    A[pv[q], pu[q]] = 1
            A[7, 0] = 1
            for u in range(1, 7):
                A[u, 7] = 1
            if not strong_landau(n, A):
                continue
            m = 0
            for u in range(1, 7):
                if A[u, 0]:
                    m += 1
            if m == 0:
                continue
            typeN_direct(n, A, 7, 0, o3)
            out[c, 0] += 1
            if o3[1] < 0:
                out[c, 1] += 1
            if o3[2] != 0:
                out[c, 2] += 1
            if o3[3] != m * o3[0] + o3[1]:
                out[c, 3] += 1
            if o3[0] == 0:
                out[c, 4] += 1
            else:
                v = o3[1] / (m * o3[0])
                if v < mn[c]:
                    mn[c] = v
    res = np.zeros(5, dtype=np.int64)
    for t in range(nchunk):
        for i in range(5):
            res[i] += out[t, i]
    return res, mn.min()


def typeN_check_one(n, A, res, mn):
    for w in range(n):
        if int(A[w].sum()) != 1:
            continue
        o = int(np.argmax(A[w]))
        m = sum(1 for u in range(n) if u not in (w, o) and A[u, o])
        if m == 0:
            continue
        o3 = np.zeros(4, dtype=np.int64)
        typeN_direct(n, A, w, o, o3)
        res[0] += 1
        if o3[1] < 0:
            res[1] += 1
        if o3[2] != 0:
            res[2] += 1
        if o3[3] != m * o3[0] + o3[1]:
            res[3] += 1
        if o3[0] == 0:
            res[4] += 1
        else:
            mn[0] = min(mn[0], o3[1] / (m * o3[0]))


def mode_typeN():
    print('MODE typeN (side observation): S_N = sum over type-N paths of (#U- before w - #U- after w); if S_N >= 0 always,'
          ' the U- conjecture would follow from N1 >= 1', flush=True)
    t0 = time.time()
    r, mn = typeN_exh8(240)
    print(f'  n = 8 exhaustive (w = 7, o = 0), U- nonempty: {int(r[0])} cases; S_N < 0 in {int(r[1])}; type-S sum != 0 in {int(r[2])}; '
          f'identity sum_U- = |U-| N1 + S_N fails in {int(r[3])}; N1 = 0 in {int(r[4])}; min S_N/(|U-| N1) = {mn:.6f}  [{time.time() - t0:.1f} s]', flush=True)
    rng = np.random.default_rng(17102026)
    for n, trials in [(9, 1500), (10, 600), (11, 200), (12, 60)]:
        t0 = time.time(); res = [0] * 5; mn = [9.0]; got = 0
        while got < trials:
            A = random_tournament(n, rng)
            w = int(rng.integers(n)); o = int(rng.integers(n - 1)); o = o if o < w else o + 1
            for x in range(n):
                if x != w:
                    A[w, x] = 0; A[x, w] = 1
            A[w, o] = 1; A[o, w] = 0
            if not strong_landau(n, A):
                continue
            got += 1
            typeN_check_one(n, A, res, mn)
        print(f'  n = {n}, {trials} random planted: U- nonempty {res[0]}; S_N < 0 in {res[1]}; type-S != 0 in {res[2]}; identity fails {res[3]}; '
              f'N1 = 0 in {res[4]}; min S_N/(|U-| N1) = {mn[0]:.6f}  [{time.time() - t0:.1f} s]', flush=True)
    for n, k in [(10, 2), (12, 2), (14, 2), (9, 3), (10, 3), (11, 3), (12, 3), (9, 4), (10, 4)]:
        t0 = time.time(); res = [0] * 5; mn = [9.0]
        cand = [(lo, hi) for lo in range(n) for hi in range(lo + 2, n)]
        for arcs in itertools.combinations(cand, k):
            A = np.zeros((n, n), dtype=np.int8)
            for a in range(n):
                for b in range(a):
                    A[a, b] = 1
            for lo, hi in arcs:
                A[hi, lo] = 0; A[lo, hi] = 1
            if not any(int(A[v].sum()) == 1 for v in range(n)):
                continue
            if not strong_landau(n, A):
                continue
            typeN_check_one(n, A, res, mn)
        print(f'  near-transitive n = {n}, k = {k}: U- nonempty {res[0]}; S_N < 0 in {res[1]}; type-S != 0 in {res[2]}; identity fails {res[3]}; '
              f'N1 = 0 in {res[4]}; min S_N/(|U-| N1) = {mn[0]:.6f}  [{time.time() - t0:.1f} s]', flush=True)


def mode_exh8():
    print('MODE exh8: all 2^21 tournaments on 8 vertices with w = 7 of out-degree 1 (o = 0)', flush=True)
    t0 = time.time()
    r, mins = exh8_sweep(240)
    print(f'  configurations: 2097152; strong: {int(r[R_STRONG])}', flush=True)
    show('n = 8 exhaustive (w = 7 only)', r, mins, t0)


if __name__ == '__main__':
    # usage: identity | reduction | random | near n k | exh8
    rng = np.random.default_rng(14102026)
    mode = sys.argv[1]
    if mode == 'identity':
        mode_identity(rng)
    elif mode == 'reduction':
        mode_reduction(np.random.default_rng(15102026))
    elif mode == 'random':
        mode_random(np.random.default_rng(16102026))
    elif mode == 'near':
        print('MODE near: all strong transitive + k reversed arcs (hi - lo >= 2)', flush=True)
        a = [int(x) for x in sys.argv[2:]]
        for n, k in zip(a[0::2], a[1::2]):
            mode_near(n, k)
    elif mode == 'exh8':
        mode_exh8()
    elif mode == 'typeN':
        mode_typeN()
