#!/usr/bin/env python3
"""
verify_B_tower_kps1.py — kind-pasteur-2026-06-09-S1 ADVERSARIAL VERIFICATION of branch B-tower.

Independently recomputes the central claims of drt_mersenne_tower_kps1 /
drt_tower31_paley_kps1 / drt_tower_followup_kps1 using DIFFERENT constructions
and DIFFERENT algorithms:

 V1. Tower matrix S_N two ways: (a) the doubling recursion (definition), (b) a CLOSED FORM
     derived by hand from the recursion: for i != j, with v = lowest set bit of i XOR j,
       S[i,j] = (-1)^( popcount((i AND j) >> (v+1)) + bit_v(i) );   S[i,i] = +1.
     Skew-Hadamard re-checked directly; row-0 all-ones follows from the closed form
     (i=0 => i AND j = 0 and bit_v(0)=0 => +1 for every j).
 V2. Cores: regularity; c3 by BRUTE-FORCE triple orientation count (not the score formula);
     doubly-regular by brute-force pairwise |N+ cap N+| (popcounts on out-masks).
 V3. H by three independent algorithms: permutation brute force (n<=8), inclusion-exclusion
     walk counts (IE: #HamPaths = sum_S (-1)^(n-|S|) walks_{n-1}(S), vectorized over all masks),
     and the library subset-DP H_count. Anchors: H(C3)=3, H(regT5)=15, H(Paley7)=189, random n=8.
     Then H(T15) by IE vs DP.
 V4. T15 link-H spectrum by PERMUTATION BRUTE FORCE (each link is a 7-tournament, 5040 perms).
 V5. T31 link-H spectrum (all 31 vertices) by the IE method (no subset-DP involved).
 V6. Triple/quad distributions of T31 vs Paley31 via bitmask-set code (not matrix products).
 V7. My own backtracking iso/aut searcher (checks BOTH arc directions, natural vertex order):
     |Aut(T7/T15/T31/Paley31)| plain-exhaustive, |Aut(T63)| with independently recomputed
     link-triple-dist colors; element orders, orbits, fixed points; B_0(T31) iso T15,
     B-_0(T31) iso T15, B_8(T31) !iso T15, B_0(Paley31) !iso T15, B_0(T63) iso T31;
     self-converse probes T15/T31 (exhaustive), T63 (color-multiset).
 V8. Circulant DRT enumerations: Z_15 (128 sets, plain loops) and Z_31 (32768 sets, np.roll).
 V9. nnz(S_tower - Walsh) for N=2..64 with Walsh built as (-1)^popcount(i AND j).

Output tee'd to 05-knowledge/results/verify_B_tower_kps1.out
"""
import sys, time, itertools
from math import comb, gcd
from collections import Counter
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import H_count  # theirs — cross-checked against perm/IE below

OUT = open('05-knowledge/results/verify_B_tower_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

try:
    (1).bit_count
    def pc(x): return x.bit_count()
except AttributeError:
    def pc(x): return bin(x).count('1')

# ---------------- V1: tower two ways ----------------
def S_block(N):
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < N:
        n = S.shape[0]; I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    return S

def S_closed(N):
    S = np.ones((N, N), dtype=np.int64)
    for i in range(N):
        for j in range(N):
            if i != j:
                x = i ^ j
                v = (x & -x).bit_length() - 1
                e = pc((i & j) >> (v + 1)) + ((i >> v) & 1)
                S[i, j] = 1 - 2 * (e & 1)
    return S

def is_SH(S):
    N = S.shape[0]
    return bool((np.abs(S) == 1).all()
                and np.array_equal(S + S.T, 2 * np.eye(N, dtype=np.int64))
                and np.array_equal(S @ S.T, N * np.eye(N, dtype=np.int64)))

def core_of(S):
    n = S.shape[0] - 1
    M = S[1:, 1:] - np.eye(n, dtype=np.int64)
    return (M > 0).astype(np.int64)

# ---------------- V2: brute structure ----------------
def c3_brute(A):
    n = A.shape[0]; L = A.tolist(); c = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if L[i][j] + L[i][k] == 1 and L[j][i] + L[j][k] == 1 \
                        and L[k][i] + L[k][j] == 1:
                    c += 1
    return c

def outmasks(A):
    n = A.shape[0]; L = A.tolist()
    return [sum(1 << j for j in range(n) if L[i][j]) for i in range(n)]

def drt_brute(A):
    n = A.shape[0]
    om = outmasks(A)
    vals = set()
    for u in range(n):
        for v in range(u + 1, n):
            vals.add(pc(om[u] & om[v]))
    reg = len(set(int(x) for x in A.sum(axis=1))) == 1
    return reg, sorted(vals)

def tdist(A):
    n = A.shape[0]; om = outmasks(A); cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            muv = om[u] & om[v]
            for t in range(v + 1, n):
                cnt[pc(muv & om[t])] += 1
    return cnt

def qdist(A):
    n = A.shape[0]; om = outmasks(A); cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            muv = om[u] & om[v]
            for x in range(v + 1, n):
                muvx = muv & om[x]
                for t in range(x + 1, n):
                    cnt[pc(muvx & om[t])] += 1
    return cnt

# ---------------- V3: H three ways ----------------
def H_perm(A):
    n = A.shape[0]; L = A.tolist()
    return sum(1 for p in itertools.permutations(range(n))
               if all(L[p[i]][p[i + 1]] for i in range(n - 1)))

def H_IE(A):
    """#HamPaths = sum over vertex subsets S of (-1)^(n-|S|) * (# walks of length n-1 in T[S])."""
    n = A.shape[0]
    NM = 1 << n
    masks = np.arange(NM, dtype=np.int64)
    Mask = ((masks[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
    V = Mask.copy()
    An = A.astype(np.int64)
    for _ in range(n - 1):
        V = (V @ An) * Mask
    tot = V.sum(axis=1)
    sign = 1 - 2 * ((n - Mask.sum(axis=1)) & 1)
    return int((sign * tot).sum())

def link(A, v, direction='out'):
    idx = np.flatnonzero(A[v] if direction == 'out' else A[:, v])
    return A[np.ix_(idx, idx)]

# ---------------- V7: my searcher (checks BOTH directions) ----------------
def my_iso(A, B, colA=None, colB=None, find_all=True, node_cap=None):
    n = A.shape[0]
    res = {'count': 0, 'perms': [], 'nodes': 0, 'complete': True}
    if B.shape[0] != n:
        return res
    LA = A.tolist()
    Bout = outmasks(B)
    LB = B.tolist()
    Bin = [sum(1 << j for j in range(n) if LB[j][i]) for i in range(n)]
    if colA is None: colA = [0] * n
    if colB is None: colB = [0] * n
    if sorted(colA) != sorted(colB):
        return res
    img = [-1] * n
    def rec(k, used):
        if k == n:
            res['count'] += 1
            res['perms'].append(img.copy())
            return not find_all
        need_out = 0; need_in = 0
        Ak = LA[k]
        for j in range(k):
            if Ak[j]:
                need_out |= 1 << img[j]
            elif LA[j][k]:
                need_in |= 1 << img[j]
        ck = colA[k]
        for t in range(n):
            bt_ = 1 << t
            if used & bt_:
                continue
            if colB[t] != ck:
                continue
            res['nodes'] += 1
            if node_cap is not None and res['nodes'] > node_cap:
                res['complete'] = False
                return True
            if (Bout[t] & used) == need_out and (Bin[t] & used) == need_in:
                img[k] = t
                if rec(k + 1, used | bt_):
                    return True
        return False
    rec(0, 0)
    return res

def perm_order(p):
    n = len(p); seen = [False] * n; o = 1
    for v0 in range(n):
        if seen[v0]:
            continue
        l = 0; x = v0
        while not seen[x]:
            seen[x] = True; x = p[x]; l += 1
        o = o * l // gcd(o, l)
    return o

def orbits(n, perms):
    par = list(range(n))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for p in perms:
        for v in range(n):
            a, b = f(v), f(p[v])
            if a != b:
                par[a] = b
    g = {}
    for v in range(n):
        g.setdefault(f(v), []).append(v)
    return sorted(g.values(), key=lambda x: (-len(x), x))

def circ(n, Sset):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in Sset:
                A[i, j] = 1
    return A

def main():
    t0_all = time.time()
    PASS = []; FAIL = []
    def verdict(name, ok, detail=''):
        (PASS if ok else FAIL).append(name)
        w(f'[{"PASS" if ok else "FAIL"}] {name}' + (f'  — {detail}' if detail else ''))

    w('=== verify_B_tower_kps1 — adversarial re-verification of branch B-tower ===')
    w('')

    # ---- V1 ----
    w('--- V1: tower construction, two independent builds ---')
    for N in (2, 4, 8, 16, 32, 64):
        Sb = S_block(N); Sc = S_closed(N)
        same = np.array_equal(Sb, Sc)
        sh = is_SH(Sc)
        row0 = bool((Sc[0] == 1).all())
        verdict(f'S_{N}: closed form == recursion', same)
        verdict(f'S_{N}: skew-Hadamard (my check)', sh)
        verdict(f'S_{N}: row 0 all ones', row0)
    S16, S32, S64 = S_closed(16), S_closed(32), S_closed(64)
    T15, T31, T63 = core_of(S16), core_of(S32), core_of(S64)
    T7 = core_of(S_closed(8))
    w('')

    # ---- V2 ----
    w('--- V2: core structure by brute force ---')
    for name, T, c3_claim, lam_claim in (('T15', T15, 140, 3), ('T31', T31, 1240, 7),
                                         ('T63', T63, 10416, 15)):
        sc = sorted(set(int(x) for x in T.sum(axis=1)))
        c3 = c3_brute(T)
        reg, vals = drt_brute(T)
        verdict(f'{name}: regular with score {(T.shape[0]-1)//2}', sc == [(T.shape[0] - 1) // 2],
                f'score set {sc}')
        verdict(f'{name}: c3 brute-force == {c3_claim}', c3 == c3_claim, f'c3={c3}')
        verdict(f'{name}: doubly regular, |N+capN+| == {lam_claim} const',
                reg and vals == [lam_claim], f'pairwise values {vals}')
    M15 = T15 - T15.T
    J15 = np.ones((15, 15), dtype=np.int64)
    verdict('T15: M^2 == J - 15I exactly',
            np.array_equal(M15 @ M15, J15 - 15 * np.eye(15, dtype=np.int64)))
    w('')

    # ---- V3 ----
    w('--- V3: H by three independent algorithms ---')
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    regT5 = circ(5, {1, 2})
    P7 = circ(7, {1, 2, 4})
    rng = np.random.default_rng(20260609)
    R8 = np.zeros((8, 8), dtype=np.int64)
    for i in range(8):
        for j in range(i + 1, 8):
            if rng.integers(2):
                R8[i, j] = 1
            else:
                R8[j, i] = 1
    for name, T, expect in (('C3', C3, 3), ('regT5', regT5, 15), ('Paley7', P7, 189),
                            ('random8', R8, None)):
        hp, hie, hdp = H_perm(T), H_IE(T), H_count(T)
        ok = (hp == hie == hdp) and (expect is None or hp == expect)
        verdict(f'H anchor {name}: perm={hp} IE={hie} DP={hdp}'
                + (f' expect={expect}' if expect is not None else ''), ok)
    t0 = time.time()
    h_ie_15 = H_IE(T15)
    t_ie = time.time() - t0
    h_dp_15 = H_count(T15)
    verdict('H(T15) == 198335025 by IE (independent algorithm)',
            h_ie_15 == 198335025, f'IE={h_ie_15} ({t_ie:.1f}s) DP={h_dp_15}')
    verdict('H(T15) IE == DP agreement', h_ie_15 == h_dp_15)
    w(f'    H(T15)/(15!/2^14) = {h_ie_15 / (1307674368000 / 16384):.4f}  (claim 2.4850)')
    w('')

    # ---- V4 ----
    w('--- V4: T15 link-H spectrum by PERMUTATION brute force ---')
    linkH15 = {}
    for v in range(15):
        B = link(T15, v)
        linkH15.setdefault(H_perm(B), []).append(v)
    w(f'    spectrum: {sorted(((h, vs) for h, vs in linkH15.items()), reverse=True)}')
    verdict('T15 link-H == {189: [0..7], 171: [8..14]}',
            linkH15 == {189: list(range(8)), 171: list(range(8, 15))})
    b0 = link(T15, 0)
    verdict('B_0(T15) iso Paley7 (my searcher)', my_iso(b0, P7, find_all=False)['count'] == 1)
    b8 = link(T15, 8)
    r = my_iso(b8, P7, find_all=False)
    verdict('B_8(T15) NOT iso Paley7 (exhaustive)', r['count'] == 0 and r['complete'],
            f'nodes={r["nodes"]}')
    w('')

    # ---- V5 ----
    w('--- V5: T31 link-H spectrum (all 31 links) by IE method ---')
    t0 = time.time()
    linkH31 = {}
    for v in range(31):
        B = link(T31, v)
        linkH31.setdefault(H_IE(B), []).append(v)
    w(f'    ({time.time()-t0:.1f}s)')
    for h, vs in sorted(linkH31.items(), reverse=True):
        w(f'    H = {h}: vertices {vs} (count {len(vs)})')
    claim31 = {198335025: [0, 1, 2, 3, 4, 5, 6, 7, 15],
               197568027: [24, 25, 26, 27, 28, 29, 30],
               197147697: [8, 9, 10, 11, 12, 13, 14],
               196908495: [23],
               196179375: [16, 17, 18, 19, 20, 21, 22]}
    verdict('T31 link-H spectrum matches claimed 5-class split', linkH31 == claim31)
    w('')

    # ---- V6 ----
    w('--- V6: T31 vs Paley31 invariants (bitmask-set code) ---')
    QR31 = sorted({(x * x) % 31 for x in range(1, 31)})
    P31 = circ(31, set(QR31))
    regP, valsP = drt_brute(P31)
    verdict('Paley31 doubly regular lam=7 (sanity)', regP and valsP == [7])
    verdict('c3(Paley31) brute == 1240', c3_brute(P31) == 1240, f'c3={c3_brute(P31)}')
    tdT, tdP = tdist(T31), tdist(P31)
    w(f'    T31   triple dist: {dict(sorted(tdT.items()))}')
    w(f'    Paley triple dist: {dict(sorted(tdP.items()))}')
    claim_tdT = {0: 28, 1: 84, 2: 252, 3: 3192, 4: 812, 5: 84, 6: 28, 7: 15}
    claim_tdP = {2: 930, 3: 2015, 4: 1550}
    verdict('T31 triple dist matches claim', dict(tdT) == claim_tdT)
    verdict('Paley31 triple dist matches claim', dict(tdP) == claim_tdP)
    verdict('T31 NOT iso Paley31 (triple dists differ — invariant proof)', tdT != tdP)
    qdT, qdP = qdist(T31), qdist(P31)
    verdict('T31 quad dist matches claim',
            dict(qdT) == {0: 3248, 1: 16296, 2: 9744, 3: 2177}, f'{dict(sorted(qdT.items()))}')
    verdict('Paley31 quad dist matches claim',
            dict(qdP) == {0: 3875, 1: 14415, 2: 11625, 3: 1550}, f'{dict(sorted(qdP.items()))}')
    w('')

    # ---- V7 ----
    w('--- V7: automorphism groups / isomorphisms with my searcher ---')
    for name, T, exp_orb in (('T7', T7, [7]), ('T15', T15, [7, 7, 1])):
        r = my_iso(T, T)
        oo = Counter(perm_order(p) for p in r['perms'])
        ob = [len(g) for g in orbits(T.shape[0], r['perms'])]
        verdict(f'|Aut({name})| == 21, orders {{1:1,3:14,7:6}}, orbits {exp_orb}',
                r['count'] == 21 and dict(oo) == {1: 1, 3: 14, 7: 6} and ob == exp_orb,
                f'count={r["count"]} orders={dict(sorted(oo.items()))} orbits={ob} nodes={r["nodes"]}')
    t0 = time.time()
    rA31 = my_iso(T31, T31)
    oo31 = Counter(perm_order(p) for p in rA31['perms'])
    orb31 = orbits(31, rA31['perms'])
    fix31 = sorted(g[0] for g in orb31 if len(g) == 1)
    verdict('|Aut(T31)| == 21 (plain exhaustive), F_21 orders, fixed {7,15,23}',
            rA31['count'] == 21 and dict(oo31) == {1: 1, 3: 14, 7: 6} and fix31 == [7, 15, 23],
            f'count={rA31["count"]} orders={dict(sorted(oo31.items()))} '
            f'orbits={[len(g) for g in orb31]} fixed={fix31} nodes={rA31["nodes"]} ({time.time()-t0:.1f}s)')
    t0 = time.time()
    rP31 = my_iso(P31, P31)
    verdict('|Aut(Paley31)| == 465 (plain exhaustive)', rP31['count'] == 465,
            f'count={rP31["count"]} nodes={rP31["nodes"]} ({time.time()-t0:.1f}s)')
    # explicit iso attempt T31 -> P31 with my searcher, exhaustive
    t0 = time.time()
    rX = my_iso(T31, P31, find_all=False, node_cap=80_000_000)
    verdict('T31 -> Paley31: my exhaustive search finds NO map',
            rX['count'] == 0 and rX['complete'], f'nodes={rX["nodes"]} ({time.time()-t0:.1f}s)')
    # links
    B0 = link(T31, 0)
    verdict('B_0(T31) literally equals T15 as labeled matrix?', np.array_equal(B0, T15),
            'informational')
    verdict('B_0(T31) iso T15', my_iso(B0, T15, find_all=False)['count'] == 1)
    Bin0 = link(T31, 0, 'in')
    verdict('B-_0(T31) (in-link) iso T15', my_iso(Bin0, T15, find_all=False)['count'] == 1)
    B8 = link(T31, 8)
    r8 = my_iso(B8, T15, find_all=False)
    reg8, vals8 = drt_brute(B8)
    a8 = my_iso(B8, B8)
    verdict('B_8(T31): NOT iso T15, NOT DRT, |Aut|=3, H=197147697',
            r8['count'] == 0 and r8['complete'] and not (reg8 and len(vals8) == 1)
            and a8['count'] == 3 and H_IE(B8) == 197147697,
            f'isoT15={r8["count"]} pairwise={vals8} |Aut|={a8["count"]} H={H_IE(B8)}')
    BP = link(P31, 0)
    rBP = my_iso(BP, T15, find_all=False)
    aBP = my_iso(BP, BP)
    regBP, valsBP = drt_brute(BP)
    verdict('B_0(Paley31): NOT iso T15, NOT DRT, |Aut|=15, H=196922055',
            rBP['count'] == 0 and rBP['complete'] and not (regBP and len(valsBP) == 1)
            and aBP['count'] == 15 and H_IE(BP) == 196922055,
            f'isoT15={rBP["count"]} pairwise={valsBP} |Aut|={aBP["count"]} H={H_IE(BP)}')
    # self-converse
    rsc15 = my_iso(T15, np.ascontiguousarray(T15.T), find_all=False)
    verdict('T15 NOT self-converse (exhaustive)', rsc15['count'] == 0 and rsc15['complete'],
            f'nodes={rsc15["nodes"]}')
    t0 = time.time()
    rsc31 = my_iso(T31, np.ascontiguousarray(T31.T), find_all=False, node_cap=80_000_000)
    verdict('T31 NOT self-converse (exhaustive)', rsc31['count'] == 0 and rsc31['complete'],
            f'nodes={rsc31["nodes"]} ({time.time()-t0:.1f}s)')
    # T63: colors from link triple dists (recomputed independently)
    t0 = time.time()
    sigs = []
    for v in range(63):
        sigs.append(tuple(sorted(tdist(link(T63, v)).items())))
    sigmap = {}
    col63 = []
    for s_ in sigs:
        sigmap.setdefault(s_, len(sigmap))
        col63.append(sigmap[s_])
    csz = sorted(Counter(col63).values(), reverse=True)
    verdict('T63 link-coloring: 12 classes sizes [10,7x7,1x4]',
            len(sigmap) == 12 and csz == [10, 7, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1],
            f'{len(sigmap)} classes {csz} ({time.time()-t0:.1f}s)')
    t0 = time.time()
    rA63 = my_iso(T63, T63, colA=col63, colB=col63)
    oo63 = Counter(perm_order(p) for p in rA63['perms'])
    orb63 = orbits(63, rA63['perms'])
    fix63 = sorted(g[0] for g in orb63 if len(g) == 1)
    verdict('|Aut(T63)| == 21, F_21 orders, fixed {7,15,...,55}',
            rA63['count'] == 21 and dict(oo63) == {1: 1, 3: 14, 7: 6}
            and fix63 == [7, 15, 23, 31, 39, 47, 55],
            f'count={rA63["count"]} orders={dict(sorted(oo63.items()))} '
            f'orbits={[len(g) for g in orb63]} fixed={fix63} nodes={rA63["nodes"]} ({time.time()-t0:.1f}s)')
    # T63 self-converse by color multisets
    sigs_op = []
    T63op = np.ascontiguousarray(T63.T)
    for v in range(63):
        sigs_op.append(tuple(sorted(tdist(link(T63op, v)).items())))
    verdict('T63 NOT self-converse (link triple-dist multisets differ)',
            Counter(sigs) != Counter(sigs_op))
    # B_0(T63) iso T31
    B0_63 = link(T63, 0)
    verdict('B_0(T63) literally equals T31 as labeled matrix?', np.array_equal(B0_63, T31),
            'informational')
    rL = my_iso(B0_63, T31, find_all=False, node_cap=50_000_000)
    verdict('B_0(T63) iso T31', rL['count'] == 1, f'nodes={rL["nodes"]}')
    w('')

    # ---- V8 ----
    w('--- V8: circulant DRT enumerations ---')
    # Z_15: all 128 antisymmetric S-sets, plain loops
    hits15 = 0
    for bits in range(128):
        Sset = set()
        for kk in range(7):
            d = kk + 1
            Sset.add(d if (bits >> kk) & 1 else 15 - d)
        Ac = circ(15, Sset)
        reg, vals = drt_brute(Ac)
        if reg and len(vals) == 1:
            hits15 += 1
    verdict('Z_15: ZERO circulant DRT(15) among all 128 antisymmetric S-sets', hits15 == 0,
            f'hits={hits15}')
    # Z_31: vectorized with np.roll
    t0 = time.time()
    NSETS = 1 << 15
    choices = np.arange(NSETS, dtype=np.int64)
    Svec = np.zeros((NSETS, 31), dtype=np.int8)
    for kk in range(15):
        d = kk + 1
        bit = ((choices >> kk) & 1).astype(np.int8)
        Svec[:, d] = bit
        Svec[:, 31 - d] = 1 - bit
    corr_min = np.full(NSETS, 99, dtype=np.int32)
    corr_max = np.full(NSETS, -1, dtype=np.int32)
    for d in range(1, 31):
        rolled = np.roll(Svec, d, axis=1)
        cd = (Svec.astype(np.int32) * rolled.astype(np.int32)).sum(axis=1)
        corr_min = np.minimum(corr_min, cd)
        corr_max = np.maximum(corr_max, cd)
    const = np.flatnonzero(corr_min == corr_max)
    QRfs = frozenset(QR31)
    NQRfs = frozenset(range(1, 31)) - QRfs
    found = [frozenset(int(x) for x in np.flatnonzero(Svec[i])) for i in const]
    verdict('Z_31: exactly 2 circulant DRT(31) S-sets = QR and NQR',
            len(found) == 2 and set(found) == {QRfs, NQRfs},
            f'{len(found)} hits ({time.time()-t0:.1f}s)')
    w('')

    # ---- V9 ----
    w('--- V9: nnz(S_tower - Walsh) ---')
    ok9 = True
    for k in range(1, 7):
        N = 1 << k
        S = S_closed(N)
        W = np.array([[1 - 2 * (pc(i & j) & 1) for j in range(N)] for i in range(N)],
                     dtype=np.int64)
        nnz = int((S != W).sum())
        w(f'    N={N}: nnz = {nnz}  (N^2/2 = {N * N // 2})')
        ok9 &= (nnz == N * N // 2)
    verdict('nnz(S_tower - Walsh) == N^2/2 at N=2..64', ok9)
    H16kron = np.array([[1]], dtype=np.int64)
    H2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    for _ in range(4):
        H16kron = np.kron(H2, H16kron)
    W16 = np.array([[1 - 2 * (pc(i & j) & 1) for j in range(16)] for i in range(16)],
                   dtype=np.int64)
    verdict('Walsh-as-popcount == Sylvester kron (N=16 sanity)', np.array_equal(W16, H16kron))
    w('')

    w(f'=== SUMMARY: {len(PASS)} PASS / {len(FAIL)} FAIL  ({time.time()-t0_all:.1f}s) ===')
    if FAIL:
        w('FAILED CHECKS:')
        for f_ in FAIL:
            w(f'  - {f_}')
    OUT.close()

if __name__ == '__main__':
    main()
