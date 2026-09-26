"""procgen_sumgraph_20260926_run.py -- re-checks every finite claim of
05-knowledge/results/procgen_sumgraph_20260926_reflection_orbits.md and prints [OK] lines.

Session collatz-procgen-20260922, lane "sumgraph" (2026-09-26).

Usage (from the worktree root):
    python3 -u 04-computation/experiments/procgen_sumgraph_20260926_run.py --full \
        > 05-knowledge/results/procgen_sumgraph_20260926.out
Without --full the window W_8 = [4374, 6561] is checked on a deterministic sample (every 16th n plus
every run boundary and its predecessor) instead of completely.
Single process, pure Python, peak RSS well under 100 MB.  Writes nothing except stdout.
"""
import os
import sys
import time
import resource
from itertools import combinations
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_sumgraph_20260926_solver as SV   # noqa: E402
import procgen_sumgraph_20260926_theory as TH   # noqa: E402

FULL = '--full' in sys.argv
N_OK = [0]
T0 = time.time()


class CheckFailed(Exception):
    pass


def check(name, cond, detail=''):
    if not cond:
        raise CheckFailed(f'{name}: {detail}')
    N_OK[0] += 1
    print(f'[OK] {name}: {detail}', flush=True)


def section(title):
    print(f'\n=== {title} (t = {time.time() - T0:.1f} s) ===', flush=True)


# ------------------------------------------------------------------------------------------------
section('A. two and three reflection matchings')

# A.0 Lemma R1
ok = all(TH.msize(n, t) == TH.msize_direct(n, t) for n in range(1, 120) for t in range(1, 2 * n + 2))
check('A.0 Lemma R1', ok, '|M_t| = floor((n - |t-n-1|)/2) for all n < 120, 1 <= t <= 2n+1')

# A.1 Theorem A2
bad = []
cnt = 0
for n in range(3, 151):
    found = set()
    for s, t in combinations(range(3, 2 * n), 2):
        if TH.msize(n, s) + TH.msize(n, t) != n - 1:
            if n <= 40 and TH.has_cycle(n, TH.union_edges(n, [s, t])):
                bad.append(('cycle', n, s, t))
            continue
        E = TH.union_edges(n, [s, t])
        if TH.has_cycle(n, E):
            bad.append(('cycle', n, s, t))
        if TH.is_path_graph(n, E):
            found.add((s, t))
    if found != TH.a2_prediction(n):
        bad.append(('pred', n, sorted(found ^ TH.a2_prediction(n))))
    cnt += len(found)
check('A.1 Theorem A2', not bad, f'two-target unions are acyclic (all pairs, n <= 40) and Hamiltonian exactly as predicted, 3 <= n <= 150 '
      f'({cnt} Hamiltonian pairs)' + (f' BAD {bad[:3]}' if bad else ''))

# A.2 Theorem A3
bad = []
tight = paths = spor = 0
for n in range(3, 81):
    for a, b, c in combinations(range(3, 2 * n), 3):
        if not TH.a3_cond_ii(n, a, b, c):
            continue
        md = TH.maxdeg_union(n, [a, b, c]) if c - a >= n - 2 else 3
        ci = TH.a3_cond_i(n, a, b, c)
        if (md <= 2) != ci:
            bad.append(('i', n, a, b, c))
        expl = (c - a >= n) or (c - a == n - 1 and b == 2 * a - 2) or \
            (n in (6, 7) and (a, b, c) in ((4, 6, n + 2), (n, 2 * n - 4, 2 * n - 2)))
        if ci != expl:
            bad.append(('i-explicit', n, a, b, c))
        if not ci:
            continue
        tight += 1
        if n in (6, 7) and (a, b, c) in ((4, 6, n + 2), (n, 2 * n - 4, 2 * n - 2)):
            spor += 1
        isp = TH.is_path_graph(n, TH.union_edges(n, [a, b, c]))
        paths += isp
        if isp != TH.a3_cond_iii(a, b, c):
            bad.append(('iii', n, a, b, c))
check('A.2 Theorem A3', not bad, f'for 3 <= n <= 80: (i) <=> max degree <= 2 <=> the explicit list; among the {tight} triples with '
      f'(i),(ii) (incl. {spor} sporadic n=6,7 cases), the union is a Hamiltonian path iff (iii) ({paths} paths)'
      + (f' BAD {bad[:3]}' if bad else ''))

# A.3 Corollary A3'
bad = []
for n in range(3, 61):
    for a in range(3, 2 * n):
        for c in range(a + n, 2 * n):
            for b in range(a + 1, c):
                if TH.a3_cond_ii(n, a, b, c) != TH.a3_generic_table(n, a, b, c):
                    bad.append((n, a, b, c))
check("A.3 Corollary A3'", not bad, 'explicit generic table == (ii) for all a<b<c with c-a >= n, n <= 60'
      + (f' BAD {bad[:3]}' if bad else ''))

# A.4 three consecutive squares
hits = []
for j in range(2, 41):
    a, b, c = (j - 1) ** 2, j * j, (j + 1) ** 2
    if a < 3:
        continue
    for n in range(max(3, (c + 2) // 2), 2001):
        if c > 2 * n - 1:
            continue
        if TH.a3_cond_i(n, a, b, c) and TH.a3_cond_ii(n, a, b, c) and TH.a3_cond_iii(a, b, c):
            hits.append((j, n))
            check_path = TH.is_path_graph(n, TH.union_edges(n, [a, b, c]))
            if not check_path:
                hits.append(('NOTPATH', j, n))
check('A.4 Corollary A3\'\' (squares)', hits == [(4, 15), (4, 16), (4, 17)],
      f'three consecutive squares form a Hamiltonian union only for (j,n) in {hits} (j <= 40, n <= 2000)')

# A.4b Pythagorean zigzags
from math import isqrt  # noqa: E402
bad = []
cnt = 0
for t in range(2, 101):
    for s_ in range(1, t):
        u2 = s_ * s_ + t * t
        u = isqrt(u2)
        if u * u != u2 or gcd(s_, t) != 1 or s_ * s_ < 3:
            continue
        a, b, c = s_ * s_, t * t, u2
        ns = [n for n in range(max(3, (c + 2) // 2), t * t + 4) if c <= 2 * n - 1
              and TH.is_path_graph(n, TH.union_edges(n, [a, b, c]))]
        expect = [t * t - 1, t * t] + ([t * t + 1] if t * t == 2 * a - 2 else [])
        E = TH.union_edges(t * t - 1, [a, b, c])
        deg = [0] * (t * t)
        for x, y in E:
            deg[x] += 1
            deg[y] += 1
        ends = sorted(v for v in range(1, t * t) if deg[v] == 1)
        exp_ends = sorted([a, b // 2]) if t % 2 == 0 else sorted([a // 2, a])
        cnt += 1
        if ns != expect or ends != exp_ends or not TH.rotation_orbit_check(t * t - 1, a, b, c):
            bad.append((s_, t, u, ns, ends))
special = [(s_, t) for t in range(2, 10001) for s_ in [isqrt((t * t + 2) // 2)] if 2 * s_ * s_ - 2 == t * t
           and isqrt(s_ * s_ + t * t) ** 2 == s_ * s_ + t * t]
check('A.4b Pythagorean zigzags', not bad and special == [(3, 4)],
      f'for all {cnt} primitive triples s^2+t^2=u^2 (s<t<=100): the three squares give a Hamiltonian union exactly for '
      f'n in {{t^2-1, t^2}} (+ n = 17 for (3,4,5)), with ends {{s^2, t^2/2}} (t even) or {{s^2/2, s^2}} (t odd) at n = t^2-1, and '
      f'rotation by s^2 mod t^2; t^2 = 2s^2-2 with a Pythagorean triple only for (3,4), t <= 10000')

# A.5 zigzag family
bad = []
for k in range(1, 61):
    n, a, b, c = 4 * k - 1, 2 * k + 1, 4 * k, 6 * k + 1
    if n < 3:
        continue
    z = TH.zigzag(k)
    E = sorted(tuple(sorted(e)) for e in TH.union_edges(n, [a, b, c]))
    Z = sorted(tuple(sorted((z[i], z[i + 1]))) for i in range(n - 1))
    if not (E == Z and TH.a3_cond_i(n, a, b, c) and TH.a3_cond_ii(n, a, b, c) and gcd(b - a, c - b) == 1
            and TH.rotation_orbit_check(n, a, b, c)):
        bad.append(k)
check('A.5 zigzag', not bad, 'for k <= 60, M_{2k+1} u M_{4k} u M_{6k+1} on [4k-1] is exactly the path Z_k, (i)-(iii) hold, '
      'and the return map is the rotation by c-b mod c-a')

# A.6 rotation check on all generic tight triples
bad = 0
cnt = 0
for n in range(3, 61):
    for a, b, c in combinations(range(3, 2 * n), 3):
        if c - a >= n and TH.a3_cond_ii(n, a, b, c) and TH.a3_cond_iii(a, b, c):
            cnt += 1
            bad += not TH.rotation_orbit_check(n, a, b, c)
check('A.6 rotation', bad == 0, f'return map to the b-matching is x -> x + (c-b) mod (c-a) on all {cnt} generic tight triples, n <= 60')

# A.7 two-target chains of C_n
hits = []
for n in range(3, 3001):
    T = SV.targets_for(n, 'pow23')
    Ts = set(T)
    if any((s, t) in TH.a2_prediction(n) for s in T for t in T if s < t):
        hits.append(n)
check("A.7 Corollary A2' (C_n)", hits == [3, 7, 8], f'two-target Hamiltonian paths of C_n, 3 <= n <= 3000: n in {hits}')

# ------------------------------------------------------------------------------------------------
section('B. the Collatz-alphabet sum graph C_n')

bad = [n for n in range(3, 1501) if not TH.top_layer_check(n)]
bad2 = []
for n in range(3, 20001):
    P, Q = TH.pow2_le(n), TH.pow3_le(n)
    M = max(P, Q)
    high = [t for t in SV.pow23(2 * n - 1) if M < t <= 2 * n - 1]
    if not set(high) <= {2 * P, 3 * Q}:
        bad2.append(n)
check('B.1 Lemma T', not bad and not bad2, 'targets in (M, 2n-1] are among {2P, 3^(b+1)} for all n <= 20000, and the neighbour '
      'lists of all top vertices are exactly {2P-x, 3^(b+1)-x} for n <= 1500')

WIN = {3: (18, 27), 4: (49, 80), 5: (162, 243), 6: (473, 728), 7: (1319, 2186), 8: (4374, 6561), 9: (11491, 19682)}
deltas = []
for a, (lo, hi) in WIN.items():
    Ps = sorted({TH.pow2_le(m) for m in range(lo, hi)})
    deltas.append((a, [abs(3 ** a - 2 * P) for P in Ps]))
check('B.1b top translations', deltas == [(3, [5]), (4, [17, 47]), (5, [13]), (6, [217, 295]), (7, [139, 1909]), (8, [1631]),
                                          (9, [3299, 13085])],
      'delta = |3^a - 2P| for the values of P met in W_a (below 3^a): ' + str(deltas))

# B.2 independent re-verification to 2200
KNOWN = ('1-3, 5-8, 18-26, 49-63, 65-66, 68-80, 179-194, 224-243, 473-575, 665-728, 1319-1418, 1536-1616, '
         '1620-1663, 1703-1713, 1792-1802, 1920-2114, 2120-2186')
status = {}
nodes_none = []
t1 = time.time()
for n in range(1, 2201):
    r = SV.decide(n, 'pow23', node_limit=400000)
    status[n] = r['status']
    if r['status'] == 'PATH':
        T = SV.targets_for(n, 'pow23')
        if not SV.verify_path(n, r['seq'], T):
            raise CheckFailed(f'bad path at {n}')
    if r['status'] == 'NONE':
        nodes_none.append(r.get('nodes', 0))
    if r['status'] == 'ABORT':
        raise CheckFailed(f'abort at {n}')
ham = [n for n in status if status[n] == 'PATH']
from collections import Counter  # noqa: E402
cnts = Counter(status.values())
check('B.2 Ham set to 2200', TH.fmt_runs(TH.runs_of(ham)) == KNOWN,
      f'independent solver: C_n has a Hamiltonian path for n <= 2200 exactly on {TH.fmt_runs(TH.runs_of(ham))} '
      f'(LOCAL {cnts["LOCAL"]}, PATH {cnts["PATH"]}, NONE {cnts["NONE"]}; {time.time() - t1:.0f} s)')
check('B.2b propagation completeness to 2200', all(x == 0 for x in nodes_none),
      f'all {len(nodes_none)} non-local NONE verdicts are reached by unit propagation alone (0 branching nodes)')

# B.2c paths forced by propagation alone (two leaves)
full = []
for n in range(3, 2201):
    if status[n] != 'PATH':
        continue
    T = SV.targets_for(n, 'pow23')
    E = SV.edges_of(n, T)
    d = SV.degrees(n, E)
    lv = [v for v in range(1, n + 1) if d[v] <= 1]
    if len(lv) != 2:
        continue
    S = SV.State(n, E, lv)
    S.propagate(list(range(1, n + 1)))
    if sum(1 for x in S.st if x == SV.FO) == n - 1:
        full.append(n)
check('B.2c unique forced paths', full == [3, 5, 8, 65, 243], f'two-leaf Hamiltonian n <= 2200 whose path is forced by propagation '
      f'alone (hence unique): {full}')

# B.3 W_8
W8 = (4374, 6561)
if FULL:
    ns8 = list(range(W8[0], W8[1] + 1))
else:
    EXPECT_BOUNDS = []   # filled from the full run in the note; the sample below always includes the window ends
    ns8 = sorted(set(range(W8[0], W8[1] + 1, 16)) | {W8[0], W8[1]})
t1 = time.time()
st8 = {}
maxnodes8 = 0
for n in ns8:
    r = SV.decide(n, 'pow23', node_limit=400000)
    if r['status'] == 'ABORT':
        raise CheckFailed(f'abort at {n}')
    if r['status'] == 'PATH' and not SV.verify_path(n, r['seq'], SV.targets_for(n, 'pow23')):
        raise CheckFailed(f'bad path at {n}')
    st8[n] = r['status']
    if r['status'] == 'NONE':
        maxnodes8 = max(maxnodes8, r.get('nodes', 0))
ham8 = [n for n in ns8 if st8[n] == 'PATH']
none8 = [n for n in ns8 if st8[n] == 'NONE']
local8 = [n for n in ns8 if st8[n] == 'LOCAL']
if FULL:
    status.update(st8)
    check('B.3 W_8', not local8, f'W_8 = [4374,6561] (all admissible): Hamiltonian exactly on {TH.fmt_runs(TH.runs_of(ham8))}; '
          f'{len(none8)} NONE, all by propagation alone: {maxnodes8 == 0}; {time.time() - t1:.0f} s')
    check('B.3b propagation completeness on W_8', maxnodes8 == 0, 'no branching node was needed for any NONE in W_8')
else:
    check('B.3 W_8 (sample)', not local8, f'{len(ns8)} sampled n: PATH {len(ham8)}, NONE {len(none8)}; {time.time() - t1:.0f} s')


# B.3c the single (resolvable) choke of W_8 at the defect 512 = P/8
K512 = [512, 1536, 1675, 2421, 2560, 2977, 3584, 4608, 4886, 5025]
w8c = {}
for n in (5024, 5025, 5100, 5214, 5215):
    T = SV.targets_for(n, 'pow23')
    E = SV.edges_of(n, T)
    S = SV.State(n, E, [4096])
    try:
        S.propagate(list(range(1, n + 1)))
        w8c[n] = None
    except SV.Contra:
        w8c[n] = sorted(S.core())
check('B.3c W_8 choke at 512', w8c[5024] is None and w8c[5215] is None and all(w8c[n] == K512 for n in (5025, 5100, 5214)),
      f'base propagation (end 4096 forced) finds the choke core {K512} at n = 5025, 5100, 5214 and no conflict at 5024, 5215')

# B.4 hand-proved conflicts: hypotheses of each instance, for every n of the claimed range
def nb_sets(n):
    nb = TH.nbrs(n)
    return [set(a) for a in nb]


def leaves_of(n, nb):
    return sorted(v for v in range(1, n + 1) if len(nb[v]) <= 1)


def sat(nb, v, ws):
    """v is adjacent to every w in ws and every such w has degree exactly 2."""
    return all(w in nb[v] and len(nb[w]) == 2 for w in ws)


inst_fail = []


def instance(name, rng, leaves, conds):
    for n in rng:
        nb = nb_sets(n)
        if leaves is not None and leaves_of(n, nb) != leaves:
            inst_fail.append((name, n, 'leaves', leaves_of(n, nb)))
            return
        for kind, v, data in conds:
            ok_ = (nb[v] == set(data)) if kind == 'N' else sat(nb, v, data)
            if not ok_:
                inst_fail.append((name, n, kind, v, sorted(nb[v])))
                return
        if status.get(n) not in (None, 'NONE'):
            inst_fail.append((name, n, 'status', status.get(n)))
            return
    check(f'B.4 {name}', True, f'hypotheses hold for every n in [{rng.start}, {rng.stop - 1}] and the solver agrees (NONE)')


instance('W5 K32 (ends forced)', range(162, 179), [64, 128],
         [('N', 32, [49, 96]), ('S', 96, [147, 160])])
instance('W5 K32+K16 (second order)', range(195, 207), [128],
         [('N', 32, [49, 96]), ('S', 96, [147, 160]), ('S', 95, [148, 161]), ('N', 33, [31, 48, 95]),
          ('N', 48, [16, 33, 80, 195]), ('S', 48, [195]), ('S', 65, [178, 191]), ('N', 16, [11, 48, 65, 112]),
          ('S', 112, [131, 144])])
instance('W5 K32 doubly choked', range(207, 211), [128],
         [('N', 32, [49, 96]), ('S', 96, [147, 160]), ('S', 49, [194, 207])])
instance('W5 K32 + K16 (first order)', range(211, 224), [128],
         [('N', 32, [49, 96, 211]), ('S', 96, [147, 160]), ('S', 49, [194, 207]), ('S', 32, [211]),
          ('N', 16, [11, 48, 65, 112]), ('S', 48, [195, 208]), ('S', 65, [178, 191]), ('S', 112, [131, 144])])
instance('W6 K64 (ends forced)', range(576, 665), [256, 512],
         [('N', 64, [17, 179, 192, 448]), ('S', 179, [333, 550]), ('S', 192, [320, 537]), ('S', 448, [281, 576])])
instance('W7 K256 (ends forced)', range(1419, 1536), [512, 1024],
         [('N', 256, [473, 768]), ('S', 768, [1280, 1419])])
instance('W7 K128(second order)+K256', range(1664, 1703), [1024],
         [('N', 256, [473, 768]), ('S', 768, [1280, 1419]), ('N', 128, [115, 384, 601, 896]),
          ('S', 601, [1447, 1586]), ('S', 896, [1152, 1291]), ('N', 384, [128, 345, 640, 1664]), ('S', 384, [1664]),
          ('N', 345, [167, 384, 679]), ('S', 679, [1369, 1508])])
instance('W7 K256 doubly choked', range(1714, 1792), [1024],
         [('N', 256, [473, 768]), ('S', 768, [1280, 1419]), ('S', 473, [1575, 1714])])
instance('W7 K128+K256 (first order)', range(1803, 1920), [1024],
         [('N', 256, [473, 768, 1792]), ('S', 256, [1792]), ('S', 768, [1280, 1419]), ('S', 473, [1575, 1714]),
          ('N', 128, [115, 384, 601, 896]), ('S', 384, [1664, 1803]), ('S', 601, [1447, 1586]), ('S', 896, [1152, 1291])])
check('B.4 all instances', not inst_fail, 'no hypothesis failed' + (f' FAIL {inst_fail[:3]}' if inst_fail else ''))

# B.4b parametric families W1-W3 (endpoints are linear forms in P = 2^p and B = 3^(a-1))
from math import log2  # noqa: E402
crit = {'W1': lambda P, B: 3 * P > 4 * B and 7 * P < 12 * B, 'W2': lambda P, B: 3 * P > 4 * B and 5 * P < 8 * B,
        'W3': lambda P, B: 3 * P > 4 * B and 2 * P < 3 * B}
fam_bad = []
fam_cnt = Counter()
listing = []
for a in range(3, 201):
    B, P = TH.level(a)
    for fam, (lo, hi) in TH.w_ranges(a).items():
        ne = lo < hi
        if ne != crit[fam](P, B):
            fam_bad.append(('criterion', a, fam))
        if not ne:
            continue
        fam_cnt[fam] += 1
        # hypotheses are monotone in n, so the two ends of the range suffice; sharpness: they fail just outside
        if not (TH.w_hypotheses(a, fam, lo) and TH.w_hypotheses(a, fam, hi - 1)):
            fam_bad.append(('hyp', a, fam))
        if TH.w_hypotheses(a, fam, lo - 1) or TH.w_hypotheses(a, fam, hi):
            fam_bad.append(('not sharp', a, fam))
        if a <= 10:
            if not all(TH.w_hypotheses(a, fam, n) for n in range(lo, hi)):
                fam_bad.append(('hyp-all', a, fam))
        if a <= 12:
            listing.append(f'{fam}@a={a}: [{lo},{hi - 1}]')
        if a <= 8:
            for n in range(lo, hi):
                if status.get(n) not in (None, 'NONE', 'LOCAL'):
                    fam_bad.append(('solver', a, fam, n))
check('B.4b families W1-W3', not fam_bad,
      'for 3 <= a <= 200: W1/W2/W3 ranges are nonempty exactly when 4/3 < rho < 12/7, 8/5, 3/2 (rho = P/3^(a-1)); their hypotheses '
      'hold on the whole range (both ends, monotone; every n for a <= 10) and fail just outside; the solver agrees for a <= 8. '
      'Ranges for a <= 12: ' + '; '.join(listing) + (f' BAD {fam_bad[:4]}' if fam_bad else ''))
# W3 disjointness for every level: a coincidence alpha1*P + beta1*B = alpha2*P + beta2*B forces rho = P/B to a fixed rational
from fractions import Fraction as Fr  # noqa: E402


def _is_23(r):
    num, den = r.numerator, r.denominator
    for pr in (2, 3):
        while num % pr == 0:
            num //= pr
        while den % pr == 0:
            den //= pr
    return num == 1 and den == 1


R4 = [(Fr(1, 4), 0), (Fr(5, 4), 0), (Fr(-3, 4), 3), (Fr(9, 4), -1), (Fr(1, 4), 2)]
R8 = [(Fr(1, 8), 0), (Fr(13, 8), 0), (Fr(-3, 8), 3), (Fr(9, 8), 0), (Fr(-7, 8), 3), (Fr(17, 8), -1), (Fr(1, 8), 2)]
coinc = []
for a1, b1 in R4:
    for a2, b2 in R8:
        if a1 == a2:
            if b1 == b2:
                coinc.append(((a1, b1), 'identical'))
            continue
        r = Fr(b2 - b1) / (a1 - a2)
        if Fr(4, 3) < r < Fr(3, 2) and _is_23(r):
            coinc.append(((a1, b1), (a2, b2), r))
check('B.4b2 W3 disjointness (all levels)', not coinc, 'no element of R(P/4) can equal an element of R(P/8) for rho = 2^p/3^(a-1) '
      'in (4/3, 3/2): each coincidence forces rho to a rational that is not of the form 2^i/3^j or lies outside the range')
check('B.4c family densities', True, f'levels 3..200 with a nonempty range: W1 {fam_cnt["W1"]}, W2 {fam_cnt["W2"]}, W3 {fam_cnt["W3"]} '
      f'of 198 (equidistribution predicts log2(9/7) = {log2(9 / 7):.3f}, log2(6/5) = {log2(6 / 5):.3f}, log2(9/8) = {log2(9 / 8):.3f} '
      f'of the levels: observed {fam_cnt["W1"] / 198:.3f}, {fam_cnt["W2"] / 198:.3f}, {fam_cnt["W3"] / 198:.3f})')
for n in (40960, 42664):
    T = SV.targets_for(n, 'pow23')
    E = SV.edges_of(n, T)
    d = SV.degrees(n, E)
    lv = [v for v in range(1, n + 1) if d[v] <= 1]
    S = SV.State(n, E, lv)
    try:
        S.propagate(list(range(1, n + 1)))
        core = None
    except SV.Contra:
        core = sorted(S.core())
    check(f'B.4d W_10 at n = {n}', lv == [16384, 32768] and core == [8192, 24576, 34473, 40960],
          f'solver: ends forced {lv}; propagation contradiction with core {core} = {{P/4, 3P/4, 3B-3P/4, 5P/4}} (Proposition W1)')

# B.5 first-order Choke Lemma: soundness and coverage
sound_bad = []
cover = Counter()
for n in sorted(status):
    if n < 5 or status[n] == 'LOCAL':
        continue
    v, conf = TH.first_order_verdict(n)
    if v == 'NONE' and status[n] != 'NONE':
        sound_bad.append(n)
    if status[n] == 'NONE':
        cover['first-order' if v == 'NONE' else 'deeper'] += 1
check('B.5 Choke Lemma soundness', not sound_bad,
      f'first-order conflicts never refute a Hamiltonian n (all admissible n checked above); NONE verdicts explained at first '
      f'order: {cover["first-order"]}, needing deeper propagation: {cover["deeper"]}')

# B.6 boundary catalogue
cat = []
for n in sorted(status):
    if n - 1 in status and n > 150 and status[n] != 'LOCAL' and status[n - 1] != 'LOCAL' and \
            (status[n] == 'PATH') != (status[n - 1] == 'PATH'):
        w = TH.boundary_word(n)
        cat.append((n, 'start' if status[n] == 'PATH' else 'stop', w))
for n, kind, w in cat:
    L, t, path = w
    print(f'    boundary {n} ({kind}): n = r_{t}(v), v = {path[-1]} reached from {path[0]} in {L} lower reflections: {path}')
check('B.6 boundary words', all(w is not None and w[0] <= 2 for _, _, w in cat),
      f'every one of the {len(cat)} inner boundaries (n > 150) is T - w(s): T a top target, s a power of 2 or 3, |w| <= 2; '
      f'|w| = 0 for {sum(1 for _, _, w in cat if w[0] == 0)}, |w| <= 1 for {sum(1 for _, _, w in cat if w[0] <= 1)}')
for a in (5, 6, 7, 8):
    lo, hi = WIN[a]
    step = 8 if a == 8 else 1
    dist = Counter()
    for n in range(lo, hi + 1, step):
        w = TH.boundary_word(n)
        dist[w[0] if w else None] += 1
    tot = sum(dist.values())
    check(f'B.6b base rate W_{a}', True, f'{"every 8th" if step > 1 else "all"} n in W_{a} ({tot}): |w| = 0: {dist[0] / tot:.3f}; '
          f'|w| <= 1: {(dist[0] + dist[1]) / tot:.3f}; |w| <= 2: {(dist[0] + dist[1] + dist[2]) / tot:.3f}')

# ------------------------------------------------------------------------------------------------
section('C. square sums near the thresholds')
sq = {}
for n in range(1, 41):
    r = SV.decide(n, 'squares', node_limit=400000)
    sq[n] = r['status']
    if r['status'] == 'PATH' and not SV.verify_path(n, r['seq'], SV.targets_for(n, 'squares')):
        raise CheckFailed(f'bad square path {n}')
hamq = [n for n in sq if sq[n] == 'PATH']
check('C.1 square-sum paths', TH.fmt_runs(TH.runs_of(hamq)) == '1, 15-17, 23, 25-40',
      f'independent solver: Q_n has a Hamiltonian path for n <= 40 exactly on {TH.fmt_runs(TH.runs_of(hamq))} (A090461, CITED)')
check('C.1b statuses 18-24', [sq[n] for n in range(18, 25)] == ['LOCAL', 'NONE', 'NONE', 'NONE', 'NONE', 'PATH', 'NONE'],
      '18 LOCAL (three leaves 16,17,18); 19-22 and 24 NONE by propagation; 23 PATH')

nb19 = [set(a) for a in TH.nbrs(19, 'squares')]
c19 = (leaves_of(19, nb19) == [16, 18] and nb19[3] == {1, 6, 13} and sat(nb19, 1, [8, 15]) and sat(nb19, 6, [10, 19]))
check('C.2 Q_19', c19, 'ends forced {16,18}; N(3) = {1,6,13}, 1 saturated by the degree-2 vertices 8,15 and 6 by 10,19: '
      '3 would have to be an end')
nb20 = [set(a) for a in TH.nbrs(20, 'squares')]
c20 = (leaves_of(20, nb20) == [18] and sat(nb20, 5, [4, 11, 20]) and nb20[3] == {1, 6, 13} and sat(nb20, 1, [8, 15])
       and sat(nb20, 6, [10, 19]))
check('C.2b Q_20', c20, 'over-saturation at 5 (degree-2 neighbours 4,11,20: free end in {4,11,20}) and the choke at 3 '
      '(free end in {3,8,10,15,19}) are disjoint')
for n in (21, 22):
    T = SV.targets_for(n, 'squares')
    E = SV.edges_of(n, T)
    S = SV.State(n, E, [18])
    try:
        S.propagate(list(range(1, n + 1)))
        K = None
    except SV.Contra:
        K = sorted(S.core())
    refuted = []
    for e in (K or []):
        if e == 18:
            continue
        S2 = SV.State(n, E, [18, e])
        try:
            S2.propagate(list(range(1, n + 1)))
        except SV.Contra:
            refuted.append((e, len(S2.core())))
    nbn = [set(a) for a in TH.nbrs(n, 'squares')]
    check(f'C.3 Q_{n}', K == [2, 7, 9, 18] and len(refuted) == 3 and nbn[2] == {7, 14} and nbn[9] == {7, 16}
          and nbn[18] == {7}, f'lock: 18 and 9 (degree 1, 2) saturate 7, so 2 (N = {{7,14}}) is choked: base core {K}; '
          f'the three candidates are refuted by propagation (core sizes {refuted})')
nb23 = [set(a) for a in TH.nbrs(23, 'squares')]
check('C.3b Q_23', 23 in nb23[2] and sq[23] == 'PATH', 'the lock opens at 23 = r_25(2): vertex 2 gains the partner 23')
T = SV.targets_for(24, 'squares')
E = SV.edges_of(24, T)
S = SV.State(24, E, [18])
try:
    S.propagate(list(range(1, 25)))
    K = None
except SV.Contra:
    K = sorted(S.core())
sizes = []
for e in K:
    if e == 18:
        continue
    S2 = SV.State(24, E, [18, e])
    try:
        S2.propagate(list(range(1, 25)))
        sizes.append(None)
    except SV.Contra:
        sizes.append(len(S2.core()))
check('C.4 Q_24', K is not None and None not in sizes, f'base core has {len(K)} vertices; all {len(sizes)} candidates refuted, '
      f'cores of sizes {sorted(sizes)}')

# ------------------------------------------------------------------------------------------------
section('D. rotation bridge (numerical companions)')
rows = [TH.rotation_orbit_check(n, 9, 16, 25) for n in (15, 16, 17)]
check('D.1 square window as rotations', all(rows), 'Q_15, Q_16, Q_17 chains: every return to the 16-matching is x -> x + 9 mod 16 '
      '(gcd(9,16) = 1, a single rotation orbit)')

ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
peak = ru / 1e6 if sys.platform == 'darwin' else ru / 1e3
print(f'\nRUN {N_OK[0]} checks passed; wall {time.time() - T0:.0f} s; peak RSS {peak:.0f} MB', flush=True)
