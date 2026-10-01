#!/usr/bin/env python3
"""procgen_edim_20261001_run.py -- runner for the edge multiset dimension lane (procgen_edim, 2026-10-01).

Re-checks every claim of 05-knowledge/results/procgen_edim_20261001_edge_multiset_dimension.md:
  S1  the cited certificates of Allikvere (arXiv:2608.09983, Table 1) for Q6, Q7, Q8      [VERIFIED]
  S2  elementary lemmas (trivial stabilizer, antipodal reversal, alternating sum, counting
      bound edim_m(Q6) >= 7, entropy lower bound table)                                    [PROVED + numeric checks]
  S3  Aut(Q5)-orbit representatives of a-subsets of Q5, a=0..15 (C), re-validated in numpy:
      distinct, minimal in orbit, count = Burnside                                          [FINITE-EXACT]
  S4  METHOD A exhaustive search, k=1..15 (C): leaf counts re-derived by a numpy DP;
      no resolving set for k<=14; resolving 15-sets found                                   [FINITE-EXACT]
  S5  METHOD B (independent normal form + implementation), same checks                      [FINITE-EXACT]
  S6  k=15: A and B give the same 229 Aut(Q6)-orbits; every rep verified in pure Python,
      stabilizer 1; the paper's 15-set is among them  => edim_m(Q6) = 15                     [FINITE-EXACT]
  S7  k=14 near-resolving cross-check: orbits with defect <= 2 agree between A and B        [FINITE-EXACT]
  S8  normal-form test on random 14- and 15-sets (the WLOG reductions of A and B)
  S9  tournament-structured landmark sets of the n=5 tiling cube (none resolving)           [FINITE-EXACT]
  S10 new upper bounds for Q7..Q12 by explicit sets (exact verification)                     [VERIFIED]
  S11 union bounds: q=1/2 reproduces the paper (U(10)>1>U(11)); sparse q gives
      edim_m(Q_d) <= M_d for d=11..18,20,24,28,32 (double precision, 1e-6 margin per factor)  [VERIFIED]
Prints to stdout only; ends with ALL CHECKS PASSED.
Usage: python3 -u procgen_edim_20261001_run.py [--quick]   (--quick = k<=12 smoke test, not a certification)
Optional env PROCGEN_EDIM_BUILD=<dir> keeps the compiled helpers / rep files there (default: temp dir, removed).
Runtime ~1 h single core (S5 k=15 dominates); peak RSS < 150 MB.
"""
import os, sys, subprocess, tempfile, random, time, re
from math import log
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_edim_20261001_lib as L

NCHECK = 0
def ok(cond, msg):
    global NCHECK
    if not cond:
        print('[FAIL]', msg, flush=True); sys.exit(1)
    NCHECK += 1; print('[OK]', msg, flush=True)

BUILD = os.environ.get('PROCGEN_EDIM_BUILD') or tempfile.mkdtemp(prefix='procgen_edim_')
os.makedirs(BUILD, exist_ok=True)
CLEANUP = 'PROCGEN_EDIM_BUILD' not in os.environ
def build(name, extra=()):
    src = os.path.join(HERE, 'procgen_edim_20261001_%s.c' % name); exe = os.path.join(BUILD, name)
    subprocess.check_call(['cc', '-O3', '-o', exe, src, '-lm'] + list(extra)); return exe
def run(args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout
def section(t): print('\n==== %s ====' % t, flush=True)
T0 = time.time()
QUICK = '--quick' in sys.argv          # smoke test only: k <= 12, skips S6/S7 (NOT a certification run)
KMAX = 12 if QUICK else 15
if QUICK: print('QUICK MODE: smoke test only (k<=12); the certification run is without --quick', flush=True)

# ======================================================================= S1
section('S1 cited certificates (Allikvere 2026, Table 1)')
PAPER = {6: '02283022a042a00a', 7: '3f904303076d2c0bc1ab7b64b3ae2ff0',
         8: '6e26bea413a06ca54c8e0782b28229de' + 'b50e4bb04481b789efe84aea02c1c88b'}
PAPER_SIZE = {6: 15, 7: 63, 8: 115}
for d, h in PAPER.items():
    S = L.mask_to_set(int(h, 16))
    ok(len(S) == PAPER_SIZE[d] and L.is_resolving(d, S) and L.is_resolving_np(d, S),
       'paper certificate Q%d: |S|=%d is edge-multiset resolving (pure-python + numpy checkers)' % (d, len(S)))

# ======================================================================= S2
section('S2 elementary lemmas')
rng = random.Random(20261001)
def hist(d, e, S):
    u, v = e; H = [0] * d
    for s in S: H[min(L.hd(u, s), L.hd(v, s))] += 1
    return H
for d in (5, 6, 7):
    N = 1 << d; E = L.edges(d); good_rev = good_alt = True
    for trial in range(20):
        S = rng.sample(range(N), rng.randint(1, N - 1))
        for (u, v) in E:
            i = (u ^ v).bit_length() - 1
            H = hist(d, (u, v), S); Hb = hist(d, (u ^ (N - 1), v ^ (N - 1)), S)
            if Hb != H[::-1]: good_rev = False
            U = (N - 1) ^ (1 << i)
            walsh = sum((-1) ** bin(U & s).count('1') for s in S)
            if sum((-1) ** r * H[r] for r in range(d)) != (-1) ** bin(U & u).count('1') * walsh: good_alt = False
    ok(good_rev, 'Q%d antipodal reversal H_{ebar}(r) = H_e(d-1-r) (20 random S, all edges)' % d)
    ok(good_alt, 'Q%d alternating sum: sum_r (-1)^r H_e(r) = chi_U(u) * S^(U), U = 1-vector xor e_i' % d)
# odd |S| in even d: colliding pairs come in pairs (no palindromic histogram of odd total)
par_ok = all(L.colliding_pairs(6, rng.sample(range(64), k)) % 2 == 0 for k in (7, 9, 11, 13, 15) for _ in range(40))
ok(par_ok, 'Q6, |S| odd: number of colliding edge pairs is even (200 random sets)')
# counting bound: >= 192-12m edges avoid S and its antipodes, their histograms live on levels 1..4
ok(all(192 - 12 * m > L.comb(m + 3, 3) for m in range(1, 7)) and 192 - 12 * 7 <= L.comb(10, 3),
   'counting bound: 192-12m > C(m+3,3) for m<=6, fails at m=7  => edim_m(Q6) >= 7')
cnt_ok = True
for _ in range(50):
    m = rng.randint(1, 20); S = rng.sample(range(64), m)
    z = sum(1 for e in L.edges(6) if hist(6, e, S)[0] == 0 and hist(6, e, S)[5] == 0)
    cnt_ok &= (z >= 192 - 12 * m)
ok(cnt_ok, 'sanity: #edges with H(0)=H(5)=0 is >= 192-12|S| (50 random sets)')
ENT = {6: 4, 7: 5, 8: 5, 10: 6, 12: 8, 16: 11, 20: 15, 32: 28, 64: 81, 128: 275, 256: 1159, 1024: 50116}
ok(all(L.entropy_mmin(d) == m for d, m in ENT.items()),
   'entropy lower bound edim_m(Q_d) >= m_ent(d): ' + ', '.join('%d:%d' % kv for kv in ENT.items()))

# ======================================================================= S3
section('S3 Aut(Q5)-orbit representatives of a-subsets of Q5, a=0..%d' % KMAX)
q5 = build('q5orbits'); A_SEARCH = build('searchA'); B_SEARCH = build('searchB'); CANON6 = build('canon6')
out = run([q5, str(KMAX), os.path.join(BUILD, 'r')])
BURN5 = L.burnside_subset_orbits(5)
ok(BURN5[:16] == [1, 1, 5, 10, 47, 131, 472, 1326, 3779, 9013, 19963, 38073, 65664, 98804, 133576, 158658],
   'Burnside counts of a-subset orbits of Q5, a=0..15: 1,1,5,10,47,131,472,1326,3779,9013,19963,38073,65664,98804,133576,158658')
REP = lambda a: os.path.join(BUILD, 'r_%d.txt' % a)
for a in range(KMAX + 1):
    n = L.check_q5_reps(a, REP(a), BURN5)
    ok(True, 'reps a=%d: %d sets, pairwise distinct, each minimal in its Aut(Q5)-orbit, count = Burnside' % (a, n))

# ======================================================================= S4/S5
def exhaustive(method, exe, kmax=15, maxdef=0):
    res = {}
    for k in range(1, kmax + 1):
        for a in range((k + 1) // 2, k + 1):
            o = run([exe, str(k), str(a), REP(a)] + ([str(maxdef)] if maxdef else []))
            summ = [l for l in o.splitlines() if l.startswith('SUMMARY')][0]
            f = dict(x.split('=') for x in summ.split()[1:])
            res[(k, a)] = (int(f['leaves']), int(f['found']), int(f['near']),
                           [l for l in o.splitlines() if l.startswith('RESOLVING')],
                           [l for l in o.splitlines() if l.startswith('NEAR')])
    return res
for method, exe in (('A', A_SEARCH), ('B', B_SEARCH)):
    section('S%d METHOD %s exhaustive search, k=1..%d' % (4 if method == 'A' else 5, method, KMAX))
    t = time.time(); R = exhaustive(method, exe, KMAX)
    globals()['RES_' + method] = R
    for k in range(1, KMAX + 1):
        lv = sum(R[(k, a)][0] for a in range((k + 1) // 2, k + 1))
        ok(all(R[(k, a)][0] == L.leafcount(k, a, method, REP(a)) for a in range((k + 1) // 2, k + 1)),
           'method %s k=%d: leaf counts per a equal the independent numpy DP count (total %d)' % (method, k, lv))
        fnd = sum(R[(k, a)][1] for a in range((k + 1) // 2, k + 1))
        if k <= 14: ok(fnd == 0, 'method %s k=%d: no edge-multiset resolving set' % (method, k))
        else: ok(fnd > 0, 'method %s k=15: %d resolving leaves found' % (method, fnd))
    print('  (method %s time %.0f s)' % (method, time.time() - t), flush=True)

# ======================================================================= canonical forms (C helper)
def canon(lines):
    o = subprocess.run([CANON6], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout
    return [tuple(l.split()) for l in o.splitlines()]

# ======================================================================= S6
def s6():
    global OA
    section('S6 k=15: orbit comparison and verification')
    orbA = canon([l for a in range(8, 16) for l in RES_A[(15, a)][3]])
    orbB = canon([l for a in range(8, 16) for l in RES_B[(15, a)][3]])
    OA = set(x[0] for x in orbA); OB = set(x[0] for x in orbB)
    ok(OA == OB, 'methods A and B find the same set of Aut(Q6)-orbits of resolving 15-sets')
    ok(len(OA) == 229, 'number of Aut(Q6)-orbits of edge-multiset resolving 15-sets of Q6 = 229')
    ok(all(x[1] == '1' for x in orbA + orbB), 'every resolving 15-set found has trivial stabilizer in Aut(Q6) (orbit size 46080)')
    ok(all(L.is_resolving(6, L.mask_to_set(int(h, 16))) for h in OA), 'all 229 canonical representatives verified resolving (pure Python)')
    pc = canon([PAPER[6]])[0][0]
    ok(pc in OA, 'the paper\'s 15-set (canonical %s) is one of the 229 orbits' % pc)
    ok(L.stabilizer_size(6, L.mask_to_set(int(PAPER[6], 16))) == 1, 'pure-Python stabilizer of the paper set in Aut(Q6) is trivial')
    print('RESULT: edim_m(Q6) = 15  (no resolving set of size <= 14 by methods A and B; 229 orbits of size 15)')
    ORBFILE = os.path.join(HERE, '..', '..', '05-knowledge', 'results', 'procgen_edim_20261001_q6_resolving15_orbits.txt')
    if os.path.exists(ORBFILE):
        filed = set(l.strip() for l in open(ORBFILE) if l.strip() and not l.startswith('#'))
        ok(filed == OA, 'recomputed orbit list equals the deposited list procgen_edim_20261001_q6_resolving15_orbits.txt')
    print('  first five canonical 15-set masks:', ' '.join(sorted(OA)[:5]))

# ======================================================================= S7
def s7():
    section('S7 k=14 near-resolving cross-check (defect = 192 - #distinct histograms <= 2)')
    def near14(exe):
        lines = []
        for a in range(7, 15):
            o = run([exe, '14', str(a), REP(a), '2'])
            lines += [l for l in o.splitlines() if l.startswith('NEAR') or l.startswith('RESOLVING')]
        ok(not any(l.startswith('RESOLVING') for l in lines), 'no resolving 14-set (re-run with defect output)')
        by = {}
        for l, c in zip(lines, canon(lines)):
            dfc = int(re.search(r'def=(\d+)', l).group(1)); by.setdefault(dfc, set()).add(c[0])
        return by
    NA = near14(A_SEARCH); NB = near14(B_SEARCH)
    ok(NA == NB, 'defect<=2 14-set orbits agree between A and B: ' + ', '.join('defect %d: %d orbits' % (k, len(v)) for k, v in sorted(NA.items())))
    ok(sorted((k, len(v)) for k, v in NA.items()) == [(1, 1), (2, 16)] and NA[1] == {'000001810690226d'},
       'k=14: minimum defect 1, attained by exactly one orbit (canonical 000001810690226d); 16 orbits of defect 2')
    mind = min(NA) if NA else None
    smp = sorted(NA[mind])[:20] if NA else []
    ok(all(L.defect(6, L.mask_to_set(int(h, 16))) == mind for h in smp), 'minimum defect at k=14 is %s (20 orbit reps rechecked in Python)' % mind)


if not QUICK:
    s6(); s7()

# ======================================================================= S8
section('S8 normal-form (WLOG) test on random sets')
G5 = L.aut_perms(5)
def canon5(mask):
    best = None
    for g in G5:
        im = sum(1 << g[v] for v in range(32) if (mask >> v) & 1)
        if best is None or im < best[0]: best = (im, g)
    return best
repsets = {a: set(int(l, 16) for l in open(REP(a))) for a in range(KMAX + 1)}
nf_ok = True
for k in (KMAX - 1, KMAX):
    for _ in range(25):
        S = rng.sample(range(64), k)
        for rule in ('max', 'min'):
            rep, Bv, beta = L.normal_form(S, rule, canon5)
            a = bin(rep).count('1'); b = len(Bv); dl = a - b
            T = L.mask_to_set(rep) + [v + 32 for v in Bv]
            bt = [sum(1 if not (s >> i) & 1 else -1 for s in T) for i in range(5)]
            good = (a + b == k and dl >= 0 and rep in repsets[a] and
                    (all(abs(x) <= dl for x in bt) if rule == 'max' else all(abs(x) >= dl for x in bt)))
            c1 = canon(['%016x' % L.set_to_mask(S)])[0][0]; c2 = canon(['%016x' % L.set_to_mask(T)])[0][0]
            nf_ok &= good and c1 == c2
ok(nf_ok, 'normal forms of 100 random %d/%d-sets' % (KMAX - 1, KMAX) + ' (both rules) lie in the searched domains and in the same Aut(Q6)-orbit')

# ======================================================================= S9
section('S9 tournament-structured landmark sets in the n=5 tiling cube Q6 (THM-474)')
tiles, TT = L.tiling_tournaments(5)
iso = {}
for t in range(64): iso.setdefault(L.tour_canon(TT[t]), []).append(t)
cls = list(iso.values()); Hs = {}
for t in range(64): Hs.setdefault(L.ham_paths(TT[t]), []).append(t)
ok(sorted(len(c) for c in cls) == [1, 1, 1, 1, 3, 5, 5, 5, 9, 9, 11, 13], '64 tilings fall into 12 tournament isomorphism classes, sizes 1,1,1,1,3,5,5,5,9,9,11,13')
def best_union(parts):
    best = None; nres = 0
    for r in range(1, len(parts) + 1):
        for comb_ in L.itertools.combinations(range(len(parts)), r):
            S = [v for i in comb_ for v in parts[i]]
            if len(S) in (0, 64): continue
            c = L.colliding_pairs(6, S); nres += (c == 0)
            if best is None or c < best[0]: best = (c, len(S))
    return nres, best
nres, best = best_union(cls)
ok(nres == 0 and best == (4, 30), 'no union of isomorphism classes is resolving (4095 unions; best: 4 colliding pairs at size 30)')
nres, best = best_union(list(Hs.values()))
ok(nres == 0, 'no union of H(T)-level sets is resolving (127 unions; best %d colliding pairs at size %d)' % best)
gt = [tiles.index((6 - b, 6 - a)) for (a, b) in tiles]   # grid transpose (x,y)->(n+1-y,n+1-x) on tiles
ok(sorted(gt) == list(range(6)) and gt != list(range(6)), 'grid transpose G is a nontrivial coordinate permutation of Q6 (an automorphism)')
Ht = lambda t: L.ham_paths(TT[t])
img = lambda t: sum(((t >> j) & 1) << gt[j] for j in range(6))
ok(all(Ht(img(t)) == Ht(t) for t in range(64)), 'H(T) is G-invariant on tilings => every union of H-levels has nontrivial stabilizer (cannot resolve)')

# ======================================================================= S10
section('S10 explicit small resolving sets for Q7..Q12 (found by simulated annealing; exact verification)')
CERT = {   # d -> explicit landmark set (vertex v = integer with binary coordinates), found by procgen_edim_20261001_sad.c
    7: [5, 57, 54, 104, 35, 109, 115, 49, 6, 39, 55, 102, 15, 32, 85, 97, 21, 113, 41],
    8: [227, 92, 87, 110, 36, 32, 218, 61, 41, 64, 31, 104, 81, 251, 197, 101, 97, 112, 59, 67,
         231, 242, 80, 212, 47, 52],
    9: [504, 67, 283, 159, 247, 106, 267, 339, 3, 411, 363, 127, 373, 321, 244, 483, 410, 208, 160, 70,
         193, 34, 481, 480, 385, 100, 303, 24, 484, 380, 507, 226, 64, 117, 447, 505, 375, 131],
    10: [923, 564, 857, 521, 984, 14, 342, 541, 736, 364, 117, 642, 194, 73, 357, 207, 753, 683, 527, 338,
         578, 901, 695, 365, 977, 469, 327, 864, 929, 804, 876, 593, 70, 148, 273, 470, 882, 515, 893, 74,
         321, 102, 773, 550, 68, 461, 589, 869],
    11: [373, 1353, 1043, 1642, 1460, 1894, 1859, 388, 642, 1890, 1509, 1251, 750, 677, 1433, 818, 180, 1591, 1892, 170,
         618, 65, 1151, 522, 110, 1884, 808, 871, 18, 629, 32, 421, 320, 1790, 2016, 609, 186, 78, 1392, 371,
         1614, 892, 950, 1048, 120, 565, 1595, 787, 1345, 1303, 794, 470, 1775, 318, 646, 1374, 494, 593, 1264, 1835,
         44, 769, 270, 550, 85],
    12: [1587, 3991, 3997, 1428, 1698, 873, 631, 3526, 3853, 798, 2068, 1901, 2598, 1766, 832, 421, 2726, 3622, 202, 335,
         2600, 388, 552, 311, 262, 3077, 1509, 3330, 2343, 164, 1091, 424, 173, 3631, 110, 599, 3286, 306, 1138, 101,
         2642, 1006, 64, 481, 83, 2586, 2322, 989, 2366, 2177, 3710, 314, 1071, 275, 2049, 3273, 891, 1439, 1014, 3918,
         2022, 796, 305, 781, 3211, 26, 2071, 3238, 3263, 2926, 709, 2627, 512, 763, 1119, 3057],
}
for d, S in sorted(CERT.items()):
    good = L.is_resolving_np(d, S) and (d > 8 or L.is_resolving(d, S))
    ok(good and len(set(S)) == len(S), 'Q%d: explicit set of size %d is edge-multiset resolving  => edim_m(Q%d) <= %d' % (d, len(S), d, len(S)))

# ======================================================================= S11
section('S11 union bounds (random-flow forest lemma, any inclusion probability q)')
u10, u11 = L.union_bound(10, 0.5), L.union_bound(11, 0.5)
ok(u10 > 1 > u11, 'q=1/2 reproduces the paper: U(10)=%.3f > 1 > U(11)=%.4f' % (u10, u11))
UB = {11: (511, 0.23421), 12: (576, 0.13357), 13: (667, 0.07666), 14: (781, 0.045528),
      15: (910, 0.02648), 16: (1056, 0.015496), 17: (1221, 0.009052), 18: (1408, 0.00519906),
      20: (1837, 0.00170376), 24: (2969, 0.000173393), 28: (4537, 1.66387e-05), 32: (6638, 1.52675e-06)}
for d, (M, q) in UB.items():
    M2, u = L.size_bound(d, q)
    ok(M2 is not None and M2 <= M and M < 2 ** (d - 1),
       'sparse random landmarks q=%.6g: U(%d)=%.3f<1 and P(|S|>%d)<1-U  => edim_m(Q%d) <= %d  (ln M / d^(1/3) = %.3f; density 1/2: ~2^%d)' % (q, d, u, M, d, M, log(M) / d ** (1 / 3), d - 1))

if CLEANUP:
    import shutil; shutil.rmtree(BUILD, ignore_errors=True)
print('\nTotal checks: %d   (elapsed %.0f s)' % (NCHECK, time.time() - T0))
print('ALL CHECKS PASSED')
