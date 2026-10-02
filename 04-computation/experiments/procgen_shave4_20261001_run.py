#!/usr/bin/env python3
"""procgen_shave4_20261001_run.py -- shave4 lane (collatz-procgen-20260922), 2026-10-01.

Re-verifies every computational claim of
    05-knowledge/results/procgen_shave4_20261001_redei_graphs.md
(Redei graphs: parity of embedding counts of spanning oriented graphs in tournaments; THM-4526(A) audit;
u(9) = 14).  Prints only to stdout and ends with ALL CHECKS PASSED.

Needs: cc (C99), nauty (gentourng, geng, directg, countg, labelg, pickg), python3.  Single-threaded, < 300 MB.
Build products and class files go to scratch/procgen_shave4/build/ (never committed).
Usage:  python3 procgen_shave4_20261001_run.py [--skip-u9] [--with-r9] [--with-k14]
   --skip-u9 : do not re-run the exhaustive K = 15 search on 9 vertices (prints the recorded result instead)
   --with-k14: also re-run the enumeration of all 14-arc shavings on 9 vertices (54 classes)
   --with-r9 : also re-run the multi-depth Redei search on 9 vertices (long; by default the recorded
               result 05-knowledge/results/procgen_shave4_20261001_r9search.out is re-checked directly)
"""
import os, sys, subprocess, itertools, random, collections, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
BUILD = os.path.join(ROOT, 'scratch', 'procgen_shave4', 'build')
os.makedirs(BUILD, exist_ok=True)
T0 = time.time()
NCHECK = 0

def check(cond, msg):
    global NCHECK
    if not cond:
        print("FAIL", msg); sys.stdout.flush(); sys.exit(1)
    NCHECK += 1
    print("PASS", msg); sys.stdout.flush()

def sh(cmd, inp=None):
    r = subprocess.run(cmd, shell=True, input=inp, capture_output=True, text=True, cwd=BUILD)
    if r.returncode != 0:
        print("FAIL command:", cmd, r.stderr[:500]); sys.exit(1)
    return r.stdout

def section(t):
    print("\n=== %s  [%.0fs]" % (t, time.time() - T0)); sys.stdout.flush()

# ---------------------------------------------------------------- build
section("0. build")
for src, exe in [("procgen_shave4_20261001_engine.c", "engine"), ("procgen_shave4_20261001_closure.c", "closure")]:
    sh("cc -O2 -o %s %s" % (exe, os.path.join(HERE, src)))
for n in (7, 8, 9):
    sh("cc -O2 -DN=%d -o u9_n%d %s" % (n, n, os.path.join(HERE, "procgen_shave4_20261001_u9.c")))
for n in range(2, 10):
    p = os.path.join(BUILD, "cls%d.txt" % n)
    if not os.path.exists(p):
        sh("gentourng -q %d > cls%d.txt" % (n, n))
for n in range(1, 8):
    p = os.path.join(BUILD, "og%d.txt" % n)
    if not os.path.exists(p):
        sh("geng -q %d | directg -q -o -G > og%d.txt" % (n, n))
CLS = {n: [l.strip() for l in open(os.path.join(BUILD, "cls%d.txt" % n))] for n in range(2, 10)}
check([len(CLS[n]) for n in range(2, 10)] == [1, 2, 4, 12, 56, 456, 6880, 191536], "tournament class counts A000568(2..9)")

# ---------------------------------------------------------------- python helpers
def adj_of(s, n):
    A = [[0] * n for _ in range(n)]; b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[b] == '1': A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A

def emb_count(n, arcs, A):
    """exact number of bijections pi with pi(S) inside T (plain backtracking)"""
    out = [[] for _ in range(n)]; inn = [[] for _ in range(n)]
    for u, v in arcs: out[u].append(v); inn[v].append(u)
    img = [-1] * n; cnt = 0
    def rec(k, used):
        nonlocal cnt
        if k == n: cnt += 1; return
        for y in range(n):
            if used >> y & 1: continue
            ok = True
            for v in out[k]:
                if v < k and not A[y][img[v]]: ok = False; break
            if ok:
                for u in inn[k]:
                    if u < k and not A[img[u]][y]: ok = False; break
            if ok:
                img[k] = y; rec(k + 1, used | 1 << y); img[k] = -1
    rec(0, 0)
    return cnt

def phi_vector(n, arcs):
    return [emb_count(n, arcs, adj_of(s, n)) & 1 for s in CLS[n]]

def has_involution(n, arcs):
    E = set(tuple(sorted(e)) for e in arcs)
    deg = [0] * n
    for a, b in E: deg[a] += 1; deg[b] += 1
    perm = [-1] * n
    def rec(i, nontriv):
        while i < n and perm[i] != -1: i += 1
        if i == n:
            return nontriv and all(tuple(sorted((perm[a], perm[b]))) in E for a, b in E)
        perm[i] = i
        if rec(i + 1, nontriv): return True
        perm[i] = -1
        for j in range(i + 1, n):
            if perm[j] == -1 and deg[j] == deg[i]:
                perm[i] = j; perm[j] = i
                if rec(i + 1, True): return True
                perm[i] = -1; perm[j] = -1
        return False
    return rec(0, False)

def is_automorphism(perm, arcs):
    S = set(arcs)
    return all((perm[u], perm[v]) in S for u, v in arcs)

def parse_arcs(tokens):
    return [tuple(map(int, t.split('>'))) for t in tokens]

# ---------------------------------------------------------------- 1. THM-4526(A)
section("1. THM-4526(A) re-verified: H_n = path + first->last arc")
expect = {3: (2, 1, 1, 1), 4: (4, 0, 4, 0), 5: (12, 1, 8, 4), 6: (56, 0, 56, 0), 7: (456, 0, 254, 202),
          8: (6880, 0, 6880, 0), 9: (191536, 0, 100178, 91358)}
zeros = {}
for n in range(3, 10):
    out = sh("./engine hp < cls%d.txt" % n)
    summ = [l for l in out.split('\n') if l.startswith('HPSUMMARY')][0]
    kv = dict(x.split('=') for x in summ.split()[1:])
    got = (int(kv['classes']), int(kv['zero']), int(kv['odd']), int(kv['even']))
    zeros[n] = [l.split()[1] for l in out.split('\n') if l.startswith('ZERO')]
    check(got == expect[n] and kv['min_nonzero'] == '1',
          "n=%d: classes=%d, classes with no copy of H_n=%d, odd/even counts %d/%d (closing HPs = n*hc checked per class)" % ((n,) + got))
check(zeros[3] == ['101'] and zeros[5] == ['1110101111'], "the only avoiders are C3 (n=3) and C3[1,C3,1] = 1110101111 (n=5)")
A5 = adj_of('1110101111', 5)
blk = [1, 2, 3]
check(all(A5[0][x] for x in blk) and all(A5[x][4] for x in blk) and A5[4][0] and
      A5[1][2] and A5[2][3] and A5[3][1], "1110101111 is a => C3 => b -> a, i.e. C3[1,C3,1]")

# ---------------------------------------------------------------- 2. path + chords
section("2. Redei graphs of the form path + chords (R_n)")
R_expect = {3: [[]], 4: [[], ['0-3']], 5: [[], ['0-3'], ['1-4'], ['0-3', '1-4']],
            6: [[], ['0-3'], ['0-5'], ['2-5'], ['0-5', '1-4'], ['0-3', '2-5'], ['0-3', '1-4', '2-5']],
            7: [[], ['0-3', '1-6'], ['0-5', '1-6'], ['0-3', '3-6'], ['0-5', '3-6']],
            8: [[], ['0-7']]}
def norm(L): return sorted(tuple(sorted(x)) for x in L)
for n in range(3, 9):
    out = sh("./engine til %d < cls%d.txt" % (n, n))
    til = [l for l in out.split('\n') if l.startswith('TIL ')][0]
    kv = dict(x.split('=') for x in til.split()[1:] if '=' in x)
    check(kv['unassigned'] == '0' and kv['conflicts'] == '0' and kv['nonuniform'] == '0' and
          int(kv['odd_tiling_classes']) == len(CLS[n]) and int(float(kv['orbit_sum'])) == 2 ** (n * (n - 1) // 2),
          "n=%d: tilings partition {0,1}^C(n-1,2) into the %d classes, every class has an odd number of tilings (Redei), "
          "multiplicities = |Aut|, orbit sum = 2^C(n,2)" % (n, len(CLS[n])))
    R = [l.split()[2:] for l in out.split('\n') if l.startswith('R ')]
    check(norm(R) == norm(R_expect[n]), "n=%d: R_n = %s" % (n, ['{' + ','.join(c) + '}' for c in R_expect[n]]))
out = sh("./engine chords 9 0-7 1-8 0-8 < cls9.txt")
lines = dict((l.split()[0], l.split()[1:]) for l in out.split('\n') if l)
check(lines['SINGLE'] == [] and lines['PAIRS'] == ['{0-7,1-8}'] and lines['CANDIDATES'] == ['{}', '{0-7,1-8}'],
      "n=9 (all 191536 classes): no single chord, the only chord pair is {07,18}; among subsets of {07,18,08} only {} and {07,18}")

# interval heredity: candidates for n = 10, 11
def heredity_candidates(n, Rsmall):
    chords = [(i, j) for i in range(n) for j in range(i + 2, n)]
    allowed_by_len = {}
    for L, sets in Rsmall.items():
        allowed_by_len[L] = [set(tuple(map(int, c.split('-'))) for c in S) for S in sets]
    cands = []
    # restrict search: chords that survive single interval tests
    pool = []
    for (i, j) in chords:
        ok = True
        for L, sets in allowed_by_len.items():
            for k in range(0, n - L + 1):
                if k <= i and j <= k + L - 1:
                    sh_ = (i - k, j - k)
                    if not any(sh_ in S for S in sets): ok = False
        if ok: pool.append((i, j))
    for r in range(len(pool) + 1):
        for C in itertools.combinations(pool, r):
            Cs = set(C); ok = True
            for L, sets in allowed_by_len.items():
                for k in range(0, n - L + 1):
                    sub = set((i - k, j - k) for (i, j) in Cs if k <= i and j <= k + L - 1)
                    if sub not in sets: ok = False; break
                if not ok: break
            if ok: cands.append(sorted(Cs))
    return cands
Rk = {L: R_expect[L] for L in range(4, 9)}
Rk[9] = [[], ['0-7', '1-8']]
cand10 = heredity_candidates(10, {8: Rk[8], 9: Rk[9]})
check(sorted(map(tuple, cand10)) == sorted([(), ((0, 9),), ((0, 7), (1, 8), (2, 9)), ((0, 7), (0, 9), (1, 8), (2, 9))]),
      "n=10: interval heredity (B5) leaves exactly the candidates {}, {09}, {07,18,29}, {07,18,29,09}")
Rk[10] = [[], ['0-9']]
cand11 = heredity_candidates(11, {9: Rk[9], 10: Rk[10]})
check(sorted(map(tuple, cand11)) == sorted(tuple(sorted(c)) for r in range(4) for c in itertools.combinations([(0, 9), (1, 10), (0, 10)], r)),
      "n=11: interval heredity leaves exactly the 8 subsets of {09, 1-10, 0-10}")

def count_hp_chords(s, n, chords):
    A = adj_of(s, n); need = collections.defaultdict(list)
    for i, j in chords: need[j].append(i)
    v = [0] * n; cnt = 0
    def rec(d, used):
        nonlocal cnt
        if d == n: cnt += 1; return
        for w in range(n):
            if used >> w & 1: continue
            if d > 0 and not A[v[d - 1]][w]: continue
            if any(not A[v[i]][w] for i in need[d]): continue
            v[d] = w; rec(d + 1, used | 1 << w)
    rec(0, 0)
    return cnt

refute = {10: [[(0, 7), (1, 8), (2, 9)], [(0, 7), (1, 8), (2, 9), (0, 9)]],
          11: [[(0, 9)], [(1, 10)], [(0, 10)], [(0, 9), (0, 10)], [(1, 10), (0, 10)], [(0, 9), (1, 10), (0, 10)]]}
for n, sets in refute.items():
    for C in sets:
        out = sh("./engine witness %d 7 400 %s" % (n, ' '.join('%d-%d' % c for c in C)))
        tok = out.split()
        check(tok[0] == 'WITNESS', "n=%d chords %s: random search finds a tournament with an even count" % (n, C))
        tour = [t for t in tok if t.startswith('tour=')][0][5:]
        c_py = count_hp_chords(tour, n, C)
        check(c_py % 2 == 0 and str(c_py) == [t for t in tok if t.startswith('count=')][0][6:],
              "   independent Python count confirms: %d (even) embeddings in %s" % (c_py, tour))
out = sh("./engine witness 11 3 300 0-9 1-10")
check(out.startswith('NOWITNESS'), "D_11 = P_11 + {09,1-10}: 300 random 11-tournaments, all counts odd (EMPIRICAL support; PROVED by D1)")
out = sh("./engine witness 13 5 60 0-11 1-12")
check(out.startswith('NOWITNESS'), "D_13: 60 random 13-tournaments, all counts odd (EMPIRICAL support; PROVED by D1)")

# ---------------------------------------------------------------- 3. D1 identities
section("3. Theorem D1: exact identities behind the proof")
exp = {5: (12, 12, None), 6: (56, 32, 56), 7: (456, 456, None), 8: (6880, 3588, 6880), 9: (191536, 191536, None)}
for n in range(5, 10):
    out = sh("./engine dcheck %d < cls%d.txt" % (n, n))
    kv = dict(x.split('=') for x in out.split()[1:])
    ok = kv['identity_failures'] == '0' and int(kv['classes']) == exp[n][0] and int(kv['emb_odd']) == exp[n][1]
    if exp[n][2] is not None: ok = ok and int(kv['matches_1+sum_hc(T-w)']) == exp[n][2]
    check(ok, "n=%d: emb(D_n)=H-A-B+AB, A=sum_w hc(T-w)d-(w), B=sum_w hc(T-w)d+(w), AB even, for every class; "
              "emb(D_n) odd in %s classes%s" % (n, kv['emb_odd'], "" if n % 2 else "; = 1+sum_w hc(T-w) mod 2 for every class"))

# ---------------------------------------------------------------- 4. census
section("4. census of parity-rigid and Redei graphs")
dag_expect = {2: 1, 3: 2, 4: 5, 5: 21, 6: 43, 7: 156}
redei_lists = {}
for n in range(2, 8):
    out = sh("geng -q %d | directg -q -a -G | ./engine dag %d cls%d.txt" % (n, n, n))
    summ = [l for l in out.split('\n') if l.startswith('DAGSUMMARY')][0]
    kv = dict(x.split('=') for x in summ.split()[1:])
    redei_lists[n] = [l for l in out.split('\n') if l.startswith('REDEI')]
    check(int(kv['redei']) == dag_expect[n], "n=%d: %s DAG classes, %s Redei (Gray-code completion parity, full class table)" % (n, kv['dags'], kv['redei']))
    rig = [l for l in out.split('\n') if l.startswith('REDEI') or l.startswith('EVEN0')]
    bad = 0
    for l in rig:
        w = l.split(); arcs = parse_arcs(w[w.index('arcs') + 1:])
        if arcs and not has_involution(n, arcs): bad += 1
    check(bad == 0, "   B4: every rigid DAG with |Aut S| odd and >= 1 arc has an involution of its underlying graph (%d checked)" % len(rig))
    bad = sum(1 for l in redei_lists[n] if int(l.split('ext=')[1].split()[0]) % 2 == 0)
    check(bad == 0, "   B4: every Redei graph has an odd number of linear extensions")
rmax = {n: max(int(l.split('m=')[1].split()[0]) for l in redei_lists[n]) for n in range(2, 8)}
check([rmax[n] for n in range(2, 8)] == [1, 2, 4, 6, 8, 9], "r(n) = max arcs of a Redei graph = 1,2,4,6,8,9 for n = 2..7 (= u(n))")
for n in range(5, 8):
    out = sh("geng -q %d | directg -q -a -G | ./engine redei %d cls%d.txt" % (n, n, n))
    L1 = sorted(l.replace(' iso=%s' % l.split('iso=')[1].split()[0], '') for l in redei_lists[n])
    L2 = sorted(l for l in out.split('\n') if l.startswith('REDEI'))
    check(L1 == L2, "n=%d: backtracking tester with the B4 filters returns the same %d Redei classes" % (n, len(L2)))
out = sh("geng -q 8 | directg -q -a -G | ./engine redei 8 cls8.txt")
L8 = [l for l in out.split('\n') if l.startswith('REDEI')]
m8 = collections.Counter(int(l.split('m=')[1].split()[0]) for l in L8)
check(len(L8) == 220 and max(m8) == 10 and m8[10] == 2, "n=8: all 20286025 DAG classes: 220 Redei, max 10 arcs (2 classes); r(8) = 10 < u(8) = 11; histogram %s" % sorted(m8.items()))
def d6_gen(n, arcs):
    bits = [0] * (n * n)
    for u, v in arcs: bits[u * n + v] = 1
    while len(bits) % 6: bits.append(0)
    o = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        o += chr(v + 63)
    return o
def canon_list(n, L):
    return sh("labelg -q", inp='\n'.join(d6_gen(n, a) for a in L) + '\n').split() if L else []
alls = {n: [parse_arcs(l.split('arcs ')[1].split()) for l in redei_lists[n]] for n in range(4, 8)}
alls[8] = [parse_arcs(l.split('arcs ')[1].split()) for l in L8]
tested = bad = 0
for n in range(5, 9):
    small = set(canon_list(n - 1, alls[n - 1])); subs = []
    for S_ in alls[n]:
        for kind in (0, 1):
            cand = [v for v in range(n) if not any((b if kind == 0 else a) == v for a, b in S_)]
            if len(cand) != 1: continue
            s = cand[0]; mp = {v: i for i, v in enumerate([x for x in range(n) if x != s])}
            subs.append([(mp[a], mp[b]) for a, b in S_ if s not in (a, b)])
    cs = canon_list(n - 1, subs); tested += len(subs); bad += sum(1 for c in cs if c not in small)
check(tested == 138 and bad == 0, "ideal formula corollary: deleting the unique source (or sink) of a Redei graph leaves a Redei graph (%d cases, n = 5..8)" % tested)
nine = [l for l in redei_lists[7] if ' m=9 ' in l]
check(len(nine) == 2, "n=7: the two 9-arc Redei graphs: " + " | ".join(l.split('arcs ')[1] for l in nine))
ten = [l for l in L8 if ' m=10 ' in l]
H4c = canon_list(4, [[(0, 1), (1, 2), (2, 3), (0, 3)]])[0]
for l in ten:
    arcs = parse_arcs(l.split('arcs ')[1].split()); Sset = set(arcs); E = set(frozenset(a) for a in arcs)
    found_tp = None
    for p_ in itertools.permutations(range(8)):
        if any(p_[p_[i]] != i for i in range(8)) or all(p_[i] == i for i in range(8)): continue
        if not all(frozenset((p_[u], p_[v])) in E for u, v in arcs): continue
        rev = [(u, v) for u, v in arcs if (p_[u], p_[v]) not in Sset]
        if len(rev) == 2: found_tp = (p_, rev); break
    ok = found_tp is not None
    if ok:
        core = [a for a in arcs if a not in found_tp[1]]
        comp = []; seen = set()
        for st in range(8):
            if st in seen: continue
            cc = {st}; fr = [st]
            while fr:
                x = fr.pop()
                for a, b in core:
                    for y in ((b,) if a == x else (a,) if b == x else ()):
                        if y not in cc: cc.add(y); fr.append(y)
            seen |= cc; comp.append(sorted(cc))
        ok = sorted(len(c) for c in comp) == [4, 4] and all(
            canon_list(4, [[(c.index(a), c.index(b)) for a, b in core if a in c]])[0] == H4c for c in comp)
    check(ok, "n=8 10-arc Redei graph %s = (H_4 u H_4) + a twisted pair %s under the involution %s" %
          (l.split('arcs ')[1], found_tp[1] if found_tp else None, [(i, found_tp[0][i]) for i in range(8) if found_tp and found_tp[0][i] > i]))

# ---------------------------------------------------------------- 5. closure
section("5. deletion-reversal closure (Horn) and the D69 derivations")
open(os.path.join(BUILD, 'targets4.txt'), 'w').write("0>1 1>2 2>3 0>3\n")
open(os.path.join(BUILD, 'targets5.txt'), 'w').write("0>1 1>2 2>3 3>4 0>3 1>4\n")
open(os.path.join(BUILD, 'targets6.txt'), 'w').write("0>1 1>2 2>3 3>4 4>5 0>3 1>4 2>5\n")
out = sh("./closure 7")
stats = [l for l in out.split('\n') if l.startswith('n=')]
exp_stats = {1: (1, 1, 1, 1, 1, 0), 2: (2, 2, 1, 2, 1, 0), 3: (7, 5, 2, 5, 2, 0), 4: (42, 23, 5, 22, 5, 1),
             5: (582, 178, 21, 174, 20, 4), 6: (21480, 2999, 43, 2976, 43, 23), 7: (2142288, 130574, 156, 130345, 151, 229)}
for l in stats:
    n = int(l.split()[0][2:])
    kv = dict(x.split('=') for x in l.split()[1:] if '=' in x and not x.startswith('(redei'))
    rig = int(l.split('rigid=')[1].split()[0]); red = int(l.split('(redei=')[1].split(')')[0])
    expl = int(l.split('explained=')[1].split()[0]); expr = int(l.split('(redei ')[1].split(')')[0])
    unex = int(l.split('unexplained=')[1].split()[0])
    got = (int(l.split('classes=')[1].split()[0]), rig, red, expl, expr, unex)
    check(got == exp_stats[n] and ' wrong=0 ' in l and 'contradictions=0' in l,
          "n=%d: %d oriented graphs, %d parity-rigid (%d Redei); Horn closure certifies %d (%d Redei), %d not certified" % ((n,) + got))
targ = out[out.index('TARGET n=4'):]
print("   derivations of path + span-3 (n = 4, 5, 6) printed by the closure:")
for l in targ.split('\n'):
    if l.startswith('TARGET') or l.startswith('  '): print("   " + l)
check(all(('TARGET n=%d' % n) in out for n in (4, 5, 6)) and out.count('true status 1, known 1') == 3,
      "path + span-3 is certified Redei by the closure for n = 4, 5, 6")

# ---------------------------------------------------------------- 6. D69 hand proof
section("6. D69 hand proof (n = 6) and the n = 7 obstruction")
P6 = [(i, i + 1) for i in range(5)]
S = P6 + [(0, 3), (1, 4), (2, 5)]
revg = [a for a in S if a != (1, 4)] + [(4, 1)]
check(is_automorphism([0, 3, 4, 1, 2, 5], revg), "step 1: iota = (1 3)(2 4) is an automorphism of S with 1>4 reversed (twisted pair)")
phiS = phi_vector(6, S); phiS2 = phi_vector(6, P6 + [(0, 3), (2, 5)])
check(phiS == [1] * 56 and phiS2 == [1] * 56, "path + span-3 and P_6 + {03,25} are Redei on all 56 classes")
okA = okAB = okIE = okB = True
for s in CLS[6]:
    A = adj_of(s, 6); Aop = [[A[j][i] for j in range(6)] for i in range(6)]
    def counts(M):
        H = a = b = ab = e = 0
        for p in itertools.permutations(range(6)):
            if all(M[p[i]][p[i + 1]] for i in range(5)):
                H += 1; x = M[p[3]][p[0]]; y = M[p[5]][p[2]]
                a += x; b += y; ab += x and y; e += (not x) and (not y)
        return H, a, b, ab, e
    H, a, b, ab, e = counts(A); _, aop, _, _, _ = counts(Aop)
    okIE &= (e == H - a - b + ab); okAB &= (ab % 2 == 0); okA &= (a % 2 == 0); okB &= (b == aop)
check(okIE and okAB and okA and okB, "step 2: emb = H - A - B + AB with AB even, A even, B(T) = A(T^op) on every 6-class")
c0, c1, c2, c3, w1, w2 = range(6)
C4tail = [(c0, c1), (c1, c2), (c2, c3), (c3, c0), (c3, w1), (w1, w2)]
T1 = [(c1, c2), (c2, c3), (c3, c0), (c3, w1), (w1, w2)]
Y2 = [(c1, c0), (c1, c2), (c2, c3), (c3, c0), (c3, w1), (w1, w2)]
check(phi_vector(6, C4tail) == [0] * 56, "step 3: phi(C_4 + tail of length 2) = 0 on all classes (so A is even)")
Y2g = [a for a in Y2 if a != (c3, c0)] + [(c0, c3)]
check(is_automorphism([c2, c1, c0, c3, w1, w2], Y2g), "   (c0 c2) is an automorphism of Y_2 with c3>c0 reversed")
path6 = [(c1, c0), (c1, c2), (c2, c3), (c3, w1), (w1, w2)]
check(phi_vector(6, Y2) == [1] * 56 and phi_vector(6, path6) == [1] * 56, "   phi(Y_2) = phi(c0<c1>c2>c3>w1>w2) = 1 (Forcade: compositions (6),(1,5))")
T1p = [(c1, c2), (c2, c3), (c0, c3), (c3, w1), (w1, w2)]
T1pp = [(c2, c1), (c2, c3), (c0, c3), (c3, w1), (w1, w2)]
check(is_automorphism([c2, c1, c0, c3, w1, w2], [(c2, c3), (c0, c3), (c3, w1), (w1, w2)]), "   T_1' minus c1>c2 has the involution (c0 c2)")
T1ppg = [(c2, c1), (c2, c3), (c0, c3), (w1, c3), (w1, w2)]
check(is_automorphism([c0, w2, w1, c3, c2, c1], T1ppg), "   (c2 w1)(c1 w2) is an automorphism of T_1'' with c3>w1 reversed")
fin = [(c2, c1), (c2, c3), (c0, c3), (w1, w2)]
check(phi_vector(6, T1) == [1] * 56 and phi_vector(6, T1p) == [1] * 56 and phi_vector(6, T1pp) == [1] * 56 and phi_vector(6, fin) == [1] * 56,
      "   phi(T_1) = phi(T_1') = phi(T_1'') = phi(antidirected P_4 + P_2) = 1; so phi(C_4 + tail) = 1 + 1 = 0")
# n = 7
c0, c1, c2, c3, w1, w2, w3 = range(7)
T1_7 = [(c1, c2), (c2, c3), (c3, c0), (c3, w1), (w1, w2), (w2, w3)]
Y3 = [(c1, c0), (c1, c2), (c2, c3), (c3, w1), (w1, w2), (w2, w3)]
H4P3 = [(i, i + 1) for i in range(6)] + [(0, 3)]
check(not has_involution(7, T1_7), "n=7: the tree T_1 of the same expansion is S(1,2,3), the smallest asymmetric tree (no involution)")
vT1 = phi_vector(7, T1_7); vY3 = phi_vector(7, Y3); vH = phi_vector(7, H4P3)
check(len(set(vT1)) == 2 and vY3 == [0] * 456 and all(h == (1 + t + y) % 2 for h, t, y in zip(vH, vT1, vY3)) and len(set(vH)) == 2,
      "n=7: phi(H_4.P_3) = 1 + phi(T_1) + phi(Y_3) with phi(Y_3) = 0 and phi(T_1) non-constant (%d of 456 classes even)" % vT1.count(0))
S7 = [(i, i + 1) for i in range(6)] + [(i, i + 3) for i in range(4)]
P7avoid = sh("./engine embedall 7 cls7.txt %s" % ' '.join('%d\\>%d' % a for a in S7))
check('avoiders=2' in P7avoid, "path + span-3 at n = 7 is avoided by exactly 2 classes (count 0 is even): not Redei, and by B5 not for any n >= 7")
def tt_flip(n, flips):
    A = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
    for i, j in flips: A[i][j], A[j][i] = A[j][i], A[i][j]
    return A
pairs7 = list(itertools.combinations(range(7), 2))
single_odd = all(emb_count(7, S7, tt_flip(7, [p])) % 2 == 1 for p in pairs7)
dbl_even = [pp for pp in itertools.combinations(pairs7, 2) if emb_count(7, S7, tt_flip(7, list(pp))) % 2 == 0]
check(emb_count(7, S7, tt_flip(7, [])) == 1 and single_odd and len(dbl_even) == 13 and emb_count(7, S7, tt_flip(7, [(0, 6), (3, 6)])) == 2,
      "near-transitive witness: TT_7 with 0>6 and 3>6 reversed has exactly 2 copies; every single reversal keeps the count odd; "
      "13 double reversals make it even (= the 13 degree-2 terms of the ANF of phi - 1)")

# ---------------------------------------------------------------- 7. trees
section("7. oriented trees: Theorem T (hereditarily symmetric => rigid)")
def tree_has_S123(n, edges):
    nb = collections.defaultdict(list)
    for a, b in edges: nb[a].append(b); nb[b].append(a)
    def depth(x, parent):
        return 1 + max([depth(y, x) for y in nb[x] if y != parent], default=0)
    for c in range(n):
        d = sorted((depth(y, c) for y in nb[c]), reverse=True)
        if len(d) >= 3 and d[0] >= 3 and d[1] >= 2 and d[2] >= 1: return True
    return False
def hereditarily_symmetric(n, edges):
    nb = collections.defaultdict(set)
    for a, b in edges: nb[a].add(b); nb[b].add(a)
    for mask in range(1, 1 << n):
        vs = [v for v in range(n) if mask >> v & 1]
        if len(vs) < 2: continue
        # connected?
        seen = {vs[0]}; st = [vs[0]]
        while st:
            x = st.pop()
            for y in nb[x]:
                if mask >> y & 1 and y not in seen: seen.add(y); st.append(y)
        if len(seen) != len(vs): continue
        idx = {v: i for i, v in enumerate(vs)}
        sub = [(idx[a], idx[b]) for a, b in edges if mask >> a & 1 and mask >> b & 1]
        if not has_involution(len(vs), sub): return False
    return True
tree_exp = {3: (3, 3, 0), 4: (8, 8, 0), 5: (27, 27, 0), 6: (91, 91, 0), 7: (350, 286, 64), 8: (1376, 1056, 320)}
for n in range(3, 9):
    out = sh("geng -q -c %d %d:%d | directg -q -o -G | ./engine rigid %d cls%d.txt" % (n, n - 1, n - 1, n, n))
    rows = [l for l in out.split('\n') if l.startswith('RIGID') or l.startswith('NONRIGID')]
    nr = sum(1 for l in rows if l.startswith('RIGID'))
    check((len(rows), nr, len(rows) - nr) == tree_exp[n], "n=%d: %d oriented trees, %d parity-rigid, %d not" % ((n,) + tree_exp[n]))
    cache = {}; bad1 = 0; paths = 0; badF = 0; tree_nonrigid = collections.defaultdict(bool)
    for l in rows:
        w = l.split(); arcs = parse_arcs(w[w.index('arcs') + 1:]); E = tuple(sorted(tuple(sorted(a)) for a in arcs))
        if E not in cache: cache[E] = (hereditarily_symmetric(n, list(E)), tree_has_S123(n, list(E)))
        hs, s123 = cache[E]
        rigid = l.startswith('RIGID')
        if hs and not rigid: bad1 += 1
        tree_nonrigid[E] |= (not rigid)
        deg = collections.Counter(x for a in arcs for x in a)
        if max(deg.values()) <= 2:   # oriented Hamiltonian path: Forcade formula
            paths += 1
            end = [v for v in range(n) if deg[v] == 1][0]; order = [end]; prev = -1
            while len(order) < n:
                x = order[-1]; nxt = [y for a in arcs for y in a if x in a and y != x and y != prev]
                nxt = [y for y in nxt if y not in order]; prev = x; order.append(nxt[0])
            pos = {v: i for i, v in enumerate(order)}
            Bpos = [min(pos[u], pos[v]) + 1 for u, v in arcs if pos[u] > pos[v]]
            val = 0
            for r in range(len(Bpos) + 1):
                for K in itertools.combinations(sorted(Bpos), r):
                    cuts = [0] + list(K) + [n]; parts = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
                    acc = 0; cf = True
                    for p_ in parts:
                        if acc & p_: cf = False
                        acc |= p_
                    val ^= cf
            if not (rigid and int(l.split('val=')[1].split()[0]) == val): badF += 1
    check(bad1 == 0, "   every orientation of every hereditarily symmetric tree is rigid (Theorem T)")
    bad2 = sum(1 for E in cache if tree_nonrigid[E] != cache[E][1])
    check(bad2 == 0, "   a tree has a non-rigid orientation iff it contains the spider S(1,2,3) (%d trees; %d such)" % (len(cache), sum(1 for E in cache if cache[E][1])))
    check(badF == 0, "   all %d oriented Hamiltonian paths are rigid with value = #{carry-free cut compositions} mod 2 (Forcade)" % paths)

section("7b. cycles (Corollary G2) and doubling (D4)")
for n in range(4, 9):
    cyc = []
    for line in sh("geng -q -c %d %d:%d" % (n, n, n)).split():
        N_ = ord(line[0]) - 63; bits = []
        for ch in line[1:]:
            v = ord(ch) - 63
            for k in range(5, -1, -1): bits.append((v >> k) & 1)
        deg = [0] * N_; b = 0
        for j in range(1, N_):
            for i in range(j):
                if bits[b]: deg[i] += 1; deg[j] += 1
                b += 1
        if all(d == 2 for d in deg): cyc.append(line)
    out = sh("directg -q -o -G | ./engine rigid %d cls%d.txt" % (n, n), inp='\n'.join(cyc) + '\n')
    rows = [l for l in out.split('\n') if l.startswith('RIGID') or l.startswith('NONRIGID')]
    nr = sum(1 for l in rows if l.startswith('RIGID'))
    check(len(cyc) == 1 and (nr == len(rows) if n % 2 == 0 else nr == 0),
          "C_%d: %d orientation classes, %s parity-rigid" % (n, len(rows), "all" if n % 2 == 0 else "none"))
c4 = [parse_arcs(l.split('arcs ')[1].split()) for l in redei_lists[4]]
tests = [S + [(u + 4, v + 4) for u, v in S] + [(x, x + 4)] for S in c4 for x in range(4)]
def d6_of8(arcs):
    bits = [0] * 64
    for u, v in arcs: bits[u * 8 + v] = 1
    o = '&' + chr(8 + 63)
    for k in range(0, 64, 6):
        v = 0
        for b in (bits[k:k + 6] + [0] * 6)[:6]: v = 2 * v + b
        o += chr(v + 63)
    return o
can8 = set(sh("labelg -q", inp='\n'.join(d6_of8(parse_arcs(l.split('arcs ')[1].split())) for l in L8) + '\n').split())
cant = sh("labelg -q", inp='\n'.join(d6_of8(a) for a in tests) + '\n').split()
check(len(tests) == 20 and all(c in can8 for c in cant), "doubling: S u S + (x,0)->(x,1) is Redei for all 5 Redei S on 4 vertices and all x (20/20 in the n=8 census)")

# ---------------------------------------------------------------- 8. halving lemma
section("8. halving lemma: emb(K + (u->v)) = emb(K)/2 when an involution of K swaps u, v")
rng = random.Random(20261001); trials = 0; good = 0
for n in (6, 7):
    iota = list(range(n)); iota[0], iota[1] = 1, 0; iota[2], iota[3] = 3, 2   # (0 1)(2 3)
    for t in range(15):
        K = set()
        for _ in range(rng.randint(3, 7)):
            u, v = rng.sample(range(n), 2)
            if {u, v} == {0, 1}: continue
            a, b = (u, v), (iota[u], iota[v])
            if (v, u) in K or (b[1], b[0]) in K or a == (b[1], b[0]): continue
            K.add(a); K.add(b)
        K = sorted(K)
        if not is_automorphism(iota, K): continue
        Kf = K + [(0, 1)]
        for s in rng.sample(CLS[n], 6):
            Am = adj_of(s, n)
            e1, e2 = emb_count(n, K, Am), emb_count(n, Kf, Am)
            trials += 1; good += (e1 == 2 * e2)
check(trials > 50 and good == trials, "exact equality emb(K)=2 emb(K+(0>1)) in %d random (K, T) pairs" % trials)

# ---------------------------------------------------------------- 9. u(n), n <= 7, independent of S15
section("9. u(n) for n <= 7 (forward arc sets, full completion test)")
def unav(n, e):
    out = sh("./engine unav %d cls%d.txt %d" % (n, n, e))
    summ = [l for l in out.split('\n') if l.startswith('UNAV')][0]
    return int(summ.split('shavings=')[1]), [int(l.split()[1]) for l in out.split('\n') if l.startswith('U ')]
for n, u in [(3, 2), (4, 4), (5, 6), (6, 8), (7, 9)]:
    a, _ = unav(n, u + 1); b, masks = unav(n, u)
    check(a == 0 and b > 0, "n=%d: no forward %d-arc shaving, %d forward %d-arc shavings: u(%d) = %d" % (n, u + 1, b, u, n, u))
    if n == 7:
        pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
        perms = list(itertools.permutations(range(7)))
        canon = set()
        for m in masks:
            arcs = [pairs[k] for k in range(21) if m >> k & 1]
            best = min(tuple(sorted((p[u_], p[v_]) for u_, v_ in arcs)) for p in perms)
            canon.add(best)
        check(len(canon) == 51, "n=7: the %d forward 9-arc shavings form 51 isomorphism classes (as S15 found)" % len(masks))
        red9 = set()
        for l in nine:
            arcs = parse_arcs(l.split('arcs ')[1].split())
            red9.add(min(tuple(sorted((p[u_], p[v_]) for u_, v_ in arcs)) for p in perms))
        check(red9 <= canon, "n=7: both 9-arc Redei graphs are among the 51 maximum shavings")

# ---------------------------------------------------------------- 10. u(9) = 14, r(n)
section("10. D70: u(9) = 14 (orderly search over a symmetric host), and r(n)")
def d6_of(n, arcs):
    bits = [0] * (n * n)
    for u, v in arcs: bits[u * n + v] = 1
    while len(bits) % 6: bits.append(0)
    o = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        o += chr(v + 63)
    return o
def classes_of(n, sets):
    if not sets: return set()
    return set(sh("labelg -q", inp='\n'.join(d6_of(n, a) for a in sets) + '\n').split())
def found_sets(out, tag):
    return [[tuple(map(int, a.split('>'))) for a in l.split('arcs:')[1].split()] for l in out.split('\n') if l.startswith(tag + ' ')]
P7 = ''.join('1' if (j - i) % 7 in (1, 2, 4) else '0' for i in range(7) for j in range(i + 1, 7))
H8 = '1' * 7 + ''.join('1' if (j - i) % 7 in (1, 2, 4) else '0' for i in range(7) for j in range(i + 1, 7))
out = sh("./u9_n7 cls7.txt 10 %s < /dev/null" % P7)
check('K=10 found=0' in out, "method check n=7 (host P_7): no 10-arc shaving")
out = sh("./u9_n7 cls7.txt 9 %s < /dev/null" % P7); F = found_sets(out, 'SHAVING')
check(len(F) == 54 and len(classes_of(7, F)) == 51, "method check n=7: 54 Aut(P_7)-orbits of 9-arc shavings = 51 classes (S15)")
out = sh("./u9_n8 cls8.txt 12 %s < /dev/null" % H8)
check('K=12 found=0' in out, "method check n=8 (host P_7 + source): no 12-arc shaving (S15: u(8) = 11)")
out = sh("./u9_n8 cls8.txt 11 %s < /dev/null" % H8); F = found_sets(out, 'SHAVING')
check(len(F) == 2290 and len(classes_of(8, F)) == 1617, "method check n=8: 2290 orbits of 11-arc shavings = 1617 classes (S15)")
out = sh("SHAVE4_PARITY=1 ./u9_n8 cls8.txt 11 %s < /dev/null" % H8)
check('K=11 found=0' in out, "parity mode n=8: no 11-arc Redei graph")
out = sh("SHAVE4_PARITY=1 ./u9_n8 cls8.txt 10 %s < /dev/null" % H8); F = found_sets(out, 'REDEI')
cen = [parse_arcs(l.split('arcs ')[1].split()) for l in L8 if ' m=10 ' in l]
check(classes_of(8, F) == classes_of(8, cen) and len(classes_of(8, F)) == 2, "parity mode n=8: the 10-arc Redei graphs are the census' 2 classes")
S14 = "0>1 0>3 0>6 1>4 1>5 1>7 2>3 2>5 3>4 3>7 4>7 5>6 5>8 6>8"
out = sh("./engine embedall 9 cls9.txt %s" % S14.replace('>', '\\>'))
check('avoiders=0' in out, "the 14-arc graph %s embeds in all 191536 9-classes: u(9) >= 14" % S14)
out = sh("./engine embedall 9 cls9.txt %s 0\\>2" % S14.replace('>', '\\>'))
check('avoiders=217' in out, "   (adding 0>2 makes it avoidable: 217 classes avoid it)")
tab = sh("gentourng -q -z 9 | countg --a 2>&1")
gs = {int(m.split('groupsize=')[1]): int(m.split()[0]) for m in tab.split('\n') if 'groupsize=' in m}
check(gs == {1: 188337, 3: 3047, 5: 40, 7: 4, 9: 94, 15: 8, 21: 4, 27: 1, 81: 1},
      "|Aut| census of the 191536 classes (nauty countg): 14 classes have |Aut| > 9 (15^8, 21^4, 27, 81)")
big = sh("gentourng -q -z 9 | pickg -q -a10:").split()
def d6_to_up(s):
    s = s[1:]; n = ord(s[0]) - 63; bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    return ''.join('1' if bits[i * n + j] else '0' for i in range(n) for j in range(i + 1, n))
def op6(s):
    s0 = s[1:]; n = ord(s0[0]) - 63; bits = []
    for ch in s0[1:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    return d6_of(n, [(i, j) for i in range(n) for j in range(n) if bits[j * n + i]])
can = sh("labelg -q", inp='\n'.join(big + [op6(s) for s in big]) + '\n').split()
sc = sum(1 for i in range(len(big)) if can[i] == can[i + len(big)])
check(sc == 4, "excess law (HYP-3817/3819): #{self-converse 9-classes with |Aut| > 9} = 4, predicting kappa(9) = ceil(log2 191536) + 4 = 22, u(9) = 14")
open(os.path.join(BUILD, 'pool0.txt'), 'w').write('\n'.join(d6_to_up(s) for s in big) + '\n')
C3C3 = [ (i, j) for i in range(9) for j in range(9) if i != j and ((j // 3) == (i // 3 + 1) % 3 or (i // 3 == j // 3 and j % 3 == (i % 3 + 1) % 3)) ]
tt5 = any(all((a, b) in set(C3C3) for a, b in itertools.combinations(X, 2)) for X in itertools.permutations(range(9), 5))
check(not tt5, "the host C3[C3] contains no transitive 5-subtournament")
if '--skip-u9' in sys.argv:
    print("   (skipped by --skip-u9) recorded: u9_n9 cls9.txt 15 < pool0.txt -> K=15 found=0 (see the .out file)")
else:
    t1 = time.time()
    out = sh("./u9_n9 cls9.txt 15 - < pool0.txt")
    print("   " + "\n   ".join(l for l in out.split('\n') if l.startswith('K=') or l.startswith('depth') or l.startswith('N=')))
    final_pool = [int(x) for x in [l for l in out.split('\n') if l.startswith('POOL')][0].split()[1:]]
    open(os.path.join(BUILD, 'poolfinal.txt'), 'w').write('\n'.join(CLS[9][i] for i in final_pool) + '\n')
    check('K=15 found=0' in out, "exhaustive orderly search (host C3[C3], |Aut| = 81): no 15-arc shaving on 9 vertices (%.0f s): u(9) = 14, kappa(9) = 22" % (time.time() - t1))

U9_RECORD = os.path.join(ROOT, '05-knowledge', 'results', 'procgen_shave4_20261001_u9search.out')
mx = [[tuple(map(int, a.split('>'))) for a in l.split()[1:]] for l in open(U9_RECORD).read().split('\n') if l.startswith('MAX14 ')]
bad = 0
for a in mx:
    out = sh("./engine embedall 9 cls9.txt %s" % ' '.join('%d\\>%d' % e for e in a))
    if 'avoiders=0' not in out: bad += 1
cm = set(canon_list(9, mx))
check(len(mx) == 54 and bad == 0 and len(cm) == 54 and canon_list(9, [parse_arcs(S14.split())])[0] in cm,
      "the 54 recorded maximum (14-arc) shavings on 9 vertices embed in all 191536 classes and are pairwise non-isomorphic; S_14 is one of them")
if '--with-k14' in sys.argv:
    t1 = time.time()
    out = sh("./u9_n9 cls9.txt 14 - < %s" % ('poolfinal.txt' if os.path.exists(os.path.join(BUILD, 'poolfinal.txt')) else 'pool0.txt'))
    F = found_sets(out, 'SHAVING')
    check(set(canon_list(9, F)) == cm, "re-run: the 14-arc shavings form exactly these 54 classes (%d orbits, %.0f s)" % (len(F), time.time() - t1))

# ---------------------------------------------------------------- 11. r(9)
section("11. r(9): the densest Redei graphs on 9 vertices")
R9_RECORD = os.path.join(ROOT, '05-knowledge', 'results', 'procgen_shave4_20261001_r9search.out')
rec = open(R9_RECORD).read()
ex = [[tuple(map(int, a.split('>'))) for a in l.split('arcs:')[1].split()] for l in rec.split('\n') if l.startswith('EXIST ')]
r9 = max(len(a) for a in ex)
lines_in = '\n'.join('9 %d 1 %s' % (len(a), ' '.join('%d %d' % x for x in a)) for a in ex) + '\n'
out = sh("./engine redei 9 cls9.txt", inp=lines_in)
check(out.count('REDEI n=9') == len(ex) and r9 == 11,
      "existence: the %d recorded %d-arc graphs are Redei (odd embedding count in all 191536 classes; engine backtracking)" % (len(ex), r9))
byd = [l for l in rec.split('\n') if l.startswith('REDEI_BY_DEPTH')][0]
kline = [l for l in rec.split('\n') if l.startswith('K=15')][0]
print("   recorded search (host C3[C3], parity mode, Redei tests at sizes 12..15): " + byd + " | " + kline)
check(all(int(t.split(':')[1]) == 0 for t in byd.split()[1:]) and 'found=0' in kline,
      "non-existence: no Redei graph with 12..15 arcs on 9 vertices, so r(9) = 11 (u(9) = 14)")
if '--with-r9' in sys.argv:
    t1 = time.time()
    pf = 'poolfinal.txt' if os.path.exists(os.path.join(BUILD, 'poolfinal.txt')) else 'pool0.txt'
    out = sh("SHAVE4_PARITY=1 SHAVE4_PARITY_MIN=12 ./u9_n9 cls9.txt 15 - < %s" % pf)
    byd2 = [l for l in out.split('\n') if l.startswith('REDEI_BY_DEPTH')][0]
    check(byd2 == byd and 'K=15 found=0' in out, "re-run of the parity search reproduces %s (%.0f s)" % (byd2, time.time() - t1))

print("\n%d checks, %.0f s" % (NCHECK, time.time() - T0))
print("ALL CHECKS PASSED")
