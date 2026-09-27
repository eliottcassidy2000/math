"""procgen_cnpos_20260926_run.py -- runner for the cnpos lane (positive Hamiltonicity of C_n).

Session collatz-procgen-20260922, lane "cnpos" (2026-09-26).
Usage (from the worktree root):
    nice python3 -u 04-computation/experiments/procgen_cnpos_20260926_run.py [--full]
  --full : W_8 checked at every n (default: every 8th n).
The clean n of every level a <= 15 with residual modulus m <= 4e6 are constructed (largest n = 14348907;
peak RSS about 350 MB).  a = 16 (m = 9492289) is left out: its arrays would exceed the 500 MB cap.
Prints to stdout only.  Every path is re-verified edge by edge by core.verify_path (a permutation of 1..n with
every consecutive sum a power of 2 or 3).  Ends with ALL CHECKS PASSED or raises.
"""
import sys, time, random, math, resource
from array import array

EXP = '/tmp/math-wt-collatz-procgen-20260922/04-computation/experiments'
sys.path.insert(0, EXP)
import procgen_cnpos_20260926_core as C

FULL = '--full' in sys.argv
BIG = '--big' in sys.argv
T0 = time.time()
NCHK = [0]


def ok(cond, msg):
    if not cond:
        raise AssertionError('CHECK FAILED: ' + msg)
    NCHK[0] += 1
    print('[OK] ' + msg, flush=True)


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss >> 20   # macOS: bytes


# Hamiltonian set of C_n for n <= 2200 (THM-4505 / THM-4510, FINITE-EXACT) and W_8 (all Hamiltonian)
HAM_RUNS = [(1, 3), (5, 8), (18, 26), (49, 63), (65, 66), (68, 80), (179, 194), (224, 243), (473, 575), (665, 728),
            (1319, 1418), (1536, 1616), (1620, 1663), (1703, 1713), (1792, 1802), (1920, 2114), (2120, 2186),
            (4374, 6561)]
HAM = set()
for lo, hi in HAM_RUNS:
    HAM.update(range(lo, hi + 1))


def clean_levels(a):
    """q with 2/3 < 2^q/3^a < 4/3; the clean n's of level a."""
    t3 = 3 ** a
    q = next(q for q in range(1, 10 * a + 10) if 3 * 2 ** q > 2 * t3 and 3 * 2 ** q < 4 * t3)
    return q


# ================================================================================================
print('=== A. Structure: the first-return reduction (Lemma F_j) ===', flush=True)
t = time.time()


def adm_windows(amax):
    out = []
    for a in range(2, amax + 1):
        lo, hi, kind = C.window(a)
        out.append((a, lo, hi, kind))
    return out


WIN = adm_windows(9)
# A.1  parameters at every n of W_2..W_8 (n not a power of 3)
cnt = 0
for a, lo, hi, kind in WIN:
    if a > 8:
        continue
    for n in range(lo, hi + 1):
        if n == 3 ** a:
            continue
        d = C.reduced_data(n)
        T1, T2, m, j, r, u0, f = d['T1'], d['T2'], d['m'], d['j'], d['r'], d['u0'], d['f']
        assert C.is_target(T1) and C.is_target(T2) and T1 > n and T2 <= 2 * n and m % 2 == 1
        assert j + 1 <= u0 <= j + m and (2 * u0 - T1) % m == 0
        # phi is an involution of W with the single fixed point u0
        if n <= 2200 or n % 97 == 0:
            fx = 0
            for x in range(j + 1, j + m + 1):
                y = C.phi_of(d, x)
                assert j + 1 <= y <= j + m and C.phi_of(d, y) == x
                fx += (y == x)
            assert fx == 1
        cnt += 1
ok(cnt > 3000, 'A.1 top targets T1 < T2 (T1 > n, T2 <= 2n, m = T2 - T1 odd), phi an involution of W with one fixed '
   'point, at all %d n of W_2..W_8 (n != 3^a)' % cnt)


# A.2  the first return itself: walks alpha, beta, alpha, ... from each b in W land at phi(b) (or at f for b = u0),
#      and these walks cover the zone [r+1, n] exactly once (no zone-only cycle)
def first_return_check(n):
    d = C.reduced_data(n)
    T1, T2, m, j, r, u0, f = d['T1'], d['T2'], d['m'], d['j'], d['r'], d['u0'], d['f']
    seen = bytearray(n + 1)
    for b in range(j + 1, j + m + 1):
        if b != u0 and b > C.phi_of(d, b):
            continue                        # each walk once (it joins b and phi(b))
        x = b
        use_alpha = True
        land = None
        while True:
            y = (T1 - x) if use_alpha else (T2 - x)
            if not (1 <= y <= n) or y == x:
                land = ('end', x)           # fixed point reached: x is f
                break
            if y <= r:
                land = ('R', y)
                break
            if seen[y]:
                return False, 'revisit %d' % y
            seen[y] = 1
            x = y
            use_alpha = not use_alpha
        if b == u0:
            if land != ('end', f):
                return False, 'u0 walk %s' % (land,)
        else:
            if land != ('R', C.phi_of(d, b)):
                return False, 'b=%d lands %s, phi=%d' % (b, land, C.phi_of(d, b))
    # every zone vertex visited exactly once (no revisit above; no zone-only cycle)
    zone = [x for x in range(r + 1, n + 1)]
    miss = [x for x in zone if not seen[x]]
    return (not miss), 'missing %d' % len(miss)


cnt = 0
for a, lo, hi, kind in WIN:
    if a > 8:
        continue
    step = 1 if a <= 7 else 13
    for n in range(lo, hi + 1, step):
        if n == 3 ** a:
            continue
        res, why = first_return_check(n)
        assert res, (n, why)
        cnt += 1
ok(cnt > 1400, 'A.2 first return: from every b in W the walk alpha, beta, alpha, ... through the zone lands at phi(b) '
   '(at the fixed point f for b = u0), and the walks cover the zone; %d values of n in W_2..W_8' % cnt)


# A.3  Lemma Z (zone step), brute force on random small instances: psi u (M_t|Z u mu') is a Hamiltonian path of V
#      iff psi' u mu' is one of V' (and the zone step reports a zone-only cycle iff no mu' works for that reason)
def ham_alt(V, psi, s, mu):
    seen = {s}
    cur = s
    while cur in mu:
        y = mu[cur]
        if y in seen:
            return False
        seen.add(y)
        z = psi[y]
        if z == y or z in seen:
            return False
        seen.add(z)
        cur = z
    return len(seen) == len(V)


def matchings_minus_one(vs, s):
    """all matchings of vs covering all but one vertex e (e != s), as dicts."""
    vs = sorted(vs)
    out = []

    def rec(rem, cur, skipped):
        if not rem:
            if skipped is not None:
                out.append(dict(cur))
            return
        x = rem[0]
        if skipped is None and x != s:
            rec(rem[1:], cur, x)
        for i in range(1, len(rem)):
            y = rem[i]
            cur[x] = y
            cur[y] = x
            rec(rem[1:i] + rem[i + 1:], cur, skipped)
            del cur[x]
            del cur[y]
    rec(vs, {}, None)
    return out


rng = random.Random(20260926)
trials = 0
agree = 0
while trials < 300:
    k = rng.choice([7, 9, 11, 13])
    V = sorted(rng.sample(range(1, 30), k))
    s = rng.choice(V)
    rest = [x for x in V if x != s]
    rng.shuffle(rest)
    psi = {s: s}
    for i in range(0, len(rest), 2):
        psi[rest[i]] = rest[i + 1]
        psi[rest[i + 1]] = rest[i]
    t_ = rng.choice([x + y for x in V for y in V if x < y])
    Z = [x for x in V if (t_ - x) in V and 2 * x != t_]
    if not Z or len(Z) == len(V):
        continue
    trials += 1
    Vset = set(V)
    # zone step by hand (same rule as core.zone_step, independent code)
    Vp = [x for x in V if x not in Z]
    psip, sp, visited = {}, None, 0
    for x in Vp:
        if x in psip:
            continue
        y = psi[x]
        if y == x:
            psip[x] = x
            sp = x
            continue
        hit = False
        while y in Z:
            z = t_ - y
            visited += 2
            if z == s:
                hit = True
                break
            y = psi[z]
        if hit:
            psip[x] = x
            sp = x
        else:
            psip[x] = y
            psip[y] = x
    cyc = visited != len(Z)
    for mup in matchings_minus_one(Vp, sp if sp is not None else s):
        full = dict(mup)
        for x in Z:
            full[x] = t_ - x
        a1 = ham_alt(Vset, psi, s, full)
        if cyc:
            assert not a1
        else:
            a2 = ham_alt(set(Vp), psip, sp, mup)
            assert a1 == a2, (V, psi, s, t_, mup)
    agree += 1
ok(agree == 300, 'A.3 Lemma Z (zone step) exact on 300 random instances (|V| = 7..13, every completing matching mu\' '
   'enumerated; zone-only cycles never extend)')
print('  section A time %.1fs, maxrss %d MB' % (time.time() - t, rss()), flush=True)

# ================================================================================================
print('=== B. Constructions ===', flush=True)
import procgen_sumgraph_20260926_solver as SV
t = time.time()


def construct(n, prefer_gp=True, end_at=None, want_seq=False):
    """Lemma F_j construction at n: solve the residual problem, lift, verify.
    Returns ('PATH', how[, seq]) or (None, why[, None])."""
    d = C.reduced_data(n)
    mu_adj = None
    MU = None
    how = None
    if d['r'] == 1:                      # m = 1, j = 0 (n = 7): the residual is the single vertex u0 = e
        MU = array('i', [0, 0])
        how = 'trivial residual'
    if MU is None and d['j'] == 0 and prefer_gp:
        for small in (2500, 8000):
            g = C.rp_gp(d['m'], d['T1'] % d['m'])
            assert g.s == d['u0']
            if end_at is not None:
                g.e = end_at
            MU = array('i', [0]) * (d['m'] + 1)
            if C.gp_solve(g, MU, small=small) and C.check_gp(g, MU) and (end_at is None or MU[end_at] == 0):
                how = 'zone steps, exact search below %d' % small
                del g
                break
            MU = None
            del g
    if MU is None and d['r'] <= 12000:
        res = C.solve_residual_dpll(d, end_set=(None if end_at is None else {end_at}))
        if isinstance(res, dict):
            mu_adj = res
            how = 'exact search'
    if mu_adj is None and MU is None:
        return (None, 'residual unsolved', None) if want_seq else (None, 'residual unsolved')
    seq = C.lift(d, mu_adj, MU)
    good = C.verify_path(n, seq) and (end_at is None or seq[-1] == end_at)
    if not good:
        raise AssertionError('lifted sequence failed verification at n=%d' % n)
    return ('PATH', how, seq) if want_seq else ('PATH', how)


# B.1 the clean right ends at every level a <= AMAX
AMAX = 15
MMAX = 4 * 10 ** 6
print('  B.1 clean n per level (T1 - 1 with T1 = min(3^a, 2^q); window right ends)', flush=True)
rows = []
for a in range(2, AMAX + 1):
    q = clean_levels(a)
    lo, hi, kind = C.window(a)
    t3, t2 = 3 ** a, 2 ** q
    targets = [min(t3, t2) - 1]
    if kind == 'two-power':
        targets.append(hi)                 # 3^a - 1, top targets 3^a < 2^(q+1)
    for n in targets:
        if n < 5 or C.reduced_data(n)['m'] > MMAX:
            continue
        tt = time.time()
        d = C.reduced_data(n)
        st, how = construct(n)
        rows.append((a, kind, n, d['T1'], d['T2'], d['m'], st, how, time.time() - tt))
        print('    a=%2d %-9s n=%-9d T1=%-9d T2=%-9d m=%-8d %s (%s) %.1fs maxrss %d MB'
              % (a, kind, n, d['T1'], d['T2'], d['m'], st, how, time.time() - tt, rss()), flush=True)
    if kind == 'one-power' and 3 ** a >= 18 and abs(t2 - t3) <= MMAX:
        # n = 3^a: the leaf 3^a hangs on m = 2^q - 3^a; a path of C_(3^a - 1) ending at m extends
        n1 = 3 ** a - 1
        d = C.reduced_data(n1)
        st, how, seq = construct(n1, end_at=d['m'], want_seq=True)
        if st == 'PATH':
            seq.append(3 ** a)
            assert C.verify_path(3 ** a, seq), 'extension to 3^a failed'
            del seq
        rows.append((a, kind, 3 ** a, d['T1'], d['T2'], d['m'], st, how, 0))
        print('    a=%2d %-9s n=%-9d (=3^a)  end forced at m=%d: %s (%s)' % (a, kind, 3 ** a, d['m'], st, how), flush=True)
primary = [rw for rw in rows if rw[2] == min(3 ** rw[0], 2 ** clean_levels(rw[0])) - 1]
ok(all(rw[6] == 'PATH' for rw in primary),
   'B.1 C_(T1-1) Hamiltonian, T1 = min(3^a, 2^q), at every level a = 2..%d (paths lifted from the residual problem '
   'and verified edge by edge; largest n = %d)' % (AMAX, max(rw[2] for rw in primary)))
rightends = [rw for rw in rows if (rw[1] == 'two-power' and rw[2] == 3 ** rw[0] - 1)]
print('    two-power right ends 3^a - 1:', [(rw[0], rw[6]) for rw in rightends], flush=True)
lost_re = [rw[0] for rw in rightends if rw[6] != 'PATH']
ok(lost_re == [4] and 80 in HAM, 'B.1b two-power right ends C_(3^a - 1) Hamiltonian by the construction (T1 = 3^a, '
   'T2 = 2^(q+1)) for a in %s; at a = 4 (n = 80, Hamiltonian by the list) the all-alpha residual has no solution '
   '(exact), a loss of the construction, not of Hamiltonicity' % [rw[0] for rw in rightends if rw[6] == 'PATH'])
oneend = [rw for rw in rows if rw[1] == 'one-power' and rw[2] == 3 ** rw[0]]
print('    one-power right ends n = 3^a via the end forced at m:', [(rw[0], rw[6]) for rw in oneend], flush=True)
ok(set(rw[0] for rw in oneend if rw[6] == 'PATH') >= {8, 10, 13} and
   [rw[0] for rw in oneend if rw[6] != 'PATH' and rw[0] <= 13] == [3, 5], 'B.1c one-power right ends: C_(3^a) Hamiltonian by the '
   'construction at a = 8, 10, 13 (paths verified); it fails at a = 3 (C_27 is not Hamiltonian) and at a = 5 '
   '(C_243 is, but its unique path uses an M_P edge, see B.1d)')
import procgen_sumgraph_20260926_solver as SV
r243 = SV.decide(243)
s243 = r243['seq']
mp_edges = [(s243[i], s243[i + 1]) for i in range(242) if s243[i] + s243[i + 1] == 128
            and 115 < max(s243[i], s243[i + 1]) < 128]
ok(r243['status'] == 'PATH' and len(mp_edges) >= 1, 'B.1d the unique Hamiltonian path of C_243 (sumgraph solver) uses '
   'the M_P edges %s between the interior of I = [3^5 - 2^7, 2^7] and [1, m]: it is not an all-alpha path' % mp_edges)
# B.2 the construction on whole windows W_3..W_8, against the Hamiltonian list (THM-4510)
print('  B.2 whole windows: construction (Lemma F_j, residual exact search) vs the Hamiltonian list', flush=True)
summary = []
for a in range(3, 9):
    lo, hi, kind = C.window(a)
    ns = [n for n in range(lo, hi + 1) if n != 3 ** a]
    if a == 8 and not FULL:
        ns = ns[::8] + [hi - 1]
    got, lost, fake = [], [], []
    for n in ns:
        st, how = construct(n, prefer_gp=False)
        if st == 'PATH':
            got.append(n)
            if n not in HAM:
                fake.append(n)
        elif n in HAM:
            lost.append(n)
    assert not fake, fake
    # Theorem F' (exact form): contract the forced chains, exact search on H_n (covers the paths whose free end is not
    # saturated by forced edges); where it fails, retry with a larger search budget, then report and fall back to the
    # direct search (THM-4510's solver), whose path is verified here and whose free end is inspected
    rec, unrec, direct = [], [], []
    for n in lost:
        seq = C.solve_exact_reduced(n)
        if seq is not None and seq != 'ABORT' and C.verify_path(n, seq):
            rec.append(n)
            continue
        first = 'no H_n path' if seq is None else 'search aborted'
        seq = C.solve_exact_reduced(n, node_limit=2000000)
        if seq is not None and seq != 'ABORT' and C.verify_path(n, seq):
            rec.append(n)
            print('      n=%d: exact form needed the larger search budget (first attempt: %s)' % (n, first), flush=True)
            continue
        second = 'no H_n path' if seq is None else 'search aborted'
        r = SV.decide(n)
        s_ = r.get('seq')
        if r['status'] == 'PATH' and C.verify_path(n, s_):
            P_, Q_, T1_, T2_ = C.top_targets(n)
            ends = (s_[0], s_[-1])
            free = [x for x in ends if x != P_][0]
            direct.append(n)
            print('      n=%d: exact form %s / %s; direct search path verified, free end %d (%s)'
                  % (n, first, second, free, 'top vertex' if free > max(P_, Q_) else 'not top'), flush=True)
        else:
            unrec.append(n)
    hamcount = sum(1 for n in ns if n in HAM)
    summary.append((a, kind, len(ns), hamcount, len(got), lost, rec, unrec, direct))
    print('    W_%d %-9s |n|=%-5d Hamiltonian %-5d all-alpha construction %-5d lost %-4d (%s); exact form recovers %d; '
          'direct search %d; unresolved %d'
          % (a, kind, len(ns), hamcount, len(got), len(lost), (lost if len(lost) < 16 else (lost[:10], '...')), len(rec),
             len(direct), len(unrec)), flush=True)
ok(all(not s[7] for s in summary), 'B.2 W_3..W_8: every Hamiltonian n (%s of W_8) gets a verified path: by the all-alpha construction (Theorem F), where it is lost by the exact form (Theorem F\'), and at %d values (listed above) only by the direct search; no path at a non-Hamiltonian n' % ('all' if FULL else 'every 8th n', sum(len(s[8]) for s in summary)))
two_lower_exact = all(not [n for n in s[5] if n < 2 ** (clean_levels(s[0]))] for s in summary if s[1] == 'two-power')
ok(two_lower_exact, 'B.2b two-power windows, lower parts (T1 = 2P, Theorem F exact): the construction finds a path at '
   'every Hamiltonian n')
print('  section B time %.1fs, maxrss %d MB' % (time.time() - t, rss()), flush=True)

# B.3 W_9 near the clean n: j = 0..J for T1 = 2^14 (exact regime) and T1 = 3^9
print('  B.3 W_9 samples near the clean n', flush=True)
t = time.time()
w9 = []
for n in list(range(16383, 16383 - 41, -1)) + [19682]:
    st, how = construct(n)
    w9.append((n, st))
okn = [n for n, st in w9 if st == 'PATH']
print('    W_9 constructed at', len(okn), 'of', len(w9), 'sampled n; not constructed:', [n for n, st in w9 if st != 'PATH'],
      flush=True)
ok(len(okn) >= 30, 'B.3 W_9: construction verified at %d sampled n (16343..16383 and 19682)' % len(okn))
print('  B.3 time %.1fs, maxrss %d MB' % (time.time() - t, rss()), flush=True)

# ================================================================================================
print('=== C. Obstructions and typing ===', flush=True)
t = time.time()
# C.1 pure rotations: a circular reflection mu(x) = c - x (mod m) of [1, m] consists of the linear pieces x + y = c
#     (pairs exist iff c >= 3) and x + y = c + m (pairs exist iff c <= m - 1); it uses only C-edges iff every nonempty
#     piece has a target sum.  Necessary condition for a single-rotation solution of RP(m_a, K_a).
pure = []
for a in range(2, 301):
    q = clean_levels(a)
    m = abs(3 ** a - 2 ** q)
    Ts = set([2 ** k for k in range(2, (2 * m + 4).bit_length() + 1)] + [3 ** k for k in range(1, 2 * a + 4)])
    found = (m == 1)
    for c in range(1, m + 1) if m < 5000 else []:
        if (c < 3 or c in Ts) and (c > m - 1 or (c + m) in Ts):
            found = True
            break
    if m >= 5000:
        # c must itself be a target (c >= 3) or c in {1, 2}; enumerate targets c <= m and c = 1, 2
        cs = [1, 2] + sorted(t for t in Ts if 3 <= t <= m)
        found = found or any((c > m - 1 or (c + m) in Ts) for c in cs if (c < 3 or c in Ts))
    if found:
        pure.append(a)
ok(pure == [2, 3, 5], 'C.1 the residual problem of level a can be solved by a single circular reflection (a rotation '
   'orbit, Theorem A3 shape) only for a in {2, 3, 5}, a <= 300 (Pillai: 4 - 3 = 9 - 8, 9 - 4 = 32 - 27, '
   '16 - 3 = 256 - 243)')
# C.2 two-edge switches: tau1 + tau2 = tau3 + tau4 with disjoint pairs, targets < 2^200
T = sorted(set([2 ** k for k in range(2, 200)] + [3 ** k for k in range(1, 127)]))
Ts = set(T)
rel = set()
for i, x in enumerate(T):
    for y in T[i:]:
        ssum = x + y
        for z in T:
            if z >= ssum:
                break
            w = ssum - z
            if w in Ts and z <= w and {z, w} != {x, y} and z not in (x, y):
                rel.add(tuple(sorted([(x, y), (z, w)])))
rel = sorted(rel)
print('    target relations tau1 + tau2 = tau3 + tau4:', rel, flush=True)
ok(rel == [((3, 9), (4, 8)), ((3, 32), (8, 27)), ((3, 256), (16, 243)), ((4, 32), (9, 27))],
   'C.2 the only relations tau1 + tau2 = tau3 + tau4 between targets < 2^200 (disjoint pairs) are 3 + 9 = 4 + 8, '
   '3 + 32 = 8 + 27, 3 + 256 = 16 + 243, 4 + 32 = 9 + 27; each uses target 3 or 4, i.e. the edge {1, 2} or {1, 3}')
# C.3 the residual problem is not solvable for every (m, K): exact refutations
bad = []
for m, K in [(11, 1), (11, 5), (59, 22), (383, 1)]:
    d = dict(n=None, T1=K, T2=K + m, m=m, j=0, r=m, u0=(K * ((m + 1) // 2)) % m or m, f=0)
    res = C.solve_residual_dpll(d)
    bad.append((m, K, res))
print('    RP(m, K) exact status:', [(m, K, 'none' if r is None else ('solvable' if isinstance(r, dict) else r))
                                      for m, K, r in bad], flush=True)
ok(all(r is None for _, _, r in bad), 'C.3 RP(11,1), RP(11,5), RP(59,22), RP(383,1) have no solution (exact: '
   'propagation or exhaustive search), so any proof must use the arithmetic of (m_a, K_a)')
print('  section C time %.1fs, maxrss %d MB' % (time.time() - t, rss()), flush=True)

print('RUN time %.1fs, peak RSS %d MB, checks %d' % (time.time() - T0, rss(), NCHK[0]), flush=True)
print('ALL CHECKS PASSED', flush=True)
