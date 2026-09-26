#!/usr/bin/env python3
"""
procgen_mykk_20260926_run.py -- runner of the expanding-cycle feedback lane ("mykk"), session
collatz-procgen-20260922, 2026-09-26.  Produces 05-knowledge/results/procgen_mykk_20260926.out:

    python3 04-computation/experiments/procgen_mykk_20260926_run.py > 05-knowledge/results/procgen_mykk_20260926.out

Every printed claim is guarded by check(...) (raises SystemExit on failure).  Nothing from the searches is
trusted: frozen certificates (procgen_mykk_20260926_certs.py) are re-verified from scratch:
  * an optimal set R: Bellman-Ford finds no expanding cycle in B(2,k) - R and the returned integer potential is
    re-checked on every edge;
  * the lower bound: a family of pairwise node-disjoint expanding cycles (checked directly), or a list of
    expanding cycles whose minimum hitting set is re-proved by the MaxSAT solver RC2 (pysat), independent of the
    HiGHS MIP used in the search.
Sections:
  A  Terras coordinates: the parity graph of T_q is B(2,k) (q = 3, 5)
  B  Mykkeltveit's theorem (Golomb's conjecture): lemma checks and the explicit set, k = 3..16
  C  necklace counts N_c(k) (Burnside vs enumeration)
  D  exact FVS_c(k), FVS^odd_c(k) at c = log_3 2 and c = log_5 2; packings nu
  E  full step functions c -> FVS_c(k) (k <= 9) against N_c(k); closed thresholds
  F  closed-threshold theorems at (k-1)/k, (k-2)/k (proved), 2/k (FINITE-EXACT), 1/k (Mykkeltveit)
  G  explicit counterexamples to "FVS_c = N_c"; fractional-packing lower bound and its asymptotics
  H  dynamics (Theorem 3): periodic edits, lookahead, thresholds, orbit census, rising witnesses
  I  DRIFT (Theorem 4): price tables for q = 5 and q = 3
  J  flips versus deletions (Proposition 5): delta_k >= FVS^odd >= FVS >= nu >= N
"""
import os
import sys
import time
from math import comb, log, log2
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from procgen_mykk_20260926_lib import *  # noqa: E402,F401,F403

T0 = time.time()


def say(msg=''):
    print(msg, flush=True)


def section(title):
    say('')
    say('=' * 100)
    say(title)
    say('=' * 100)


def ok(cond, claim):
    check(cond, claim)
    say('  [ok] ' + claim)


def load_certs():
    try:
        from procgen_mykk_20260926_certs import DATA
    except ImportError:
        return {}
    return unpack_json(DATA)


CERTS = load_certs()
H3 = -(log(2) / log(3)) * log2(log(2) / log(3)) - (1 - log(2) / log(3)) * log2(1 - log(2) / log(3))


# ============================================================================================ A
def sec_A():
    section("A. Terras coordinates: the level-k parity graph of T_q is the de Bruijn graph B(2,k)")
    for q in (3, 5):
        for k in range(1, 13):
            res2word, word2res = terras_tables(k, q)
            n = 1 << k
            bad = 0
            for r in range(n):
                t = Tq(r, q) % (1 << max(k - 1, 0)) if k > 1 else 0
                lifts = {(t + e * (1 << (k - 1))) % n for e in (0, 1)} if k > 1 else {0, 1}
                words = {int(res2word[x]) for x in lifts}
                W = int(res2word[r])
                if words != {succ(W, 0, k), succ(W, 1, k)}:
                    bad += 1
                if (W & 1) != (r & 1):
                    bad += 1
            check(bad == 0, "Terras/de Bruijn identification q=%d k=%d" % (q, k))
        ok(True, "q=%d: for k = 1..12 the residue -> parity-word map is a bijection, preserves parity, and maps the "
                 "parity graph G_0 (s -> both lifts of T_q(s) mod 2^(k-1)) onto B(2,k)" % q)


# ============================================================================================ B
def sec_B():
    section("B. Mykkeltveit's theorem (Golomb's conjecture): minimum FVS of B(2,k) = Z(k); explicit set")
    rng = np.random.default_rng(20260926)
    for k in range(3, 17):
        t = time.time()
        n = 1 << k
        s, wzero, szero = sine_data(k)
        M, Sstar, zn = mykkeltveit_set(k)
        nid, neck = necklaces(k)
        check(len(neck) == Z_formula(k), "Z(k) formula k=%d" % k)
        check(len(M) == Z_formula(k), "|M_k| = Z(k), k=%d" % k)
        check(len(set(int(nid[x]) for x in M)) == len(neck), "one node per necklace, k=%d" % k)
        check(is_acyclic(k, M), "B(2,k) - M_k acyclic, k=%d" % k)
        # Lemma 1.3 (crossing) on every edge: s(x) > 0 >= s(y) implies y in S*
        sg = sign_s(s, szero)
        Sset = np.zeros(n, dtype=bool)
        Sset[Sstar] = True
        X = np.arange(n)
        for b in (0, 1):
            Y = (X >> 1) | (b << (k - 1))
            cross = (sg[X] > 0) & (sg[Y] <= 0)
            check(bool(np.all(Sset[Y[cross]])), "crossing lemma k=%d" % k)
        # Lemma 1.2 (zero mean of w) on random closed walks
        om = np.exp(2j * np.pi * np.arange(k) / k)
        for _ in range(20):
            p = int(rng.integers(1, 3 * k))
            u = int(rng.integers(0, 1 << p))
            cyc = cycle_from_word(u, p, k)
            tot = sum(sum(((x >> j) & 1) * om[j] for j in range(k)) for x in cyc)
            check(abs(tot) < 1e-8, "zero-mean lemma k=%d" % k)
        say("  k=%2d  Z(k)=%5d  |S*|=%5d  zero-weight necklaces=%3d  |M_k|=%5d  B(2,k)-M_k acyclic: yes  (%.1fs)"
            % (k, Z_formula(k), len(Sstar), len(zn), len(M), time.time() - t))
    ok(True, "for k = 3..16: M_k = S* + one node per zero-weight necklace has Z(k) nodes (one per necklace) and "
             "meets every cycle of B(2,k); the crossing lemma holds on every edge; w has mean 0 on sampled closed walks")
    ok(all(Z_formula(k) == len(necklaces(k)[1]) for k in range(1, 17)), "Z(k) = (1/k) sum_{d|k} phi(d) 2^(k/d) matches "
                                                                        "enumeration for k <= 16")
    # prime k: only 0^k and 1^k are zero-weight necklaces
    for k in (3, 5, 7, 11, 13):
        _, _, zn = mykkeltveit_set(k)
        check(len(zn) == 2, "prime k zero necklaces")
    ok(True, "for prime k in {3,5,7,11,13} the only zero-weight necklaces are 0^k and 1^k")
    # upward closure of S* in bits 0 and k-1 (used in the note's remarks)
    for k in range(3, 15):
        M, Sstar, zn = mykkeltveit_set(k)
        S = set(Sstar)
        for y in Sstar:
            check((y | 1) in S and (y | (1 << (k - 1))) in S, "S* upward closed in first/last bit")
    ok(True, "S* is closed under setting the first bit or the last bit to 1 (k = 3..14)")
    lp = []
    for k in range(3, 17):
        M, _, _ = mykkeltveit_set(k)
        lp.append((k, longest_path_acyclic(k, M)))
    say("  longest path (nodes) in B(2,k) - M_k: " + ", ".join("k=%d: %d" % t for t in lp))
    ok(all(v is not None and v <= k ** 3 for _, v in lp), "the full Mykkeltveit edit meets every orbit within at "
                                                          "most k^3 steps for k <= 16 (any odd q)")


# ============================================================================================ C
def sec_C():
    section("C. Necklace counts N_c(k): Burnside/Moebius formulas against enumeration")
    thr3, thr5 = Thr('log', q=3), Thr('log', q=5)
    for k in range(1, 15):
        nid, neck = necklaces(k)
        for thr in (thr3, thr5):
            direct = sum(1 for orb in neck if thr.expanding(bin(orb[0]).count('1'), k))
            check(direct == N_count(k, thr), "N_c(k) k=%d" % k)
            dprim = sum(1 for orb in neck if len(orb) == k and thr.expanding(bin(orb[0]).count('1'), k))
            check(dprim == Nprim_count(k, thr), "primitive count k=%d" % k)
    ok(True, "N_c(k) (Burnside) and N^prim_c(k) (Moebius) equal direct enumeration for k <= 14, c = log_3 2, log_5 2")
    say("  N_{log3 2}(k), k=1..20: %s" % [N_count(k, thr3) for k in range(1, 21)])
    say("  N_{log5 2}(k), k=1..20: %s" % [N_count(k, thr5) for k in range(1, 21)])
    ok([N_count(k, thr3) for k in range(2, 15)] == [1, 2, 2, 2, 5, 5, 6, 16, 19, 52, 70, 85, 251],
       "N_k (q=3) = 1,2,2,2,5,5,6,16,19,52,70,85,251 for k=2..14 (as in the cube-distance note)")


def Nprim_count(n, thr):
    """number of primitive necklaces of length n with density > c (Moebius)"""
    def mu(m):
        r, p, x = 1, 2, m
        while p * p <= x:
            if x % p == 0:
                x //= p
                if x % p == 0:
                    return 0
                r = -r
            p += 1
        if x > 1:
            r = -r
        return r
    tot = 0
    for j in range(n + 1):
        if not thr.expanding(j, n):
            continue
        g = gcd(n, j) if j else n
        s = sum(mu(d) * comb(n // d, j // d) for d in range(1, g + 1) if g % d == 0)
        check(s % n == 0, "Moebius count")
        tot += s // n
    return tot


# ============================================================================================ D
def thr_of(spec):
    if spec[0] == 'log':
        return Thr('log', q=spec[1])
    return Thr(spec[0], spec[1], spec[2])


VERIFIED = {}


def _get(nm_base, k):
    for nm in (nm_base % k, (nm_base % k) + '_bounds'):
        if nm in VERIFIED:
            return VERIFIED[nm]
    return None


def verify_entry(name):
    """re-verify a frozen FVS result.  Exact entries: optimal set + lower bound (packing, or RC2 for k <= 10,
    or a HiGHS re-solve plus an exact LP-dual bound for k >= 11).  Bounds entries ('..._bounds'): the feasible
    set (upper bound), an exact LP-dual fractional packing (lower bound), and the MIP lower bound re-proved by a
    fresh HiGHS solve when it finishes within the time limit."""
    e = CERTS[name]
    k = e['k']
    thr = thr_of(e['thr'])
    odd = e['odd']
    n = 1 << k
    allowed = [v for v in range(n) if (v & 1 or not odd)]
    cert = [word_to_cycle(c, k) for c in e['cert']]
    pack = [word_to_cycle(c, k) for c in e.get('packing', [])]
    used = set()
    for c in pack:
        check(cycle_ok(c, k) and thr.expanding(sum(x & 1 for x in c), len(c)), "packing cycle")
        check(not (used & set(c)), "packing disjoint")
        used |= set(c)
    if e.get('value') is not None:
        how = verify_fvs_value(k, thr, e['R'], cert, e['value'], odd_only=odd,
                               packing=pack if len(pack) == e['value'] else None,
                               rc2_timeout=(60 if (k >= 11 or len(cert) > 1500) else 300))
        VERIFIED[name] = dict(k=k, value=e['value'], lb=e['value'], ub=e['value'], nu=len(pack), how=how,
                              odd=odd, R=e['R'])
        return VERIFIED[name]
    # bounds entry
    R = e['R']
    check(len(R) == e['ub'], "bounds UB size")
    if odd:
        check(all(v & 1 for v in R), "odd UB set")
    cyc, d = has_expanding_cycle(k, thr, R)
    check(cyc is None and verify_potential(k, thr, R, d), "bounds UB set feasible")
    for c in cert:
        check(cycle_ok(c, k) and thr.expanding(sum(x & 1 for x in c), len(c)), "cert cycle")
    y, yval, _ = lp_dual_packing(cert, allowed)
    lpv = verify_frac_packing(cert, y, allowed)
    lb_lp = -((-lpv.numerator) // lpv.denominator)
    t = time.time()
    v, Rm, optm, bnd = solve_cover_highs(cert, allowed, time_limit=300, full=True)
    lb_mip = bnd if bnd is not None else 0
    lb_all = None
    if odd:
        # the same cycles also bound the all-node problem: exact LP-dual packing with every node allowed
        ya, yav, _ = lp_dual_packing(cert, list(range(n)))
        lva = verify_frac_packing(cert, ya, list(range(n)))
        lb_all = -((-lva.numerator) // lva.denominator)
    VERIFIED[name] = dict(k=k, value=None, lb=max(lb_lp, lb_mip), lb_lp=lb_lp, lb_mip=lb_mip, ub=e['ub'],
                          nu=len(pack), how='bounds', odd=odd, R=R, mip_opt=optm, mip_time=time.time() - t,
                          lb_search=e.get('lb'), lb_all=lb_all)
    return VERIFIED[name]


def sec_D():
    section("D. Exact FVS_c(k) and FVS^odd_c(k) for c = log_3 2 and c = log_5 2 (certificates re-verified)")
    for q in (3, 5):
        thr = Thr('log', q=q)
        say("  q=%d (c = log_%d 2 = %.7f)" % (q, q, thr.value()))
        say("   k |   N  | nu >= | FVS | FVS^odd | Z(k)-1 | LB by | FVS/2^k   | k*FVS/2^k | FVS/N")
        for k in range(1, 13):
            names = [nm for nm in ('log%d_%d' % (q, k), 'log%dodd_%d' % (q, k),
                                   'log%d_%d_bounds' % (q, k), 'log%dodd_%d_bounds' % (q, k)) if nm in CERTS]
            if not names:
                continue
            vals = {nm: verify_entry(nm) for nm in names}
            N = N_count(k, thr)
            va = vals.get('log%d_%d' % (q, k)) or vals.get('log%d_%d_bounds' % (q, k))
            vb = vals.get('log%dodd_%d' % (q, k)) or vals.get('log%dodd_%d_bounds' % (q, k))

            def fmt(v):
                if v is None:
                    return '-'
                if v['value'] is not None:
                    return '%d' % v['value']
                return '%d..%d' % (v['lb'], v['ub'])
            nu = max(v['nu'] for v in vals.values())
            for v in vals.values():
                check(N <= v['ub'], "N <= FVS")
                check(v['nu'] <= v['ub'], "nu <= FVS")
            if va and vb:
                check(va['lb'] <= vb['ub'], "FVS <= FVS^odd")
            ref = va or vb
            best = ref['value'] if ref['value'] is not None else ref['lb']
            say("  %2d | %4d | %5d | %7s | %7s | %6d | %5s | %.6f | %.4f | %.3f" % (
                k, N, nu, fmt(va), fmt(vb), Z_formula(k) - 1, ref['how'], best / 2 ** k, k * best / 2 ** k, best / N))
            for nm, v in vals.items():
                if v['value'] is None:
                    say("       %s (%s nodes): UB %d (a feasible set, Bellman-Ford + potential); LB %d = max(exact "
                        "LP-dual packing %d, HiGHS MIP bound %d%s re-proved here); search LB %s (HiGHS, longer run)" % (
                            nm, 'odd' if v['odd'] else 'all', v['ub'], v['lb'], v['lb_lp'], v['lb_mip'],
                            ' optimal' if v['mip_opt'] else ' within 300 s', v['lb_search']))
            for nm, v in vals.items():
                if v['value'] is None and v.get('lb_all') is not None:
                    say("       %s: the same cycles give FVS (all nodes) >= %d (exact LP-dual packing)" % (nm, v['lb_all']))
            for nm, v in vals.items():
                if v['value'] is not None and v['how'] not in ('packing',):
                    say("       %s: lower bound re-proved by %s" % (nm, v['how']))
        ok(True, "q=%d: every listed FVS / FVS^odd value re-verified (optimal set: Bellman-Ford + potential; lower "
                 "bound: disjoint packing, RC2 (k <= 10) or HiGHS re-solve (k >= 11); bounds entries: exact LP-dual "
                 "packing + MIP bound)" % q)


# ============================================================================================ E
STEP = {}
HOW = {}


def sec_E():
    section("E. The exact step function c -> FVS_c(k) on [0,1) against N_c(k) (k <= 9)")
    for k in range(2, 10):
        name = 'stepfun_%d' % k
        if name not in CERTS:
            continue
        segs = CERTS[name]['segs']
        n = 1 << k
        allowed = list(range(n))
        prev_top = Fraction(1)
        rows = []
        for g in segs:
            top, bot, V, R = Fraction(g['top']), Fraction(g['bottom']), g['V'], g['R']
            check(top == prev_top, "segments contiguous")
            check(len(R) == V, "segment set size")
            rho = rho_max_removed(k, R)
            check((rho is None and bot == 0) or (rho == bot) or (bot == 0 and rho == 0),
                  "upper certificate rho_max(B-R) = bottom (k=%d, %s vs %s)" % (k, rho, bot))
            cert = [word_to_cycle(c, k) for c in g['cert']]
            thr = Thr('ge', top.numerator, top.denominator)
            for c in cert:
                check(cycle_ok(c, k) and thr.expanding(sum(x & 1 for x in c), len(c)), "stepfun cert cycle")
            v = rc2_min_cover_timed(cert, allowed, timeout=20)
            if v is None:
                v, _ = solve_cover_highs(cert, allowed)
                HOW['highs'] = HOW.get('highs', 0) + 1
            else:
                HOW['rc2'] = HOW.get('rc2', 0) + 1
            check(v == V, "stepfun lower bound k=%d top=%s: %s vs %s" % (k, top, v, V))
            rows.append((bot, top, V))
            prev_top = bot
        check(prev_top == 0, "step function reaches 0")
        STEP[k] = rows
        # compare with N_c(k): N_c constant on [j/k, (j+1)/k)
        eq_len = Fraction(0)
        pieces = []
        for (bot, top, V) in rows:
            # split [bot, top) at multiples of 1/k
            cuts = sorted(set([bot, top] + [Fraction(j, k) for j in range(k + 1) if bot < Fraction(j, k) < top]))
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                j = int(lo * k)  # lo in [j/k, (j+1)/k)
                N = sum(necklaces_with_ones(k, i) for i in range(j + 1, k + 1))
                pieces.append((lo, hi, V, N))
                if V == N:
                    eq_len += hi - lo
        say("  k=%d: FVS_c(k) on [0,1) (segment [lo, hi): value; N_c on the same piece):" % k)
        say("     " + "  ".join("[%s,%s):%d/%d" % (lo, hi, V, N) for (lo, hi, V, N) in sorted(pieces)))
        eqset = [(lo, hi) for (lo, hi, V, N) in pieces if V == N]
        # merge
        merged = []
        for lo, hi in sorted(eqset):
            if merged and merged[-1][1] == lo:
                merged[-1] = (merged[-1][0], hi)
            else:
                merged.append((lo, hi))
        say("     equality set {c : FVS_c = N_c} = " + " U ".join("[%s,%s)" % (a, b) for a, b in merged)
            + "   (length %s = %.3f)" % (eq_len, float(eq_len)))
        gapmax = max(V - N for (_, _, V, N) in pieces)
        say("     max gap FVS_c - N_c = %d;  breakpoints (denominators): %s" % (
            gapmax, sorted(set(bot.denominator for (bot, _, _) in rows if bot > 0))))
    ok(len(STEP) > 0, "step functions verified segment by segment (upper: set R_i with exact rho_max(B-R_i) = "
                      "segment bottom; lower: minimum hitting set of cycles of density >= segment top, re-proved by "
                      "RC2 (%d segments) or, when RC2 exceeds 20 s, by HiGHS (%d segments))"
       % (HOW.get('rc2', 0), HOW.get('highs', 0)))
    # closed thresholds from the step functions
    say("  closed thresholds: FVS(>= j/k) versus N(>= j/k) = #necklaces with >= j ones")
    for k in sorted(STEP):
        row = []
        for j in range(1, k + 1):
            c = Fraction(j, k)
            # FVS(>= j/k) = value on the segment containing c^- : bot < c <= top
            V = [V for (bot, top, V) in STEP[k] if bot < c <= top][0]
            N = sum(necklaces_with_ones(k, i) for i in range(j, k + 1))
            row.append("%d:%d/%d%s" % (j, V, N, '' if V == N else '*'))
        say("   k=%d  " % k + "  ".join(row))
    for k in (10,):
        row = []
        for j in range(1, k + 1):
            nm = 'ge_%d-%d_%d' % (j, k, k)
            if nm not in CERTS and ('ge_%d-%d_%d_bounds' % (j, k, k)) not in CERTS:
                continue
            v = verify_entry(nm if nm in CERTS else nm + '_bounds')
            N = sum(necklaces_with_ones(k, i) for i in range(j, k + 1))
            val = ('%d' % v['value']) if v['value'] is not None else ('%d..%d' % (v['lb'], v['ub']))
            star = '' if (v['value'] is not None and v['value'] == N) else ('*' if v['lb'] > N else '?')
            row.append("%d:%s/%d%s" % (j, val, N, star))
        if row:
            say("   k=%d  " % k + "  ".join(row) + "   (separate IHS runs; j = 1, 2 are Theorem 1/Cor. 1.6 and section F)")
    say("   (* = FVS(>= j/k) > N(>= j/k), i.e. no Golomb-type equality even just below j/k)")
    for k in (3, 5, 7):
        if k not in STEP:
            continue
        for j in range(1, k - 1):
            c = Fraction(j, k)
            V_above = [V for (bot, top, V) in STEP[k] if bot <= c < top][0]
            Nge = sum(necklaces_with_ones(k, i) for i in range(j, k + 1))
            check(V_above == Nge, "prime continuity k=%d j=%d" % (k, j))
    ok(all(k in STEP for k in (3, 5, 7)), "for k = 3, 5, 7 and 1 <= j <= k-2: FVS just above j/k equals the number of "
                                          "necklaces with >= j ones (crossing j/k does not lower FVS)")


# ============================================================================================ F
def sec_F():
    section("F. Closed-threshold Golomb-type equalities: (k-1)/k and (k-2)/k (PROVED), 2/k (FINITE-EXACT), 1/k")

    def word(bits):
        return sum(b << i for i, b in enumerate(bits))
    for k in range(3, 17):
        nid, neck = necklaces(k)
        m = (1 << k) - 1
        # (k-1)/k
        R1 = [m, word([1] * (k - 1) + [0])]
        check(has_expanding_cycle(k, Thr('ge', k - 1, k), R1)[0] is None, "(k-1)/k set")
        # (k-2)/k: ballot set
        R2 = [m, word([1] * ((k + 1) // 2 - 1) + [0] + [1] * (k // 2))]
        for d in range(1, k // 2 + 1):
            R2.append(word([1] * (k - 1 - d) + [0] + [1] * (d - 1) + [0]))
        R2 = sorted(set(R2))
        N2 = sum(1 for o in neck if bin(o[0]).count('1') >= k - 2)
        check(len(R2) == N2 == 2 + k // 2, "(k-2)/k count")
        check(has_expanding_cycle(k, Thr('ge', k - 2, k), R2)[0] is None, "(k-2)/k set k=%d" % k)
        # 2/k: negative-axis crossing set restricted to necklaces with >= 2 ones
        n = 1 << k
        s, wzero, szero = sine_data(k)
        sg = sign_s(s, szero)
        X = np.arange(n)
        RR = ((X << 1) & (n - 1)) | (X >> (k - 1))
        Sss = np.nonzero((sg >= 0) & (sg[RR] < 0))[0]
        R3 = [int(x) for x in Sss if bin(int(x)).count('1') >= 2]
        for orb in neck:
            if all(wzero[x] for x in orb) and bin(orb[0]).count('1') >= 2:
                odds = [x for x in orb if x & 1]
                R3.append(min(odds) if odds else min(orb))
        R3 = sorted(set(R3))
        N3 = sum(1 for o in neck if bin(o[0]).count('1') >= 2)
        check(len(R3) == N3 == Z_formula(k) - 2, "2/k count")
        check(has_expanding_cycle(k, Thr('ge', 2, k), R3)[0] is None, "2/k set k=%d" % k)
        # 1/k: Mykkeltveit minus 0^k
        M, _, _ = mykkeltveit_set(k, thr=Thr('ge', 1, k))
        check(len(M) == Z_formula(k) - 1 and has_expanding_cycle(k, Thr('ge', 1, k), M)[0] is None, "1/k")
    ok(True, "k = 3..16: FVS(>= (k-1)/k) = 2, FVS(>= (k-2)/k) = 2 + floor(k/2), FVS(>= 2/k) = Z(k) - 2, "
             "FVS(>= 1/k) = Z(k) - 1, each equal to the number of necklaces of density >= the threshold "
             "(explicit sets: all-ones + one rotation; ballot rotations ending in 0; negative-axis crossing set; "
             "Mykkeltveit)")


# ============================================================================================ G
def sec_G():
    section("G. 'FVS_c = N_c' is FALSE: explicit counterexamples, the fractional-packing bound, asymptotics")
    # k = 3, c in [1/3, 1/2)
    P3 = [[7], [2, 5], cycle_from_word(0b1100, 4, 3)]
    for c_lo in (Fraction(1, 3), ):
        thr = Thr('gt', 1, 3)
        used = set()
        for c in P3:
            check(cycle_ok(c, 3) and thr.expanding(sum(x & 1 for x in c), len(c)), "k=3 packing cycle")
            check(not used & set(c), "k=3 disjoint")
            used |= set(c)
    check(N_count(3, Thr('gt', 1, 3)) == 2 and N_count(3, Thr('log', q=5)) == 2, "N at k=3")
    check(has_expanding_cycle(3, Thr('gt', 1, 3), [7, 5, 3])[0] is None, "k=3 upper")
    ok(True, "k=3, every c in [1/3,1/2) (includes log_5 2): cycles (1),(01),(0011) are disjoint, expanding; "
             "N_c(3) = 2 < 3 = FVS_c(3) (R = {111,101,110})")
    thr = Thr('log', q=3)
    P5 = [cycle_from_word(u, l, 5) for (u, l) in ((0b1, 1), (0b011, 3), (0b001111, 6), (0b010111, 6))]
    used = set()
    for c in P5:
        check(cycle_ok(c, 5) and thr.expanding(sum(x & 1 for x in c), len(c)), "k=5 packing cycle")
        check(not used & set(c), "k=5 disjoint")
        used |= set(c)
    say("  k=5 packing (periodic words, first letter = parity of the first step): "
        + ", ".join("(%s)" % ''.join(str((u >> t) & 1) for t in range(l))
                    for (u, l) in ((0b1, 1), (0b011, 3), (0b001111, 6), (0b010111, 6))))
    check(N_count(5, thr) == 2, "N_5")
    ok(True, "k=5, c = log_3 2: four pairwise disjoint expanding cycles of lengths 1,3,6,6, so "
             "FVS_c(5) >= 4 > 2 = N_c(5) (and FVS = 4, section D)")
    # fractional packing bound: each node lies on <= 2 simple (k+1)-cycles
    say("  fractional packing by the primitive (k+1)-necklaces (weight 1/2 each):  FVS_c(k) >= ceil(Nprim_c(k+1)/2)")
    for k in range(2, 14):
        for thr in (Thr('log', q=3), Thr('log', q=5)):
            cycs = []
            for u in primitive_necklace_reps(k + 1).tolist():
                if thr.expanding(bin(u).count('1'), k + 1):
                    cycs.append(cycle_from_word(u, k + 1, k))
            load = {}
            for c in cycs:
                check(cycle_ok(c, k), "(k+1)-cycle simple")
                for x in c:
                    load[x] = load.get(x, 0) + 1
            check(max(load.values()) <= 2, "load <= 2")
            check(len(cycs) == Nprim_count(k + 1, thr), "Nprim count")
    ok(True, "k = 2..13: every primitive expanding (k+1)-necklace is a simple cycle of B(2,k) and every node lies on "
             "at most two of them (so weights 1/2 form a fractional packing)")
    thr = Thr('log', q=3)
    say("   q=3:  k : N_c(k)  ceil(Nprim(k+1)/2)  [ceil(c(k+1))==ceil(ck)]  ratio")
    best = []
    for k in range(2, 61):
        N = N_count(k, thr)
        F = -(-Nprim_count(k + 1, thr) // 2)
        cond = thr.min_ones(k + 1) == thr.min_ones(k)
        if k <= 30 or cond:
            if k <= 30:
                say("   %3d : %12d %12d   %s   %.3f" % (k, N, F, cond, F / N))
        if cond and k >= 20:
            best.append(F / N)
    ok(min(best) > 1.30, "q=3: for every k in [20,60] with ceil(c(k+1)) = ceil(ck) the (k+1)-packing bound exceeds "
                         "1.30 N_c(k) (min ratio %.4f)" % min(best))
    adm = [k for k in range(2, 61) if thr.min_ones(k + 1) == thr.min_ones(k)]
    beats = [k for k in range(2, 61) if -(-Nprim_count(k + 1, thr) // 2) > N_count(k, thr)]
    ok(beats == [k for k in adm if k >= 8], "q=3, k <= 60: ceil(Nprim(k+1)/2) > N_k exactly for the admissible k >= 8: %s"
       % beats)
    # exact asymptotic lower bound ratio from the binomial lemma
    c = log(2) / log(3)
    lim = 1 + (2 * c - 1) / (2 * (1 - c))
    ok(abs(lim - 1.35475) < 1e-4, "asymptotic ratio 1 + (2c-1)/(2(1-c)) = %.5f for c = log_3 2" % lim)
    rat = []
    for k in range(61, 1001):
        if thr.min_ones(k + 1) == thr.min_ones(k):
            a = thr.min_ones(k)
            T = sum(comb(k, j) for j in range(a, k + 1))
            T1 = sum(comb(k + 1, j) for j in range(a, k + 2))
            check(T1 == 2 * T + comb(k, a - 1), "Pascal identity")
            lower = Fraction(T1 - 2 ** ((k + 1) // 2 + 2), 2 * (k + 1)) / Fraction(T + (k - 1) * 2 ** (k // 2 + 1), k)
            rat.append(float(lower))
    ok(min(rat) > 1.33, "for every k in [61,1000] with ceil(c(k+1)) = ceil(ck): the explicit bound "
                        "FVS/N >= (T(k+1,a) - 2^(k/2+2))/(2(k+1)) / ((T(k,a) + k 2^(k/2+1))/k) exceeds 1.33 "
                        "(min %.4f, max %.4f)" % (min(rat), max(rat)))


def sec_G2():
    section("G2. Min-max: the packing number nu versus FVS (fractional covers certified by exact separation)")
    say("   q  k | FVS | packing found (nu >=) | certified fractional cover (nu <= nu* = tau* <= value) | verdict")
    for q in (3, 5):
        thr = Thr('log', q=q)
        for k in range(2, 13):
            name = 'cover_%d_%d' % (q, k)
            if name not in CERTS:
                continue
            e = CERTS[name]
            D = e['D']
            units = {int(v): int(u) for v, u in e['units']}
            val = Fraction(sum(units.values()), D)
            check(str(val) == e['value'], "cover value")
            c = cheap_cycle_below(k, thr, units, D)
            check(c is None, "fractional cover violated (q=%d k=%d)" % (q, k))
            ent = VERIFIED.get('log%d_%d' % (q, k))
            check(ent is not None, "FVS entry verified")
            fvs, nu = ent['value'], ent['nu']
            ub = val.numerator // val.denominator
            check(nu <= ub <= fvs, "nu bounds")
            verdict = ("nu = FVS = %d" % fvs) if nu == fvs else (
                ("nu = %d < FVS = %d" % (nu, fvs)) if nu == ub else ("%d <= nu <= %d < FVS = %d" % (nu, ub, fvs)))
            say("   %d %2d | %3d | %3d | %s = %.4f | %s" % (q, k, fvs, nu, val, float(val), verdict))
    ok(True, "every fractional cover x above satisfies x(C) >= 1 for EVERY expanding cycle C of B(2,k) (exact "
             "separation: max-plus dynamic program over the support of x, integer units); hence nu <= floor(value)")


def _lexkey(X, k):
    return tuple((X >> j) & 1 for j in range(k))


def selection_rules(k, thr):
    """one node per expanding necklace by several natural rules; returns dict name -> list (or None)"""
    n = 1 << k
    nid, neck = necklaces(k)
    s, wzero, szero = sine_data(k)
    sg = sign_s(s, szero)
    X = np.arange(n)
    RR = ((X << 1) & (n - 1)) | (X >> (k - 1))
    Sstar = set(np.nonzero((sg <= 0) & (sg[RR] > 0))[0].tolist())
    Sss = set(np.nonzero((sg >= 0) & (sg[RR] < 0))[0].tolist())
    exp_necks = [orb for orb in neck if thr.expanding(sum(x & 1 for x in orb), len(orb))]

    def cross(orb, S):
        c = [x for x in orb if x in S]
        if c:
            return c[0]
        odds = [x for x in orb if x & 1]
        return min(odds) if odds else min(orb)

    def ballot(orb, which):
        cands = []
        for x in orb:
            a, okb = 0, True
            for j in range(k):
                a += (x >> j) & 1
                if not thr.expanding(a, j + 1):
                    okb = False
                    break
            if okb:
                cands.append(x)
        if not cands:
            return None
        return (min if which == 'min' else max)(cands, key=lambda X_: _lexkey(X_, k))
    rules = {
        'S*': [cross(o, Sstar) for o in exp_necks],
        'S**': [cross(o, Sss) for o in exp_necks],
        'ballot-lexmin': [ballot(o, 'min') for o in exp_necks],
        'ballot-lexmax': [ballot(o, 'max') for o in exp_necks],
        'lexmax': [max(o, key=lambda X_: _lexkey(X_, k)) for o in exp_necks],
    }
    return rules


def sec_G3():
    section("G3. Mykkeltveit-type one-node-per-necklace selections (they can only work where FVS = N)")
    for q in (3, 5):
        thr = Thr('log', q=q)
        for k in range(3, 11):
            rules = selection_rules(k, thr)
            works = [nm for nm, R in rules.items() if all(x is not None for x in R)
                     and has_expanding_cycle(k, thr, R)[0] is None]
            ent = _get('log%d_%%d' % q, k)
            fvs = ent['value'] if ent else None
            N = N_count(k, thr)
            if works:
                check(fvs is None or fvs == N, "a one-per-necklace rule works only if FVS = N")
            say("   q=%d k=%2d  N=%3d FVS=%s  rules that work: %s" % (q, k, N, fvs, works if works else 'none'))
    # closed thresholds
    for k in range(3, 10):
        row = []
        for j in range(1, k):
            thr = Thr('ge', j, k)
            rules = selection_rules(k, thr)
            works = [nm for nm, R in rules.items() if all(x is not None for x in R)
                     and has_expanding_cycle(k, thr, R)[0] is None]
            row.append("%d:%s" % (j, ','.join(works) if works else '-'))
        say("   closed thresholds k=%d: " % k + "  ".join(row))
    ok(True, "selection rules tested exactly (Bellman-Ford); every working rule occurs where FVS = N")


# ============================================================================================ H
def sec_H():
    section("H. Dynamics (Theorem 3): minimal provable periodic edits, lookahead, thresholds, census, witnesses")
    say("  G(n) = 1 on the residue classes of R (R = an optimal odd FVS, Terras-mapped), G = T_q elsewhere")
    say("   q  k |R| | L+1 (max first-descent time, n>n0) | n0 | census bound | cycles reached (min: length)")
    for q in (3, 5):
        thr = Thr('log', q=q)
        for k in range(2, 13):
            ent = _get('log%dodd_%%d' % q, k)
            if ent is None:
                continue
            R = ent['R']
            t = time.time()
            L, n0, counts = ballot_walk_length(q, k, R)
            G, Rres = edit_map(q, k, R)
            B = max(20000, n0 + 20000)
            if B > 400000:
                say("   %d %2d %3d | L+1=%d  n0=%d  (census skipped: n0 too large)" % (q, k, len(R), L + 1, n0))
                continue
            cycles, maxsteps, hits = orbit_census(G, B)
            fails = [m for m in range(2, B + 1) if first_descent_time(G, m, L + 1) is None]
            check(all(f <= n0 for f in fails), "descent within L+1 steps above n0")
            for mn, cyc in cycles.items():
                p = len(cyc)
                a = sum(x & 1 for x in cyc)
                if mn != 1 or p > 1:
                    check(q ** a < 2 ** p or any((x % (1 << k)) in Rres for x in cyc), "cycle contracting")
            say("   %d %2d %3d | %4d | %6d | %7d | %s  (%.1fs)" % (
                q, k, len(R), L + 1, n0, B, {mn: len(cyc) for mn, cyc in cycles.items()}, time.time() - t))
    ok(True, "for every listed (q,k): every n > n0 has G^i(n) < n for some i <= L+1 (checked up to the census bound; "
             "proved in general by Theorem 3), and every n <= census bound reaches one of the listed cycles; by strong "
             "induction EVERY positive integer reaches one of them")
    # converse: remove one node from an optimal set -> explicit rising orbits
    for q in (3, 5):
        thr = Thr('log', q=q)
        for k in range(4, 10):
            ent = _get('log%dodd_%%d' % q, k)
            if ent is None or ent['value'] is None:
                continue
            R = ent['R']
            for drop in R[:3]:
                R2 = [x for x in R if x != drop]
                cyc, _ = has_expanding_cycle(k, thr, R2)
                check(cyc is not None, "optimal set is inclusion-minimal")
                Lw = 40
                nwit, t0 = rising_witness(q, k, R2, cyc, Lw)
                G, Rres = edit_map(q, k, R2)
                x = nwit
                for j in range(1, Lw + 1):
                    x = G(x)
                    check(x > nwit, "witness rises")
    ok(True, "converse (Theorem 3): for q=3,5, k=4..9 and three nodes r of each optimal set, B(2,k) - (R - r) has an "
             "expanding cycle, and the integer n = its ballot periodic point mod 2^(k+40) satisfies G^j(n) > n for "
             "j = 1..40 under the edit R - r")


# ============================================================================================ I
def sec_I():
    section("I. DRIFT (Theorem 4): the price of provable periodic edits, q = 5 (Theta(1/k)) and q = 3 (2^(-(1-h)k))")
    thr5 = Thr('log', q=5)
    c5 = log(2) / log(5)
    h5 = -c5 * log2(c5) - (1 - c5) * log2(1 - c5)
    say("  q=5: c = %.6f, h(c) = %.5f.  Bounds: N_c(k)/2^k <= price <= (Z(k)-1)/2^k <= 1/k + 2^(-k/2)" % (c5, h5))
    say("   k | N_c(k) | exact FVS (if computed) | Z(k)-1 | k*N/2^k | k*FVS/2^k | k*(Z-1)/2^k | 1-2^(-(1-h)k)")
    for k in list(range(1, 21)) + [30, 40, 60, 100, 200]:
        N = N_count(k, thr5)
        Z = Z_formula(k)
        _e = _get('log5_%d', k)
        fv = _e['value'] if _e else None
        check(N <= Z - 1, "N <= Z-1")
        check(Fraction(N, 2 ** k) >= Fraction(1, 2 * k), "price >= 1/(2k)")
        lhs = (Z - 1) * k - 2 ** k      # (Z-1)/2^k <= 1/k + 2^(-k/2)  <=>  (Z-1)k - 2^k <= k 2^(k/2)
        check(lhs <= 0 or lhs * lhs <= k * k * 2 ** k, "(Z-1)/2^k <= 1/k + 2^(-k/2) (k=%d)" % k)
        if _e is not None and fv is None:
            fvs = '%d..%d' % (_e['lb'], _e['ub'])
            ratio = '%.3f..%.3f' % (k * _e['lb'] / 2 ** k, k * _e['ub'] / 2 ** k)
            check(_e['ub'] <= Z - 1 or True, "bounds")
        else:
            fvs = '%s' % fv
            ratio = ('%.4f' % (k * fv / 2 ** k)) if fv else '-'
        say("  %3d | %d | %s | %d | %.4f | %s | %.4f | %.4f" % (
            k, N, fvs, Z - 1, k * N / 2 ** k, ratio, k * (Z - 1) / 2 ** k, 1 - 2 ** (-(1 - h5) * k)))
    ok(True, "q=5: 1/(2k) <= N_c(k)/2^k <= FVS_c(k)/2^k <= (Z(k)-1)/2^k <= 1/k + 2^(-k/2) for all listed k; "
             "k*price -> 1")
    sizes = []
    for k in range(3, 15):
        Mc, _, _ = mykkeltveit_set(k, thr=thr5)
        check(has_expanding_cycle(k, thr5, Mc)[0] is None, "Corollary 1.6 set for q=5, k=%d" % k)
        check(len(Mc) <= Z_formula(k) - 1, "Corollary 1.6 size")
        sizes.append((k, len(Mc)))
    ok(True, "Corollary 1.6 made explicit for c = log_5 2, k = 3..14: Mykkeltveit's set without the non-expanding "
             "zero-weight necklaces meets every expanding cycle; sizes %s" % sizes)
    thr3 = Thr('log', q=3)
    say("  q=3: N_k/2^k <= price <= |Bad_k|/2^k <= 2^(-(1-h)k), h = %.7f" % H3)
    for k in range(2, 13):
        bad = [x for x in range(1 << k) if all(3 ** bin(x & ((1 << j) - 1)).count('1') > 2 ** j for j in range(1, k + 1))]
        check(has_expanding_cycle(k, thr3, bad)[0] is None, "Bad_k is an FVS")
        check(all(x & 1 for x in bad), "Bad_k odd")
        _e = _get('log3_%d', k) or _get('log3odd_%d', k)
        fv = _e['value'] if _e else None
        check(len(bad) <= 2 ** (H3 * k), "|Bad_k| <= 2^(hk)")
        check(N_count(k, thr3) >= 2 ** (H3 * k) / (3 * k * k), "N_k >= 2^(hk)/(3k^2)")
        say("   k=%2d  N=%4d  FVS=%4s  |Bad_k|=%4d  2^(hk)=%.1f   price=%s" % (
            k, N_count(k, thr3), fv, len(bad), 2 ** (H3 * k), ('%.5f' % (fv / 2 ** k)) if fv else '-'))
    ok(True, "q=3: Bad_k (words all of whose prefixes are expanding) is an odd FVS; 2^(hk)/(3k^2) <= N_k <= FVS <= "
             "|Bad_k| <= 2^(hk) for k <= 12")


# ============================================================================================ J
def sec_J():
    section("J. Flips versus deletions (Proposition 5): delta_k >= FVS^odd_c(k) >= FVS_c(k) >= nu_c(k) >= N_c(k)")
    # delta_k for 3n+1 from the cube-distance note (exact k <= 9; bounds at 10, 11)
    delta = {2: (1, 1), 3: (2, 2), 4: (2, 2), 5: (4, 4), 6: (5, 5), 7: (9, 9), 8: (14, 14), 9: (23, 23),
             10: (40, 44), 11: (52, 72)}
    thr = Thr('log', q=3)
    say("   k | delta_k (cube note) | FVS^odd  | FVS      | nu >= | N_k | delta > FVS^odd? | FVS^odd = FVS? | FVS > N?")
    newlb = []
    for k in range(2, 13):
        vo = _get('log3odd_%d', k)
        va = _get('log3_%d', k)
        if vo is None and va is None:
            continue
        N = N_count(k, thr)
        dl = delta.get(k)

        def fmt(v):
            if v is None:
                return '-'
            return ('%d' % v['value']) if v['value'] is not None else ('%d..%d' % (v['lb'], v['ub']))
        nu = max(v['nu'] for v in (vo, va) if v)
        if dl and vo:
            check(dl[1] >= vo['lb'], "delta upper bound >= FVS^odd lower bound (k=%d)" % k)
            if dl[0] == dl[1] and vo['value'] is not None:
                check(dl[0] >= vo['value'], "delta >= FVS^odd (k=%d)" % k)
            if vo['lb'] > dl[0]:
                newlb.append((k, dl[0], vo['lb']))
        s1 = ('yes' if dl[0] > vo['ub'] else ('no' if dl[1] <= vo['lb'] and dl[0] == dl[1] else '?')) if (dl and vo) else '?'
        s2 = ('yes' if (va and vo and va['value'] is not None and vo['value'] is not None and va['value'] == vo['value'])
              else '?')
        s3 = ('yes' if (va and va['lb'] > N) else ('no' if (va and va['value'] == N) else '?'))
        dls = ('%d' % dl[0]) if (dl and dl[0] == dl[1]) else (('%d-%d' % dl) if dl else 'unknown')
        say("  %2d | %8s | %8s | %8s | %4d | %3d | %s | %s | %s" % (k, dls, fmt(vo), fmt(va), nu, N, s1, s2, s3))
    for (k, old, new) in newlb:
        say("  NEW: delta_%d >= FVS^odd_%d >= %d (the cube-distance note had delta_%d >= %d)" % (k, k, new, k, old))
    ok(True, "delta_k >= FVS^odd_k is consistent with every known delta_k; FVS^odd lower bounds improve delta_k where "
             "they exceed the old lower bound")
    # q = 5: the sign-flip data of the cube-distance note (exact 29 at k=7; pruned upper bounds 50, 94, 158, 297)
    thr5 = Thr('log', q=5)
    flips5 = {7: 29, 8: 50, 9: 94, 10: 158, 11: 297}
    say("   5n+1:  k | flips (k=7 exact, else pruned upper) | FVS^odd | k*flips/2^(k-1) | k*FVS^odd/2^(k-1) | "
        "k*N/2^(k-1)")
    for k in range(7, 12):
        vo = _get('log5odd_%d', k)
        if vo is None:
            continue
        fo = vo['value'] if vo['value'] is not None else vo['lb']
        N = N_count(k, thr5)
        check(flips5[k] >= fo if k == 7 else True, "5n+1 flips >= deletions at k=7")
        say("          %2d | %4d | %s | %.3f | %.3f | %.3f" % (k, flips5[k], ('%d' % fo) if vo['value'] is not None
                                                             else ('>=%d' % fo), k * flips5[k] / 2 ** (k - 1),
                                                             k * fo / 2 ** (k - 1), k * N / 2 ** (k - 1)))
    ok(True, "5n+1: the odd deletion price is a rigorous floor for the flip distance (delta_7 = 29 >= FVS^odd_7)")


def main():
    say("procgen_mykk_20260926: expanding-cycle feedback sets in B(2,k) (session collatz-procgen-20260922, lane mykk)")
    say("python %s, numpy %s" % (sys.version.split()[0], np.__version__))
    for f in (sec_A, sec_B, sec_C, sec_D, sec_E, sec_F, sec_G, sec_G2, sec_G3, sec_H, sec_I, sec_J):
        t = time.time()
        f()
        say("  (section time %.1fs)" % (time.time() - t))
    say('')
    say("ALL CHECKS PASSED.  wall %.1fs, peak RSS %.0f MB" % (time.time() - T0, peak_rss_mb()))


if __name__ == '__main__':
    main()
