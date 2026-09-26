#!/usr/bin/env python3
"""procgen_sporadic_20260926_run.py -- runner of the lane "sporadic" (free and sporadic cycles;
the Belaga-Mignotte off-by-one).  Session collatz-procgen-20260922, 2026-09-26.

Usage:  python3 -u 04-computation/experiments/procgen_sporadic_20260926_run.py [--quick]
        > 05-knowledge/results/procgen_sporadic_20260926.out

Every printed claim is a check(...) that raises on failure.  --quick skips the long C sweeps
(T2.6-T2.8); the recorded .out is a full run.

Parts:
  T1  free / partially free / sporadic words and cycles of T_{q,d}(y) = y/2, (q y + d)/2
  T2  the Belaga-Mignotte counts omega(14303) = 944, omega(17021) = 258 (the two long cycles)
  T3  Stern-Brocot positions of the known cycles
  T4  5x+1 and 3x-1
C programs (compiled here with cc -O2): procgen_sporadic_20260926_traj.c (least-element search for one
map), procgen_sporadic_20260926_sweep.c (least-element search over many d).
"""
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_sporadic_20260926_lib import (  # noqa: E402
    HERE, ROOT, SCRATCH, N_CHECKS, check, peak_rss_mb, child_peak_rss_mb, sha256_file, carry, gap,
    Tmap, periodic_point, parity_word, words_of_shape, primitive_period, lyndon_count, necklace_rep,
    cycle_of, cycle_data, gersonides_solutions, gersonides_predicted, small_gaps_q3, below_theta,
    sb_path, on_sb_path, ordmod, lattice_basis, cone_min_clock, compile_c, run_parallel, parse_sweep,
    parse_traj)

QUICK = "--quick" in sys.argv
T0 = time.time()
TRAJ = os.path.join(SCRATCH, "traj")
SWEEP = os.path.join(SCRATCH, "sweep")

BM_TABLE = {7463: 162, 18359: 164, 7727: 198, 15655: 207, 10289: 214, 9823: 241, 17021: 258,
            14197: 329, 13085: 335, 6487: 534, 14303: 944}
BM_TOTAL, BM_N1, BM_N2, BM_N3 = 42765, 1481, 1507, 1005


def header(s):
    print()
    print("=" * 110)
    print(s)
    print("=" * 110, flush=True)


def traj_run(q, d, ylo, yhi, tag, steplim=3000):
    out = os.path.join(SCRATCH, f"traj_{tag}.txt")
    run_parallel([[TRAJ, str(q), str(d), str(ylo), str(yhi), "1", "0", str(steplim)]], [out], 1)
    return parse_traj([out])


def cycles_on_Z(q, d, X, tag):
    """cycles of T_{q,d} on Z with least |element| <= X, from the least-element search on the positive
    integers for d and for -d (T_{q,d}(-y) = -T_{q,-d}(y)).  Returns list of (sign, m, p, a)."""
    out = [(0, 0, 1, 0)]
    for sgn, dd in ((1, d), (-1, -d)):
        cyc, ent, bad, n = traj_run(q, dd, 1, X, f"{tag}_{'p' if sgn > 0 else 'n'}")
        assert not bad, bad[:3]
        for (m, p, a, mx) in cyc:
            out.append((sgn, m, p, a))
    return out


# ================================================================================================
# T1
# ================================================================================================
def t1():
    header("T1. Free, partially free and sporadic words of T_{q,d}(y) = y/2, (q y + d)/2  (q odd, d odd, "
           "gcd(q,d) = 1)")
    rnd = random.Random(20260926)

    # T1.1 affine form and itinerary of y_w
    n = 0
    for q in (3, 5, 7, 9):
        for d in (1, -1, 5, -7, 13, 25, -139, 1001):
            if math.gcd(q, d) != 1:
                continue
            for _ in range(60):
                p = rnd.randint(1, 40)
                w = tuple(rnd.randint(0, 1) for _ in range(p))
                a = sum(w)
                D = gap(p, a, q)
                c = carry(w, q)
                # y_w mod 2^p (D odd, invertible mod 2^p); every y in that class has itinerary w
                r = (d * c * pow(D, -1, 1 << p)) % (1 << p)
                for k in range(3):
                    y = r + (rnd.randint(-10 ** 6, 10 ** 6) << p)
                    ww, yp = parity_word(y, q, d, p)
                    assert ww == w and (yp << p) == q ** a * y + d * c
                    n += 1
    check(n > 0, f"T1.1 affine form: 2^p T^p(y) = q^a y + d c_w and itinerary w on the class y = d c_w/D "
                 f"(mod 2^p), {n} random (q, d, w, y), q in {{3,5,7,9}}, p <= 40")

    # T1.2 shift criterion (Proposition 1), exhaustive on small shapes
    tested = free_found = 0
    for q in (3, 5, 7):
        for d in (1, -1, 5, -5, 7, -7, 11, 13, -13, 23, 25, -25, 47, 115, -139, 485, 455):
            if math.gcd(q, d) != 1:
                continue
            for p in range(2, 13):
                for a in range(1, p):
                    D = gap(p, a, q)
                    allint = all((d * carry(w, q)) % D == 0 for w in words_of_shape(p, a))
                    assert allint == (d % D == 0), (q, d, p, a)
                    u = tuple([1] * a + [0] * (p - a))
                    u2 = tuple([1] * (a - 1) + [0, 1] + [0] * (p - a - 1))
                    assert carry(u2, q) - carry(u, q) == 1 << (a - 1)
                    tested += 1
                    free_found += allint
    check(tested > 0 and free_found > 0,
          f"T1.2 shift criterion: for every mixed shape (1 <= a <= p-1, p <= 12), q in {{3,5,7}} and the d "
          f"in {{+-1,+-5,+-7,11,+-13,23,+-25,47,115,-139,455,485}} prime to q, 'every word integral' <=> D | d ({tested} (q,d,shape) cases, {free_found} free); the two "
          f"words 1^a 0^(p-a) and 1^(a-1) 0 1 0^(p-a-1) differ in carry by exactly 2^(a-1)")

    # T1.3 single-word shapes
    ok = True
    for q in (3, 5, 7, 9, 11):
        for d in (1, -1, 3, 5, 7, 9, 15, 21, 35):
            if math.gcd(q, d) != 1:
                continue
            for p in range(1, 12):
                ok &= periodic_point((0,) * p, q, d) == 0
                ok &= periodic_point((1,) * p, q, d) == Fraction(-d, q - 2)
    check(ok, "T1.3 single-word shapes: 0^p gives 0, 1^p gives -d/(q-2) (integral iff (q-2) | d, i.e. "
              "D(1,1) = 2 - q divides d); for q = 3 this is the fixed point -d")

    # T1.4 elementary Gersonides for every odd q
    bad = [q for q in range(3, 1000, 2) if gersonides_solutions(q, 64) != gersonides_predicted(q, 64)]
    check(not bad, "T1.4 |2^p - q^a| = 1 (p, a >= 1): exactly a = 1 with q = 2^p +- 1, and (q,p,a) = (3,3,2); "
                   "checked for all odd q < 1000, p <= 64 (hand proof in the note)")
    mixed = {q: [s for s in gersonides_solutions(q, 64) if s[1] <= s[0] - 1] for q in (3, 5, 7, 9, 15, 17)}
    check(mixed[3] == [(2, 1), (3, 2)] and mixed[5] == [(2, 1)] and mixed[7] == [(3, 1)] and
          mixed[9] == [(3, 1)] and mixed[15] == [(4, 1)] and mixed[17] == [(4, 1)],
          f"T1.4 free mixed shapes of T_(q,+-1): q=3 {mixed[3]}, q=5 {mixed[5]}, q=7 {mixed[7]}, q=9 {mixed[9]}, "
          f"q=15 {mixed[15]}, q=17 {mixed[17]}")

    # T1.5 complete list of small gaps for q = 3 (Ellison 1971)
    N = 20000
    lst, x0 = small_gaps_q3(N)
    brute = [(p, a) for p in range(1, 401) for a in range(1, 400) if abs((1 << p) - 3 ** a) <= N]
    check(x0 == 17 and lst == brute,
          f"T1.5 all (p,a), p,a >= 1, with |2^p - 3^a| <= {N}: Ellison's bound 2^x e^(-x/10) > {N} for x >= "
          f"{x0} leaves x <= 16 or x in S = {{19, 27}}; the {len(lst)} solutions agree with brute force to p = 400")
    print("      |2^p - 3^a| <= 20000, mixed shapes (a <= p-1):")
    row = [f"({p},{a}):{(1 << p) - 3 ** a}" for (p, a) in lst if a <= p - 1]
    for i in range(0, len(row), 10):
        print("        " + "  ".join(row[i:i + 10]))

    # T1.6 free shapes and free cycles for a few d, confirmed by the least-element search on Z
    X = 10 ** 6
    table = []
    ok16 = True
    for d in (1, 5, 7, 11, 13, 23, 25, 29, 47, 139):
        free = sorted(((p, a) for (p, a) in lst if a <= p - 1 and d % abs((1 << p) - 3 ** a) == 0))
        # predicted free cycles: {0}, {-d}, and every Lyndon word of every free mixed shape
        pred = set()
        for (p, a) in free:
            for w in words_of_shape(p, a):
                if primitive_period(w) == p and w == necklace_rep(w):
                    y = periodic_point(w, 3, d)
                    assert y.denominator == 1
                    orb = cycle_of(int(y), 3, d)
                    pred.add((1 if y > 0 else -1, min(abs(v) for v in orb), p, a))
        pred.add((0, 0, 1, 0))
        pred.add((-1, d, 1, 1))
        found = set(cycles_on_Z(3, d, X, f"t1_{d}"))
        ok16 &= pred <= found
        spor = found - pred
        table.append((d, free, len(pred), sorted(spor)))
        # every other cycle is on a non-free clock
        ok16 &= all(d % gap(p, a, 3) != 0 for (s, m, p, a) in spor)
    check(ok16, f"T1.6 free cycles of T_(3,d) on Z for d in {{1,5,7,11,13,23,25,29,47,139}}: the predicted set "
                f"{{0}}, {{-d}} and one cycle per Lyndon word of each free mixed shape (D | d) is contained in "
                f"the least-element census (|least| <= {X:.0e}, both signs), and every other census cycle lies on "
                f"a clock with D not dividing d")
    for d, free, npred, spor in table:
        print(f"      d = {d:>3}: free mixed shapes {free}; free cycles {npred}; others (sign, least, p, a): "
              + (", ".join(str(t) for t in spor[:8]) + (" ..." if len(spor) > 8 else "") if spor else "none"))

    # T1.7 word-level partition by primitive level f (Proposition 3) on small clocks
    cases = 0
    for d in (5 * 7 * 13, 5 * 5 * 7, 11 * 13 * 23, 139 * 5, 23 * 47):
        for p in range(2, 15):
            for a in range(1, p):
                D = gap(p, a, 3)
                g = math.gcd(D, d)
                Dp = D // g
                integral = [w for w in words_of_shape(p, a) if (d * carry(w, 3)) % D == 0]
                assert all(carry(w, 3) % Dp == 0 for w in integral)
                # partition by f = d / gcd(y_w, d): f | gcd(d, D), (D/f) | c_w, gcd(c_w f / D, f) = 1
                for w in integral:
                    y = d * carry(w, 3) // D
                    f = d // math.gcd(y, d)
                    assert g % f == 0 and carry(w, 3) % (D // f) == 0
                    assert math.gcd(abs(carry(w, 3) * f // D), f) == 1
                    assert y == (d // f) * (f * carry(w, 3) // D)
                cases += 1
    check(cases > 0, f"T1.7 on {cases} (d, clock) cases (composite d, p <= 14): the integral words are exactly "
                     f"those with (D/g) | c_w, g = gcd(D, d), and each is the (d/f)-multiple of a primitive "
                     f"T_(3,f)-periodic point for the unique f = d/gcd(y_w, d), which divides g")

    # T1.8 illustration: the three kinds for T_(3,5)
    d = 5
    kinds = []
    for (s, m, p, a) in sorted(set(cycles_on_Z(3, d, X, "t1_kinds"))):
        D = gap(p, a, 3)
        g = math.gcd(D, d)
        kind = "free" if d % D == 0 else ("partially free" if g > 1 else "sporadic")
        y = s * m
        f = d // math.gcd(y, d) if y else 1
        kinds.append((s * m, p, a, D, g, kind, f))
    got = {(y, kind) for (y, p, a, D, g, kind, f) in kinds}
    check({(-85, "sporadic"), (187, "partially free"), (347, "partially free"), (19, "free"), (23, "free"),
           (1, "free"), (5, "free"), (-5, "free"), (-25, "free"), (0, "free")} == got,
          "T1.8 T_(3,5) on Z (census |least| <= 1e6): free {0}, {-5}, {5,10}, {-25,-35,-50}, {1,4,2}, {19,..}, "
          "{23,..}; partially free {187,..}, {347,..} on (27,17) (g = 5, D/g = 1015513); sporadic 5*{-17,..} "
          "on (11,7) (g = 1)")
    for row in kinds:
        print("      least %6d  clock (%d,%d)  D = %d  g = %d  %-15s primitive level f = %d" % row)


# ================================================================================================
# T2
# ================================================================================================
LONG = {14303: 101, 17021: 5}


def window(d, p, a):
    import mpmath
    mpmath.mp.dps = 50
    G = gap(p, a, 3)
    lo = -(-(d * 3 ** (a - 1)) // G)
    hi = int(mpmath.floor(mpmath.mpf(d) / (mpmath.power(2, mpmath.mpf(p) / a) - 3)))
    return lo, hi


def clock_census_np(d, p, a, chunk=1 << 18):
    """independent per-clock count (gates-lane method): odd y in the perigee window, iterate p steps,
    keep y with first return at step p, a odd steps, y minimal and gcd(y, d) = 1."""
    import numpy as np
    lo, hi = window(d, p, a)
    found = []
    LIM = 1 << 61
    y0 = max(1, lo) | 1
    while y0 <= hi:
        y1 = min(hi, y0 + 2 * chunk)
        Y0 = np.arange(y0, y1 + 1, 2, dtype=np.int64)
        Y = Y0.copy()
        odd = np.zeros(len(Y), dtype=np.int64)
        first = np.zeros(len(Y), dtype=np.int64)
        alive = np.ones(len(Y), dtype=bool)
        for s in range(1, p + 1):
            b = Y & 1
            odd += b
            Y = np.where(b == 1, (3 * Y + d) >> 1, Y >> 1)
            alive &= (Y >= Y0) & (Y < LIM)
            Y = np.where(alive, Y, Y0 + 1)
            hit = alive & (Y == Y0) & (first == 0)
            first[hit] = s
        ok = alive & (first == p) & (odd == a)
        found.extend(int(v) for v in Y0[ok] if math.gcd(int(v), d) == 1)
        y0 = y1 + 2
    return found, lo, hi


def t2():
    header("T2. The Belaga-Mignotte counts omega(14303) = 944 and omega(17021) = 258")

    # T2.1 the two long cycles, exactly
    for d, m in LONG.items():
        orb = cycle_of(m, 3, d)
        p, a, mm, w, D = cycle_data(orb, 3, d)
        u, v, idx, e, a0 = lattice_basis(d)
        # coordinates of (p, a) in the reduced basis
        det = u[0] * v[1] - u[1] * v[0]
        i = Fraction(p * v[1] - a * v[0], det)
        j = Fraction(u[0] * a - u[1] * p, det)
        lo, hi = window(d, p, a)
        check(mm == m and primitive_period(w) == p and all(math.gcd(x, d) == 1 for x in orb) and D > 0 and
              D % d == 0 and m * D == d * carry(w, 3) and lo <= m <= hi and i.denominator == 1 and
              j.denominator == 1,
              f"T2.1 d = {d}: the orbit of {m} is a primitive T_d-cycle with least element {m}, period p = {p}, "
              f"a = {a} odd steps (a/p = {a / p:.4f}), max element {max(orb)}; D = 2^{p} - 3^{a} is divisible "
              f"by d (D/d has {(D // d).bit_length()} bits), m D = d c_w; perigee window [{lo}, {hi}]; "
              f"(p,a) = {int(i)}*{u} + {int(j)}*{v} in the reduced basis of the clock lattice (index {idx})")

    # T2.2 the main clocks
    G27 = gap(27, 17, 3)
    G65 = gap(65, 41, 3)
    check(G27 == 5 * 71 * 14303 and G65 % 17021 == 0 and G65 // 17021 == 19 * 29 * 44835377399,
          "T2.2 2^27 - 3^17 = 5077565 = 355 * 14303 (355 = 5*71) and 2^65 - 3^41 = 420491770248316829 = "
          "17021 * 19 * 29 * 44835377399; 14303 and 17021 are prime")
    from sympy import isprime
    check(isprime(14303) and isprime(17021) and isprime(44835377399),
          "T2.2 primality: 14303, 17021, 44835377399 are prime (sympy)")

    # T2.3 what the gates lane scanned, and where the cycles are (least-element search, both d)
    base = {}
    for d in LONG:
        X = 1200 * d
        cyc, ent, bad, n = traj_run(3, d, 1, X, f"t2_{d}")
        check(not bad and n == (X + 1) // 2,
              f"T2.3 d = {d}: least-element search over all {n} odd y <= {X}: no overflow / unresolved start")
        prim = [(m, p, a) for (m, p, a, mx) in cyc if math.gcd(m, d) == 1]
        non = [(m, p, a) for (m, p, a, mx) in cyc if math.gcd(m, d) != 1]
        byclock = Counter((p, a) for (m, p, a) in prim)
        base[d] = (prim, byclock)
        check(len(prim) == BM_TABLE[d] and non == [(d, 2, 1)],
              f"T2.3 d = {d}: {len(prim)} primitive cycles with least element <= {X} (= Belaga-Mignotte "
              f"omega(d)) and one non-primitive, d*{{1,2}}; per clock: "
              + ", ".join(f"{k}: {c}" for k, c in sorted(byclock.items())))
        (bp, ba), e = cone_min_clock(d, X)
        check((1 << bp) > 3 ** ba and pow(2, bp, d) == pow(3, ba, d) and bp > 2200,
              f"T2.3 d = {d}: every lattice clock whose perigee window reaches {X + 1} has p >= {bp} (least such clock "
              f"({bp},{ba})), so these are all primitive T_d-cycles with period p < {bp}")
        # gates-lane scan domain: lattice clocks p <= 250 and the main family to p <= 600
        u, v, idx, e2, a0 = lattice_basis(d)
        clocks = []
        t3a, bl = 1, 1
        for a in range(1, 1400):
            t3a = t3a * 3
            bl = t3a.bit_length()          # least p with 2^p > 3^a
            r3 = t3a % d
            for p in range(max(bl, a + 1), 2201):
                if pow(2, p, d) == r3:
                    clocks.append((p, a))
        clocks.sort()
        main = u
        scanned = [c for c in clocks if c[0] <= 250 or (c[0] % main[0] == 0 and c[1] * main[0] == c[0] * main[1]
                                                          and c[0] <= 600)]
        unscanned_with = [c for c in byclock if c not in scanned]
        check(unscanned_with == [max(byclock, key=lambda c: c[0])] and byclock[unscanned_with[0]] == 1,
              f"T2.3 d = {d}: of the {len(clocks)} lattice clocks with p <= 2200 the gates lane scanned "
              f"{len(scanned)} (p <= 250, and the main family {main} to p <= 600); the only clock carrying a "
              f"cycle outside that set is {unscanned_with[0]} with exactly one cycle")
        rows = []
        for c in clocks:
            lo, hi = window(d, *c)
            if hi >= 1:
                pred = lyndon_count(*c) * d / gap(c[0], c[1], 3) * (1 - 1 / d)
                rows.append((c, byclock.get(c, 0), c in scanned, hi, pred))
        s_uns = sum(r[4] for r in rows if not r[2])
        s_sc = sum(r[4] for r in rows if r[2])
        check(len(rows) > 0, f"T2.3 d = {d}: {len(rows)} lattice clocks with p <= 2200 have a nonempty window; "
                             f"equidistribution prediction L(p,a) (d/D)(1-1/d) summed over the scanned ones "
                             f"{s_sc:.1f}, over the unscanned ones {s_uns:.3f}")
        print(f"      d = {d}: lattice clocks p <= 2200 with a cycle or prediction >= 0.01 (clock, cycles, "
              f"scanned by gates, window top, prediction):")
        for r in rows:
            if r[1] or r[4] >= 0.01:
                print("        (%4d,%4d)  %4d  %-5s  %9d  %10.3g" % (r[0][0], r[0][1], r[1], r[2], r[3], r[4]))

    # T2.4 independent per-clock check (window method, numpy) on the clocks carrying cycles
    for d in LONG:
        prim, byclock = base[d]
        ok24 = True
        for c, cnt in sorted(byclock.items()):
            f, lo, hi = clock_census_np(d, *c)
            got = sorted(f)
            want = sorted(m for (m, p, a) in prim if (p, a) == c)
            ok24 &= got == want
        check(ok24, f"T2.4 d = {d}: the perigee-window method (every odd y in the window, p steps) reproduces "
                    f"the least-element search clock by clock: " + ", ".join(f"{k}: {v}" for k, v in
                                                                          sorted(byclock.items())))

    # T2.5 the other nine table (20) entries at the same bound
    for d in sorted(BM_TABLE):
        if d in LONG:
            continue
        X = 1200 * d
        cyc, ent, bad, n = traj_run(3, d, 1, X, f"t2_{d}")
        prim = [(m, p, a) for (m, p, a, mx) in cyc if math.gcd(m, d) == 1]
        (bp, ba), e = cone_min_clock(d, X)
        check(not bad and len(prim) == BM_TABLE[d],
              f"T2.5 d = {d}: {len(prim)} primitive cycles with least element <= 1200 d = omega(d); complete for "
              f"periods p < {bp}; longest p = {max(p for (m, p, a) in prim)}")

    if QUICK:
        print("  (--quick: T2.6-T2.8 skipped)")
        return

    # T2.6 the full sweep d <= 19999, least element <= 1200 d
    outs = [os.path.join(SCRATCH, f"sweep_base_s{s}.txt") for s in (1, 3)]
    t = time.time()
    run_parallel([[SWEEP, "1", "19999", "1200", "4", str(s)] for s in (1, 3)], outs, 2)
    summ, cyc, other = parse_sweep(outs)
    t_base = time.time() - t
    ds = [d for d in range(1, 20000) if d % 2 and d % 3]
    check(sorted(summ) == ds and not other and all(len(v) == 1 and v[0][4] == 0 and v[0][5] == 0
                                                   for v in summ.values()),
          f"T2.6 sweep: all {len(ds)} d <= 19999 with gcd(d,6) = 1, every odd y <= 1200 d tested "
          f"({t_base:.0f} s wall, 2 processes); no overflow, no unresolved start, no orbit entering a cycle "
          f"with least element > 1200 d")
    om = {d: summ[d][0][2] for d in ds}
    tot = sum(om.values())
    dist = Counter(om.values())
    # independent re-verification of every cycle in Python
    ok26 = True
    for (d, m, p, a, pr) in cyc:
        orb = cycle_of(m, 3, d, maxlen=p)
        ok26 &= orb is not None and len(orb) == p and sum(x & 1 for x in orb) == a and min(orb) == m
        ok26 &= (math.gcd(m, d) == 1) == bool(pr)
    check(ok26 and len(cyc) > 0, f"T2.6 all {len(cyc)} cycles re-verified by exact Python iteration (period, odd steps, least "
                f"element, primitivity)")
    # multiplicativity (Proposition 3): cycles of T_d = union over e | d of e * primitive cycles of T_(d/e)
    allc = {d: summ[d][0][2] + summ[d][0][3] for d in ds}
    ok = all(allc[d] == sum(om[d // e] for e in range(1, d + 1) if d % e == 0) for d in ds)
    check(ok, "T2.6 Proposition 3 on the sweep: for every d, #cycles(T_d, least <= 1200 d) = sum over e | d of "
              "#primitive cycles(T_(d/e), least <= 1200 d/e)")
    big = {d: om[d] for d in ds if om[d] > 160}
    check(big == BM_TABLE, "T2.6 table (20) reproduced exactly: the eleven d with omega(d) > 160 and their "
                           "omega: " + ", ".join(f"{d}:{big[d]}" for d in sorted(big, key=big.get)))
    check(min(om.values()) >= 1 and dist[1] == BM_N1 and dist[2] == BM_N2,
          f"T2.6 every T_d has a primitive cycle (Lagarias's Conjecture 2(1) in range); omega = 1 for "
          f"{dist[1]} d and omega = 2 for {dist[2]} d, as in Belaga-Mignotte")
    longc = sorted((c for c in cyc if c[4]), key=lambda c: -c[2])
    n600 = sum(1 for c in cyc if c[4] and c[2] > 600)
    only_long = sum(1 for d in ds if om[d] == 1 and any(c[0] == d and c[4] and c[2] > 600 for c in longc[:n600]))
    th = [c[3] / c[2] for c in longc[:n600]]
    check(n600 > 0, f"T2.6 long cycles: {n600} primitive cycles have p > 600 (a/p in [{min(th):.3f}, "
                    f"{max(th):.3f}], mean {sum(th) / len(th):.3f}); for {only_long} d the only primitive cycle is "
                    f"such a long one; the longest: " + ", ".join(f"d={c[0]} least {c[1]} (p,a)=({c[2]},{c[3]})"
                                                                 for c in longc[:4]))
    # EMPIRICAL scaling of the least element of long cycles, and the proved lower bound m > d / 2^(K+1)
    import statistics
    sel = [c for c in cyc if c[4] and c[2] >= 100]
    mpd = [c[1] * c[2] / c[0] for c in sel]
    qs = statistics.quantiles(mpd, n=4)
    lowb = True
    for (d, m, p, a, pr) in sel[:5000]:
        orb = cycle_of(m, 3, d, maxlen=p)
        w = [x & 1 for x in orb]
        K, run = 0, 0
        for b in w + w:
            run = run + 1 if b == 0 else 0
            K = max(K, min(run, p))
        lowb &= m * (1 << (K + 1)) > d
    check(lowb and 1 < qs[1] < 20 and max(c[2] / c[0] for c in sel) < 2,
          f"T2.6 EMPIRICAL: for the {len(sel)} primitive cycles with p >= 100, m p / d has quartiles "
          f"{qs[0]:.2f}, {qs[1]:.2f}, {qs[2]:.2f} (least element m ~ 5 d/p), a/p has median "
          f"{statistics.median(c[3] / c[2] for c in sel):.3f} and p/d <= {max(c[2] / c[0] for c in sel):.2f}; "
          f"PROVED bound m > d/2^(K+1), K = longest run of even steps, checked on 5000 of them")
    check(dist[3] == BM_N3 - 1 and tot == BM_TOTAL - 8,
          f"T2.6 RESIDUAL: omega = 3 for {dist[3]} d (Belaga-Mignotte: {BM_N3}) and {tot} primitive cycles in "
          f"all (Belaga-Mignotte: {BM_TOTAL})")

    # T2.7 extension: least element up to 2e7 (d < 16667) and every near-critical clock with p < 4000
    import mpmath
    mpmath.mp.dps = 40
    Pstar = 4000
    cl = []
    for a in range(1, Pstar):
        pl = (3 ** a).bit_length()
        if pl >= Pstar:
            break
        p = pl
        while p < Pstar and p * 41 < 65 * a:
            cl.append((p, a))
            p += 1
    fac = {c: 1 / (mpmath.power(2, mpmath.mpf(c[0]) / c[1]) - 3) for c in cl}
    rng = os.path.join(SCRATCH, "ranges_ext.txt")
    Xfin = {}
    span = 0
    with open(rng, "w") as fo:
        for d in ds:
            f = mpmath.mpf(1200)
            for (p, a) in cl:
                if (pow(2, p, d) - pow(3, a, d)) % d == 0:
                    f = max(f, fac[(p, a)])
            X = max(int(mpmath.floor(f * d)), 2 * 10 ** 7 if 1200 * d < 2 * 10 ** 7 else 0, 1200 * d)
            Xfin[d] = X
            if X > 1200 * d:
                fo.write(f"{d} {1200 * d + 1} {X}\n")
                span += X - 1200 * d
    outs2 = [os.path.join(SCRATCH, f"sweep_ext_s{s}.txt") for s in (0, 1)]
    t = time.time()
    run_parallel([[SWEEP, "-f", rng, "2", str(s)] for s in (0, 1)], outs2, 2)
    summ2, cyc2, other2 = parse_sweep(outs2)
    t_ext = time.time() - t
    check(not cyc2 and not other2 and all(v[0][4] == 0 and v[0][5] == 0 for v in summ2.values()),
          f"T2.7 extension ({len(summ2)} d, {span:.3e} further integers, i.e. {span / 2:.2e} odd starts, "
          f"{t_ext:.0f} s wall): least elements in "
          f"(1200 d, X_d], X_d = max(1200 d, 2e7 for d < 16667, the perigee windows of all {len(cl)} clocks "
          f"(p,a) with p < {Pstar} and log2 3 < p/a < 65/41 that satisfy d | 2^p - 3^a): NO further cycle")
    # the resulting completeness bound for every d: no lattice clock with p < Pstar has its window reaching X_d + 1
    unc = [d for d in ds if cone_min_clock(d, Xfin[d], pmax=Pstar)[0] is not None]
    check(not unc,
          f"T2.7 hence for every d <= 19999 the {tot} primitive cycles contain every primitive T_d-cycle of "
          f"period p < {Pstar} (no clock of L_d with p < {Pstar} has perigee window reaching X_d + 1); the longest "
          f"cycle found has p = {max(c[2] for c in cyc if c[4])}")
    for d in sorted(BM_TABLE):
        (bp, ba), e = cone_min_clock(d, Xfin[d])
        check(bp >= Pstar, f"T2.7 d = {d}: X_d = {Xfin[d]}, least uncovered clock ({bp},{ba}): omega(d) = "
                           f"{om[d]} counts every primitive cycle of period p < {bp}")

    # T2.8 the residual is two-sided: omega_BM >= omega_ours >= 1 for all d would force
    # {BM = 1} = {ours = 1}, {BM = 2} = {ours = 2} (equal sizes) and {BM = 3} inside {ours = 3}
    check(dist[1] == BM_N1 and dist[2] == BM_N2 and BM_N3 > dist[3],
          f"T2.8 the residual is not explained by Belaga-Mignotte having found extra cycles only: "
          f"omega_BM(d) >= omega_ours(d) >= 1 for all d, with n1 = {BM_N1} and n2 = {BM_N2} equal on both sides, "
          f"would force #{{omega_BM = 3}} <= #{{omega_ours = 3}} = {dist[3]} < {BM_N3}")
    # statistics used in T3
    return cyc, om


# ================================================================================================
# T3
# ================================================================================================
def t3(sweep_cycles):
    header("T3. Stern-Brocot positions (theta_q = log_q 2; a/p < theta_q <=> 2^p > q^a)")
    path3 = sb_path(3, 700)
    path5 = sb_path(5, 700)
    print("      SB path of log_3 2:", " ".join(f"{a}/{p}{s}" for a, p, s in path3[:16]))
    print("      SB path of log_5 2:", " ".join(f"{a}/{p}{s}" for a, p, s in path5[:16]))
    check([(a, p) for a, p, s in path3[:12]] == [(1, 1), (1, 2), (2, 3), (3, 5), (5, 8), (7, 11), (12, 19),
                                                (17, 27), (29, 46), (41, 65), (53, 84), (94, 149)],
          "T3.1 the Stern-Brocot path of log_3 2 begins 1/1 1/2 2/3 3/5 5/8 7/11 12/19 17/27 29/46 41/65 53/84 94/149")
    check([(a, p) for a, p, s in path5[:6]] == [(1, 1), (1, 2), (1, 3), (2, 5), (3, 7), (4, 9)],
          "T3.1 the Stern-Brocot path of log_5 2 begins 1/1 1/2 1/3 2/5 3/7 4/9")

    maps = [("3x+1", 3, 1), ("3x-1", 3, -1), ("5x+1", 5, 1), ("5x-1", 5, -1), ("3x+5", 3, 5),
            ("3x+7", 3, 7), ("3x+13", 3, 13), ("3x+23", 3, 23)]
    print("      map    least    (p,a)    a/p       SB-path  D=2^p-q^a   |D| vs |d|     kind")
    rows = {}
    for name, q, d in maps:
        if q == 3:
            cyc = cycles_on_Z(q, d, 10 ** 6, f"t3_{name}")
        else:
            cyc = census_q(q, d, 60)
        out = []
        for (s, m, p, a) in sorted(set(cyc)):
            D = gap(p, a, q)
            onp, side = on_sb_path(a, p, q)
            g = math.gcd(D, d)
            kind = "free" if d % D == 0 else ("partially free" if g > 1 else "sporadic")
            y = s * m
            f = abs(d) // math.gcd(y, d) if y else 1
            out.append((y, p, a, onp, side, D, kind, f))
            print("      %-6s %7d  (%3d,%3d)  %-8s  %-7s  %10d   %-12s  %s%s" % (
                name, y, p, a, f"{a}/{p}", ("yes-" + side) if onp else "no", D,
                "|D|<=|d|" if abs(D) <= abs(d) else "|D|>|d|", kind,
                "" if f == abs(d) else f" (= {abs(d) // f} x primitive T_(q,{d // (abs(d) // f)}))"))
        rows[name] = out
    r1 = rows["3x+1"]
    check(sorted((y, p, a) for (y, p, a, *_ ) in r1) == [(-17, 11, 7), (-5, 3, 2), (-1, 1, 1), (0, 1, 0), (1, 2, 1)]
          and all(t[3] for t in r1),
          "T3.2 3x+1 on Z (census |least| <= 1e6): the five cycles 0, {1,2}, {-1}, {-5,..}, {-17,..}; all five "
          "densities 0/1, 1/2, 1/1, 2/3, 7/11 lie on the Stern-Brocot path of log_3 2")
    check([t[6] for t in sorted(r1)] == ["sporadic", "free", "free", "free", "free"],
          "T3.2 3x+1: {-17,..} is the only sporadic cycle (|D| = 139); the other four are free (|D| = 1)")
    check(sorted((y, p, a) for (y, p, a, *_ ) in rows["3x-1"]) ==
          [(-1, 2, 1), (0, 1, 0), (1, 1, 1), (5, 3, 2), (17, 11, 7)],
          "T3.2 3x-1 on Z: 0, {1}, {5,7,10}, {17,..}, {-1,-2} -- the negatives of the 3x+1 cycles")
    # the next upper approximant 12/19 carries no integer cycle of 3x+-1
    D19 = gap(19, 12, 3)
    n19 = sum(1 for w in words_of_shape(19, 12) if carry(w, 3) % D19 == 0)
    check(D19 == -7153 and n19 == 0, "T3.3 clock (19,12), 3^12 - 2^19 = 7153: none of the 50388 words has "
                                     "7153 | c_w, so no integer cycle of 3x+-1 sits at the next upper approximant")
    # sweep statistics: best approximations versus gap cofactor
    if sweep_cycles is not None:
        prim = [c for c in sweep_cycles if c[4]]
        on = Counter()
        cof = Counter()
        cache = {}
        for (d, m, p, a, pr) in prim:
            g = math.gcd(a, p)
            key = (a // g, p // g)
            if key not in cache:
                cache[key] = on_sb_path(*key, 3)
            onp, side = cache[key]
            on[onp] += 1
            D = gap(p, a, 3)
            k = (D // d).bit_length() - 1
            cof["D = d (free)" if D == d else ("D/d < 2^10" if k < 10 else ("D/d < 2^100" if k < 100 else
                                                                          "D/d >= 2^100"))] += 1
        check(sum(on.values()) == len(prim),
              f"T3.4 EMPIRICAL (all {len(prim)} primitive cycles of T_d, d <= 19999): density a/p on the "
              f"Stern-Brocot path of log_3 2 for {on[True]} ({100 * on[True] / len(prim):.1f}%), off it for "
              f"{on[False]}; gap cofactor D/d: " + ", ".join(f"{k}: {v}" for k, v in sorted(cof.items())))


def census_q(q, d, pmax):
    """exact per-clock census of the cycles of T_{q,d} on Z with period p <= pmax: for each clock and
    sign, every odd y in the perigee window is iterated p steps.  Returns (sign, least |.|, p, a)."""
    import mpmath
    mpmath.mp.dps = 40
    out = {(0, 0, 1, 0)}
    for sgn, dd in ((1, d), (-1, -d)):
        # positive cycles of T_{q,dd}: y_w = dd c_w / D > 0
        for p in range(1, pmax + 1):
            for a in range(1, p + 1):
                D = gap(p, a, q)
                if dd * D <= 0:
                    continue
                r = mpmath.power(2, mpmath.mpf(p) / a)
                if dd > 0:
                    if r - q <= 0:
                        continue
                    hi = int(mpmath.floor(dd / (r - q)))
                else:
                    if q - r <= 0:
                        continue
                    hi = int(mpmath.floor(-dd / (q - r)))
                lo = max(1, -(-(abs(dd) * q ** (a - 1)) // abs(D)))
                cmax = (1 << (p - a)) * (q ** a - (1 << a)) // (q - 2) if q > 2 else None
                B = abs(dd) * cmax // abs(D) + 1
                for y in range(lo | 1, hi + 1, 2):
                    x = y
                    k = 0
                    ok = True
                    for s in range(1, p + 1):
                        if x & 1:
                            k += 1
                            x = (q * x + dd) >> 1
                        else:
                            x >>= 1
                        if x < y or x > B:
                            ok = False
                            break
                        if x == y and s < p:
                            ok = False
                            break
                    if ok and x == y and k == a:
                        out.add((sgn, y, p, a))
    return sorted(out)


# ================================================================================================
# T4
# ================================================================================================
def t4():
    header("T4. 5x+1 and 3x-1")
    c5 = census_q(5, 1, 60)
    check(c5 == [(-1, 1, 2, 1), (0, 0, 1, 0), (1, 1, 5, 2), (1, 13, 7, 3), (1, 17, 7, 3)],
          "T4.1 5x+1 on Z, exact per-clock census of all periods p <= 60 (both signs): exactly 0, {-1,-2}, "
          "{1,3,8,4,2}, {13,..}, {17,..}")
    info = []
    for (s, m, p, a) in c5:
        D = gap(p, a, 5)
        words = list(words_of_shape(p, a))
        necks = {necklace_rep(w) for w in words if primitive_period(w) == p}
        integral = {necklace_rep(w) for w in words if primitive_period(w) == p and carry(w, 5) % D == 0}
        info.append((s * m, p, a, D, len(integral), len(necks)))
    check(info == [(-1, 2, 1, -1, 1, 1), (0, 1, 0, 1, 1, 1), (1, 5, 2, 7, 1, 2), (13, 7, 3, 3, 2, 5),
                   (17, 7, 3, 3, 2, 5)],
          "T4.2 5x+1 shapes: {-1,-2} on (2,1), D = -1, FREE (1 of 1 necklace); {1,3,8,4,2} on (5,2), D = 7, "
          "SPORADIC (1 of 2); {13,..} and {17,..} on (7,3), D = 3, SPORADIC (2 of 5)")
    check(Fraction(-1, 3) == periodic_point((1,), 5, 1) and
          [s for s in gersonides_solutions(5, 400) if s[1] <= s[0] - 1] == [(2, 1)],
          "T4.3 free cycles of 5x+1 are exactly {0} and {-1,-2}: the only free mixed shape is (2,1) "
          "(5 - 4 = 1; T1.4), and the single-word shape (1,1) gives -1/3")
    # 3x-1 is conjugate to 3x+1 by y -> -y
    ok = all(Tmap(-y, 3, -1) == -Tmap(y, 3, 1) for y in range(-10 ** 5, 10 ** 5))
    check(ok, "T4.4 T_(3,-1)(-y) = -T_(3,1)(y) for |y| < 1e5 (and identically: -y/2, (-3y-1)/2), so the cycles "
              "of 3x-1 are the negatives of those of 3x+1")
    c3p = census_q(3, 1, 90)
    check(c3p == [(-1, 1, 1, 1), (-1, 5, 3, 2), (-1, 17, 11, 7), (0, 0, 1, 0), (1, 1, 2, 1)],
          "T4.5 3x+1 on Z, exact per-clock census of all periods p <= 90 (both signs): 0, {1,2}, {-1}, "
          "{-5,-7,-10}, {-17,..}; hence 3x-1 has 0, {-1,-2}, {1}, {5,7,10}, {17,..} with the same shapes "
          "(1,1) D=-1 free, (3,2) D=-1 free, (11,7) D=-139 sporadic")
    if not QUICK:
        for d, want in ((1, [(1, 2, 1)]), (-1, [(1, 1, 1), (5, 3, 2), (17, 11, 7)])):
            cyc, ent, bad, n = traj_run(3, d, 1, 10 ** 9, f"t4_{d}")
            check(not bad and [(m, p, a) for (m, p, a, mx) in cyc] == want,
                  f"T4.6 least-element search, all odd y <= 1e9: the positive cycles of 3x{'+' if d > 0 else '-'}1 "
                  f"are exactly {want} (least element, p, a)")


def main():
    print("procgen_sporadic_20260926_run.py  (session collatz-procgen-20260922, lane sporadic, 2026-09-26)")
    print("python", sys.version.split()[0], " platform", sys.platform, " quick" if QUICK else "")
    compile_c(os.path.join(HERE, "procgen_sporadic_20260926_traj.c"), TRAJ)
    compile_c(os.path.join(HERE, "procgen_sporadic_20260926_sweep.c"), SWEEP)
    t1()
    print(f"  [T1 done {time.time() - T0:.1f}s]", flush=True)
    res = t2()
    print(f"  [T2 done {time.time() - T0:.1f}s]", flush=True)
    t3(res[0] if res else None)
    print(f"  [T3 done {time.time() - T0:.1f}s]", flush=True)
    t4()
    header("Reproduction")
    for f in ("procgen_sporadic_20260926_lib.py", "procgen_sporadic_20260926_run.py",
              "procgen_sporadic_20260926_traj.c", "procgen_sporadic_20260926_sweep.c"):
        print(f"  sha256 {sha256_file(os.path.join(HERE, f))}  04-computation/experiments/{f}")
    print(f"  wall time {time.time() - T0:.0f} s;  peak RSS: python {peak_rss_mb():.0f} MB, C children "
          f"{child_peak_rss_mb():.1f} MB")
    print(f"  {N_CHECKS[0]} checks passed")


if __name__ == "__main__":
    main()
