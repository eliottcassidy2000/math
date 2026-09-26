#!/usr/bin/env python3
"""
procgen_drift_20260926_run.py -- runner of the drift lane (session collatz-procgen-20260922, 2026-09-26):
is q n + 1 (mainly q = 5) in the Haar closure of the bounded-lookahead-provable sign strategies?

    python3 04-computation/experiments/procgen_drift_20260926_run.py > 05-knowledge/results/procgen_drift_20260926.out

Every claim printed is guarded by check(...), which raises on failure.  Sections:
  A  engine self-tests (C engine vs pure-Python Karp; incremental toggles vs Karp)
  B  Theorem 1 (stationary gain identity) and the gain tail bound
  C  Theorem 2 (entropy / merge inequality), the universal floor p0, emptiness for q >= 23, flip-mass constants
  D  Lemma 3 (carry lemma for 5n+1)
  E  exact distances of 5n+1 (k <= 8) with re-derived no-good certificates
  F  the min-max density rho*(q,k), q in {3,5,7,9,11,13}
  G  upper bounds for 5n+1 up to k = KMAX (lift + prune), stationary diagnostics, fits
  H  controls: q = 3 (exponential rate), q = 7, 9 (necklace bounds, emptiness)
  I  what fails: simple explicit rules
Caches and temporary files live in scratch/procgen_drift/ (not committed).
"""
import os
import sys
import time
import random
import resource
from fractions import Fraction
from math import log, log2, exp, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_drift_20260926_lib as L       # noqa: E402
import procgen_drift_20260926_ihs as I       # noqa: E402
import procgen_drift_20260926_local as LS    # noqa: E402
import procgen_drift_20260926_minmax as MM   # noqa: E402

check = L.check
T0 = time.time()
KMAX = int(os.environ.get('DRIFT_KMAX', '18'))
K8_TIME = float(os.environ.get('DRIFT_K8_TIME', '100000'))
K8_ITERS = int(os.environ.get('DRIFT_K8_ITERS', '10'))
OPT7 = None


def out(*a):
    print(*a, flush=True)


def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == 'darwin' else r / 1024


def section(name):
    out('')
    out('=' * 100)
    out(name)
    out('=' * 100)


# ----------------------------------------------------------------------------- A
def karp_py(k, q, flip):
    """independent pure-Python Karp maximum cycle density (exact Fractions)"""
    N = 1 << k
    H = N >> 1
    t0 = L.targets(k, q, flip)
    NEG = -10 ** 9
    D = [[0] * N]
    for j in range(1, N + 1):
        prev = D[-1]
        cur = [NEG] * N
        for u in range(N):
            if prev[u] > NEG:
                val = prev[u] + (u & 1)
                for v in (int(t0[u]), int(t0[u]) + H):
                    if val > cur[v]:
                        cur[v] = val
        D.append(cur)
    best = None
    for v in range(N):
        if D[N][v] <= NEG:
            continue
        worst = None
        for j in range(N):
            if D[j][v] > NEG:
                f = Fraction(D[N][v] - D[j][v], N - j)
                if worst is None or f < worst:
                    worst = f
        if best is None or worst > best:
            best = worst
    return best


def section_A():
    section("A. Engine self-tests")
    rng = random.Random(1)
    n = 0
    for q in (3, 5, 7, 9):
        for k in (3, 4, 5):
            F = L.best_lower_approx(1 << k, q)
            for _ in range(25):
                flip = np.zeros(1 << k, dtype=np.uint8)
                p = rng.choice([0.2, 0.5, 0.8])
                for r in range(1, 1 << k, 2):
                    flip[r] = 1 if rng.random() < p else 0
                kp = karp_py(k, q, flip)
                E = L.Engine(k, q, flip, F)
                check(E.karp() == kp, "C Karp differs from Python Karp")
                st = L.classify(k, q, flip)
                check((st[0] == 'I') == (q ** kp.numerator < 2 ** kp.denominator), "classify differs from Karp")
                check(L.rho_max(k, q, flip)[0] == kp, "Dinkelbach rho_max differs from Karp")
                n += 1
    out(f"A1. C engine (Karp, potential certificate, Dinkelbach) agrees with an independent pure-Python Karp on {n} "
        f"random strategies (q in {{3,5,7,9}}, k in {{3,4,5}}): OK")
    # incremental toggles from a class-(i) set, against pure-Python Karp
    k, q = 7, 5
    S = set(OPT7)
    F = L.best_lower_approx(1 << k, q)
    E = L.Engine(k, q, L.flip_array(k, S), F)
    check(E.solve(), "k=7 optimum not class (i)")
    acc = rej = 0
    for _ in range(200):
        r = rng.randrange(1, 1 << k, 2)
        before = E.flip().copy()
        res = E.toggle(r)
        after = E.flip()
        psi = E.psi()
        cand = before.copy()
        cand[r] ^= 1
        d = karp_py(k, q, cand)
        check(res == (q ** d.numerator < 2 ** d.denominator), "incremental toggle decision differs from Karp")
        check(np.array_equal(after, cand if res else before), "toggle state wrong")
        check(L.verify_potential(k, q, after, psi, F), "potential invalid after toggle")
        acc += res
        rej += (not res)
    out(f"A2. incremental toggles (q=5, k=7, 200 random toggles from the optimum: {acc} accepted, {rej} rejected) "
        f"agree with pure-Python Karp, and the maintained potential is valid after each: OK")


# ----------------------------------------------------------------------------- B
def section_B():
    section("B. Theorem 1 (stationary gain identity): pi(odd) = theta - sum_{s in R} pi(s) g(s)")
    rng = random.Random(2)
    worst = 0.0
    cnt = 0
    for q in (3, 5, 7, 9):
        for k in (4, 5, 6, 7, 8):
            g = L.gains(k, q)
            for _ in range(6):
                flip = np.zeros(1 << k, dtype=np.uint8)
                p = rng.choice([0.1, 0.3, 0.5])
                for r in range(1, 1 << k, 2):
                    flip[r] = 1 if rng.random() < p else 0
                R = np.nonzero(flip)[0]
                for theta in (0.3, 0.5, 0.8):
                    for nodes, pi in L.stationary_theta(k, q, flip, theta):
                        lhs = pi[1::2].sum()
                        rhs = theta - float((pi[R] * g[R]).sum())
                        worst = max(worst, abs(lhs - rhs))
                        cnt += 1
    check(worst < 1e-9, "gain identity violated")
    out(f"B1. identity checked on {cnt} (strategy, theta, closed class) triples (q in {{3,5,7,9}}, k = 4..8, "
        f"theta in {{0.3,0.5,0.8}}): max |error| = {worst:.1e} (float64 solves): OK")
    # q = 5 exact optimum at k = 7
    k, q = 7, 5
    g = L.gains(k, q)
    fl = L.flip_array(k, OPT7)
    R = np.array(sorted(OPT7))
    for theta in (0.5, 0.9):
        for nodes, pi in L.stationary_theta(k, q, fl, theta):
            po = pi[1::2].sum()
            spg = float((pi[R] * g[R]).sum())
            check(abs(po - (theta - spg)) < 1e-9, "identity on the optimum")
            check(spg > theta - L.log_c(q), "class (i) consequence sum pi g > theta - c violated")
            out(f"B2. k=7 optimum of 5n+1, theta={theta}: pi(odd) = {po:.6f} = theta - sum pi g = {theta - spg:.6f}; "
                f"sum pi g = {spg:.4f} > theta - log_5 2 = {theta - L.log_c(q):.4f}; pi(R) = {pi[R].sum():.4f} "
                f"> (theta - c)/(k-1) = {(theta - L.log_c(q)) / (k - 1):.4f}")
    # gain distribution and the Hoeffding tail bound
    out("B3. gain distribution over odd residues (q = 5) and the tail bound #{g >= x} <= 2^(k-1) * 2 exp(-x^2/(2(k-1))):")
    for k in (8, 12, 16, 20):
        g = L.gains(k, 5)[1::2]
        vals, counts = np.unique(g, return_counts=True)
        tail = np.cumsum(counts[::-1])[::-1]
        for v, t in zip(vals, tail):
            if v > 0:
                check(t <= (1 << (k - 1)) * 2 * exp(-v * v / (2 * (k - 1))) + 1e-9, "gain tail bound violated")
        check(int(g.max()) <= k - 1 and int(g.min()) >= -(k - 1), "|g| <= k-1")
        out(f"    k={k:2d}: mean {g.mean():+.4f}, variance {g.var():.3f} (= {g.var() / k:.3f} k), max {int(g.max())}, "
            f"min {int(g.min())}; tail bound holds for every x > 0")


# ----------------------------------------------------------------------------- C
def section_C(pruned):
    section("C. Theorem 2 (entropy / merge inequality) and its corollaries")
    rng = random.Random(3)
    cnt = 0
    tight = []
    for q in (3, 5, 7, 9, 11):
        for k in (4, 5, 6, 7, 8, 9):
            for _ in range(5):
                flip = np.zeros(1 << k, dtype=np.uint8)
                p = rng.choice([0.1, 0.3, 0.5, 0.7])
                for r in range(1, 1 << k, 2):
                    flip[r] = 1 if rng.random() < p else 0
                R = [int(v) for v in np.nonzero(flip)[0]]
                for nodes, pi in L.stationary(k, q, flip):
                    piv = np.zeros(1 << k)
                    piv[nodes] = pi
                    po = piv[1::2].sum()
                    M, mR, mRs = L.merge_entropy(k, q, R, piv)
                    check(1 - L.h2(po) <= M + 1e-9, "entropy inequality 1 - h(pi(odd)) <= merge entropy violated")
                    check(M <= mR + mRs + 1e-12 and mR + mRs <= po + 1e-12, "merge entropy <= pi(R u R*) <= pi(odd)")
                    check(po + L.h2(po) >= 1 - 1e-9, "universal floor violated")
                    tight.append(po)
                    cnt += 1
    out(f"C1. 1 - h(pi(odd)) <= sum_merge (a+b) h(a/(a+b)) <= pi(R u R*) <= pi(odd) checked on {cnt} closed classes of "
        f"random strategies (q in {{3,5,7,9,11}}, k = 4..9); smallest pi(odd) seen {min(tight):.4f}: OK")
    gl = []
    for q in (5, 7, 9):
        for k in (8, 10, 12):
            gr = [r for r in range(1, 1 << k, 2) if (q * r + 1) % 4 != 0]
            for nodes, pi in L.stationary(k, q, L.flip_array(k, gr)):
                gl.append(float(pi[nodes % 2 == 1].sum()))
    check(all(abs(v - 1 / 3) < 1e-6 for v in gl), "greedy stationary odd frequency 1/3")
    out(f"C1b. the max-halving (greedy) strategy has stationary odd frequency 1/3 (to 1e-6) for q = 5, 7, 9 at k = 8, 10, "
        f"12 ({len(gl)} closed classes)")
    lo, hi = L.p0_root()
    p0 = float(lo)
    out(f"C2. universal floor: p0 = root of p + h(p) = 1 in (0,1/2) lies in [{float(lo):.15f}, {float(hi):.15f}]; every "
        f"stationary law of every strategy has pi(odd) >= p0, hence rho_max >= p0 (sandwich)")
    emp = []
    for q in range(3, 200, 2):
        c = log(2) / log(q)
        if c <= p0:
            emp.append(q)
    check(emp[0] == 23 and all(q in emp for q in range(23, 200, 2)), "emptiness threshold")
    check(log(2) / log(21) > float(hi), "q = 21 is not covered")
    out(f"    => class (i) is EMPTY at every level for every odd q >= 23 (log_q 2 <= log_23 2 = {log(2) / log(23):.5f} "
        f"< p0); q = 21 (log_21 2 = {log(2) / log(21):.5f}) is not covered by this bound")
    out("C3. flip-mass constants for class-(i) strategies (Theorem 2 with pi(odd) < c = log_q 2):")
    for q in (5, 7, 9, 11, 13, 15, 17, 19, 21):
        c = log(2) / log(q)
        B = 1 - L.h2(c)
        a_lo, a_hi = 1e-9, 0.5
        for _ in range(200):
            mid = (a_lo + a_hi) / 2
            if mid * log2(exp(1) / mid) < B:
                a_lo = mid
            else:
                a_hi = mid
        out(f"    q={q:2d}: c = {c:.5f}, drift (1/2)log(q/4) = {0.5 * log(q / 4):+.4f}; pi(R u R*) > 1 - h(c) = {B:.5f}; "
            f"pi(R) > A_q = {a_lo:.5f} (A log2(e/A) = 1 - h(c)); gain form: sum pi g > 1/2 - c = {0.5 - c:.4f}")
    # q = 5 sets: diagnostics
    out("C4. 5n+1 pruned class-(i) sets: entropy inequality and concentration diagnostics (f = dpi/dU):")
    out("      k   Haar(R)  pi(odd)  1-h(pi(odd))  merge-entropy  pi(R)   pi(R*)  U(R u R*)  max f on R u R*  ||f||_2^2")
    for k in sorted(pruned):
        if k > 14:
            continue
        S = pruned[k]
        fl = L.flip_array(k, S)
        N = 1 << k
        for nodes, pi in L.stationary(k, 5, fl):
            piv = np.zeros(N)
            piv[nodes] = pi
            po = piv[1::2].sum()
            M, mR, mRs = L.merge_entropy(k, 5, S, piv)
            Rs = set(L.partner(s, k, 5) for s in S) - set(S)
            A = sorted(set(S) | Rs)
            check(1 - L.h2(po) <= M + 1e-9 and M > 1 - L.h2(L.log_c(5)), "entropy inequality on the pruned set")
            f = piv * N
            check(po < L.log_c(5), "sandwich")
            out(f"     {k:2d}  {len(S) / 2 ** (k - 1):.4f}   {po:.4f}   {1 - L.h2(po):.4f}        {M:.4f}         "
                f"{mR:.4f}  {mRs:.4f}  {len(A) / N:.4f}     {f[A].max():.3f}           {np.sum(f * f) / N:.3f}")


# ----------------------------------------------------------------------------- D
def section_D():
    section("D. Lemma 3 (carry lemma, q = 5): the orbit of y - 1 against the orbit of y under 5n+1")
    T = lambda n: n // 2 if n % 2 == 0 else (5 * n + 1) // 2
    cnt = 0
    for y in range(-200000, 200001):
        a = (y % 4)
        p_y = (y & 1, T(y) & 1)
        p_m = ((y - 1) & 1, T(y - 1) & 1)
        t2y, t2m = T(T(y)), T(T(y - 1))
        if a == 1:
            check(p_y == (1, 1) and p_m == (0, 0) and t2y == 25 * t2m + 8, "case 11")
            check(t2y % 8 == t2m % 8, "25z+8 = z mod 8")
        elif a == 3:
            check(p_y == (1, 0) and p_m == (0, 1) and t2m == t2y - 1, "case 10")
        elif a == 2:
            check(p_y == (0, 1) and p_m == (1, 1) and t2m == 5 * t2y - 7, "case 01")
        else:
            check(p_y == (0, 0) and p_m == (1, 0) and t2m == 5 * t2y - 1, "case 00")
        cnt += 1
    out(f"D1. for all {cnt} integers y in [-2e5, 2e5] (the identities are polynomial, hence 2-adic): "
        f"y=1 mod 4: 11 -> 00, T^2 y = 25 T^2(y-1) + 8 (next 3 parities agree); y=3 mod 4: 10 -> 01, "
        f"T^2(y-1) = T^2 y - 1; y=2 mod 4: 01 -> 11, T^2(y-1) = 5 T^2 y - 7; y=0 mod 4: 00 -> 10, "
        f"T^2(y-1) = 5 T^2 y - 1: OK")
    # consequence: flipping r (odd) = applying y -> y-1 at y = T(r); the carry passes '10' blocks and stops at the first
    # block that is not '10'; gain +2 iff that block is '11'.  Check against the exact gains for k = 12.
    k = 12
    g = L.gains(k, 5)
    W = L.parity_words_all(k, 5)
    agree = 0
    tot = 0
    for s in range(1, 1 << k, 2):
        w = W[s, 1:]          # parities of y = T(s) (length k-1)
        i = 0
        while i + 1 < len(w) and w[i] == 1 and w[i + 1] == 0:
            i += 2
        if i + 1 < len(w):
            blk = (int(w[i]), int(w[i + 1]))
            pred = 2 if blk == (1, 1) else -1
            # the first 2-letter block where the carry stops decides the local gain; later letters can change g
            tot += 1
            agree += (np.sign(g[s]) == np.sign(pred)) or g[s] == 0
    out(f"D2. sign of the exact gain g(s) matches the carry prediction (+ at a '11' stop, - at a '0x' stop) or g = 0 "
        f"for {agree}/{tot} odd residues at k = {k} (EMPIRICAL: later disagreements can reverse the sign)")


# ----------------------------------------------------------------------------- E
def section_E():
    global OPT7
    section("E. Exact distances delta_k(5) of 5n+1 to class (i) (implicit hitting set, RC2 / SAT)")
    for k in range(2, 7):
        S = I.IHS(k, 5, log=lambda *a: None)
        st = S.run(time_limit=600, verbose=False)
        check(st == 'infeasible', f"k={k}: expected empty class (i)")
        for c, m in S.nogoods:
            check(I.nogood_ok(k, 5, c, m), "no-good not genuine")
        from pysat.solvers import Glucose4
        var = {v: i + 1 for i, v in enumerate(S.odd)}
        with Glucose4() as g:
            for lits in S.clauses:
                g.add_clause([var[v] if val == 1 else -var[v] for v, val in lits])
            check(not g.solve(), "no-goods satisfiable")
        out(f"E1. k={k}: class (i) EMPTY: {len(S.clauses)} re-derived expanding-cycle no-goods are UNSAT (Glucose4)")
    k = 7
    S = I.IHS(k, 5, log=lambda *a: None)
    st = S.run(time_limit=600, verbose=False)
    check(st == 'optimal' and S.UB == 29, "k=7 optimum")
    OPT7 = sorted(S.best)
    for c, m in S.nogoods:
        check(I.nogood_ok(k, 5, c, m), "no-good not genuine")
    st2, _ = S.solve_hs(cutoff=28)
    check(st2 == 'cutoff', "a hitting set of size 28 exists")
    rho = L.rho_max(k, 5, L.flip_array(k, OPT7))[0]
    kp = karp_py(k, 5, L.flip_array(k, OPT7))
    check(rho == kp == Fraction(3, 7), "k=7 optimum rho_max")
    out(f"E2. k=7: delta_7(5) = 29 (Haar {29 / 64:.4f}): no hitting set of size 28 over {len(S.clauses)} re-derived "
        f"no-goods (RC2); optimum rho_max = {rho} (Dinkelbach = pure-Python Karp); set {OPT7}")
    # k = 8: RC2 IHS with a time budget (certified lower bound) and the pruned upper bound
    k = 8
    ub = LS.prune(k, 5, LS.lift(set(OPT7), 7), random.Random(8), rounds=3)
    for sd in range(9, 29):
        cand = LS.prune(k, 5, LS.lift(set(OPT7), 7), random.Random(sd), rounds=3)
        if len(cand) < len(ub):
            ub = cand
    S = I.IHS(k, 5, log=lambda *a: None)
    S.offer_ub(ub)
    t = time.time()
    st = S.run(time_limit=K8_TIME, verbose=False, max_iter=K8_ITERS)
    for c, m in S.nogoods:
        check(I.nogood_ok(k, 5, c, m), "no-good not genuine")
    # the lower bound is the exact RC2 optimum of the last hitting-set round (its solution was checked against every
    # no-good); every stored no-good is re-derived above
    lb = min(S.LB, S.UB)
    check(S.verify_class_i(set(S.best)), "UB set not class (i)")
    out(f"E3. k=8: {lb} <= delta_8(5) <= {S.UB} (Haar {lb / 128:.4f} .. {S.UB / 128:.4f}); lower bound = exact RC2 "
        f"optimum of the hitting-set problem over the no-goods of the last round ({len(S.clauses)} no-goods in total, "
        f"each re-derived as a genuine expanding cycle) after {K8_ITERS} IHS rounds ({time.time() - t:.0f}s); upper "
        f"bound = best of 21 seeded prunings of the lifted k=7 optimum, class (i) verified")
    return lb, S.UB


# ----------------------------------------------------------------------------- F
def section_F():
    section("F. Min-max cycle density rho*(q,k) = min over level-k strategies of max cycle density (exact)")
    table = {}
    plan = {3: range(2, 9), 5: range(2, 9), 7: range(2, 10), 9: range(2, 9), 11: range(2, 9), 13: range(2, 9)}
    for q, ks in plan.items():
        row = []
        for k in ks:
            res = MM.minmax_density(k, q, time_limit=900)
            check(res['status'] == 'optimal', f"rho* q={q} k={k} not finished")
            MM.verify_minmax(k, q, res)
            row.append((k, res['rho'], len(res['nogoods'])))
            table[(q, k)] = res['rho']
        c = L.log_c(q)
        out(f"F. q={q:2d} (log_q 2 = {c:.5f}): " + ", ".join(f"k={k}: {r}" for k, r, _ in row))
        nonempty = [k for k, r, _ in row if q ** r.numerator < 2 ** r.denominator]
        out(f"      class (i) nonempty at levels {nonempty if nonempty else 'none'} of those computed; each value "
            f"certified by a strategy with that exact rho_max and a UNSAT set of re-derived no-goods (Glucose4)")
    check(all(table[(q, 7)] == Fraction(3, 7) for q in (5, 7, 9, 11)), "3/7 at k=7")
    check(table[(3, 8)] == Fraction(1, 2) and table[(13, 8)] == Fraction(1, 2), "1/2 values")
    return table


# ----------------------------------------------------------------------------- G
def section_G():
    section(f"G. Upper bounds for 5n+1: lift + prune from the k=7 optimum, k = 7..{KMAX} (every set certified)")
    rng = random.Random(2026)
    S = set(OPT7)
    pruned = {7: sorted(S)}
    rows = []
    for k in range(8, KMAX + 1):
        t = time.time()
        S = LS.lift(S, k - 1)
        best = None
        reps = 3 if k <= 14 else 1
        for _ in range(reps):
            S2 = LS.prune(k, 5, S, rng, rounds=3 if k <= 14 else 1)
            if best is None or len(S2) < len(best):
                best = S2
        S = best
        fl = L.flip_array(k, S)
        st = L.classify(k, 5, fl)
        check(st[0] == 'I', "pruned set not class (i)")
        pruned[k] = sorted(S)
        N5 = L.necklace_lower_bound(k, 5)
        H = 1 << (k - 1)
        rows.append((k, len(S), N5))
        out(f"G1. k={k:2d}: |R| = {len(S):6d}, Haar = {len(S) / H:.5f}, k*Haar = {k * len(S) / H:.3f}; necklace lower "
            f"bound N_k(5) = {N5} (Haar {N5 / H:.4f}, k*Haar {k * N5 / H:.3f}); certificate at F = {st[1]} "
            f"(max psi {int(st[2].max())}) checked on all {2 << k} edges [{time.time() - t:.0f}s]")
    # fits
    ks = np.array([r[0] for r in rows if r[0] >= 10], dtype=float)
    hv = np.array([r[1] / 2 ** (r[0] - 1) for r in rows if r[0] >= 10])
    if len(ks) >= 4:
        def fit(A, y):
            coef = np.linalg.lstsq(A, y, rcond=None)[0]
            return coef, float(np.sqrt(np.mean((A @ coef - y) ** 2)))
        coef, e1 = fit(np.vstack([np.ones_like(ks), ks, np.log2(ks)]).T, np.log2(hv))
        coef3, e3 = fit(np.vstack([np.ones_like(ks), np.log2(ks)]).T, np.log2(hv))
        coef2, e2 = fit(np.vstack([np.ones_like(ks), 1 / ks]).T, hv)
        coef4, e4 = fit(np.vstack([1 / ks]).T, hv)
        out(f"G2. fits over k = 10..{KMAX} (EMPIRICAL; heuristic upper bounds, not delta_k):")
        out(f"    log2 Haar = {coef[0]:+.3f} {coef[1]:+.4f} k {coef[2]:+.3f} log2 k   (rms {e1:.4f})")
        out(f"    log2 Haar = {coef3[0]:+.3f} {coef3[1]:+.3f} log2 k   (power law k^{coef3[1]:.3f}; rms {e3:.4f})")
        out(f"    Haar = {coef2[0]:+.4f} + {coef2[1]:.3f}/k   (rms {e2:.5f});   Haar = {coef4[0]:.3f}/k   (rms {e4:.5f})")
    return pruned


# ----------------------------------------------------------------------------- H
def section_H():
    section("H. Controls")
    # q = 3: the undecided construction of the cube-distance lane (THM-4479): exponential rate, no stationary constraint
    for k in (8, 10, 12, 14):
        W = L.parity_words_all(k, 3)
        bad = []
        for s in range(1, 1 << k, 2):
            a = 0
            ok = True
            for j in range(k):
                a += int(W[s, j])
                if not 3 ** a > 2 ** (j + 1):
                    ok = False
                    break
            if ok:
                bad.append(s)
        st = L.classify(k, 3, L.flip_array(k, bad))
        check(st[0] == 'I', "q=3 undecided flip not class (i)")
        extra = ''
        if k <= 12:
            for nodes, pi in L.stationary(k, 3, L.flip_array(k, bad)):
                piv = np.zeros(1 << k)
                piv[nodes] = pi
                po = piv[1::2].sum()
                M, mR, mRs = L.merge_entropy(k, 3, bad, piv)
                check(po < 0.5, "q=3 flips remove odd steps")
                extra = (f"; pi(odd) = {po:.4f} < 1/2 < log_3 2, pi(R) = {mR:.4f}, merge entropy {M:.4f} "
                         f"(for q = 3 Theorems 1-2 demand nothing: pi(odd) = 1/2 is already below log_3 2)")
        out(f"H1. q=3, k={k}: flipping Bad_k (|Bad_k| = {len(bad)}, Haar {len(bad) / 2 ** (k - 1):.4f}) is class (i) "
            f"(certificate at {st[1]}){extra}")
    for q in (7, 9):
        for k in (8, 10, 12, 14, 16):
            N = L.necklace_lower_bound(k, q)
            out(f"H2. q={q}, k={k}: expanding necklaces N_k({q}) = {N} (Haar {N / 2 ** (k - 1):.4f}, "
                f"k*Haar {k * N / 2 ** (k - 1):.3f}) -- a lower bound for the distance if class (i) were nonempty")


# ----------------------------------------------------------------------------- I
def section_I():
    section("I. What fails: simple explicit flip rules for 5n+1 (each leaves an expanding cycle)")
    q = 5
    for k in (8, 10, 12, 14):
        N = 1 << k
        W = L.parity_words_all(k, q)
        odd = list(range(1, N, 2))
        hts = []
        for s in odd:
            h = 0.0
            m = 1e9
            for x in W[s]:
                h += log(2.5) if x else -log(2)
                m = min(m, h)
            hts.append(m)
        bad = [s for s, m in zip(odd, hts) if m > 0]
        idx = {tuple(int(x) for x in W[s]): s for s in range(N)}
        seen = set()
        low, hi = [], []
        for s in odd:
            w = tuple(int(x) for x in W[s])
            rots = [w[i:] + w[:i] for i in range(k)]
            key = min(rots)
            if key in seen:
                continue
            seen.add(key)
            if not q ** sum(w) > 2 ** k:
                continue
            cand = []
            for r in rots:
                if r[0] != 1:
                    continue
                h = 0.0
                ok = True
                for x in r:
                    h += log(2.5) if x else -log(2)
                    if h <= 0:
                        ok = False
                        break
                if ok:
                    cand.append(r)
            check(len(cand) > 0, "cycle lemma: an expanding necklace has a ballot rotation")
            low.append(idx[min(cand)])
            hi.append(idx[max(cand)])
        rules = [('all 111-windows (r = 5 mod 8)', [s for s in odd if s % 8 == 5]),
                 ('max-halving (r = 1 mod 4)', [s for s in odd if s % 4 == 1]),
                 ('undecided Bad_k (Theorem-1 analogue)', bad),
                 ('Bad_k and 111', [s for s in bad if s % 8 == 5]),
                 ('one flip per expanding necklace (lowest ballot rotation)', low),
                 ('one flip per expanding necklace (highest ballot rotation)', hi)]
        g = L.gains(k, q)
        rules.append(('all positive-gain residues (g >= 1)', [s for s in odd if g[s] >= 1]))
        for name, R in rules:
            st = L.classify(k, q, L.flip_array(k, R))
            check(st[0] == 'X', f"rule unexpectedly class (i): {name}")
            a, p = st[2]
            out(f"I. k={k:2d} {name:58s} Haar {len(R) / 2 ** (k - 1):.3f}: expanding cycle of length {p} with {a} odd "
                f"nodes (5^{a} > 2^{p}), verified walk")


# ----------------------------------------------------------------------------- J
def section_J(pruned):
    section("J. Structure of the certified 5n+1 sets (critical flips, prefix classes, stationary visits)")
    q = 5
    for k in (10, 12, 14):
        S = set(pruned[k])
        F = L.best_lower_approx(1 << k, q)
        lens = []
        ncrit = 0
        for f in sorted(S):
            st = L.classify(k, q, L.flip_array(k, S - {f}), F)
            if st[0] == 'X':
                ncrit += 1
                lens.append(st[2][1])
        rm, cyc, _ = L.rho_max(k, q, L.flip_array(k, S))
        out(f"J1. k={k}: {ncrit}/{len(S)} flips are critical (removing one leaves a verified expanding cycle; lengths "
            f"median {int(np.median(lens))}, min {min(lens)}, max {max(lens)}); rho_max = {rm}, attained by a verified "
            f"cycle of length {len(cyc)}")
        N = 1 << k
        W = L.parity_words_all(k, q)
        cls = {5: '111', 1: '110', 7: '101', 3: '100'}
        frac = {m: sum(1 for s in S if s % 8 == m) / (N // 8) for m in cls}
        (nodes, pi), = L.stationary(k, q, L.flip_array(k, S))
        f = np.zeros(N)
        f[nodes] = pi * N
        vis = {}
        for m in cls:
            a = [f[s] for s in range(m, N, 8) if s in S]
            b = [f[s] for s in range(m, N, 8) if s not in S]
            vis[m] = (float(np.mean(a)) if a else float('nan'), float(np.mean(b)) if b else float('nan'))
        out("    flipped fraction of each prefix class: " + ", ".join(f"{cls[m]}: {frac[m]:.3f}" for m in (5, 1, 7, 3)))
        out("    mean dpi/dU (flipped / unflipped) per class: " +
            ", ".join(f"{cls[m]}: {vis[m][0]:.2f} / {vis[m][1]:.2f}" for m in (5, 1, 7)))
        if k == 14:
            w4 = {}
            for s in S:
                key = ''.join(str(int(x)) for x in W[s, :4])
                w4[key] = w4.get(key, 0) + 1
            tot4 = {}
            for s in range(1, N, 2):
                key = ''.join(str(int(x)) for x in W[s, :4])
                tot4[key] = tot4.get(key, 0) + 1
            out("    flipped fraction by 4-letter prefix: " + ", ".join(f"{kk}: {w4.get(kk, 0) / tot4[kk]:.2f}"
                                                                    for kk in sorted(tot4)))


def main():
    out("procgen_drift_20260926 -- drift lane: distance of q n + 1 (q = 5 mainly) to bounded-lookahead provability")
    out(f"python {sys.version.split()[0]}, numpy {np.__version__}; KMAX = {KMAX}, K8_TIME = {K8_TIME:.0f}s")
    section_E()
    section_A()
    section_B()
    section_D()
    table = section_F()
    pruned = section_G()
    section_C(pruned)
    section_H()
    section_I()
    section_J(pruned)
    section("Z. Resources")
    out(f"wall time {time.time() - T0:.0f}s; peak RSS of the runner {rss_mb():.0f} MB")
    for f in ('procgen_drift_20260926_engine.c', 'procgen_drift_20260926_lib.py', 'procgen_drift_20260926_ihs.py',
              'procgen_drift_20260926_local.py', 'procgen_drift_20260926_minmax.py', 'procgen_drift_20260926_run.py'):
        out(f"sha256 {f} {L.sha256(os.path.join(HERE, f))}")
    out("ALL CHECKS PASSED")


if __name__ == '__main__':
    main()
