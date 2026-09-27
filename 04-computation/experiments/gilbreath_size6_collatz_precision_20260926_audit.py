#!/usr/bin/env python3
"""gilbreath_size6_collatz_precision_20260926_audit.py -- independent auditor verification, written from scratch, of
 (1) 05-knowledge/results/gilbreath_size6_extinction_20260926.md + HYP-9163 + gilbreath_size6_exact_20260926.py/.out
 (2) 05-knowledge/results/collatz_precision_residual_20260926.md + THM-4512 + collatz_precision_residual_20260926.py/.out
     + collatz_coefficient_stopping_20260926.py/.out
 (3) the wiring (INDEX, PROBLEM-LEDGER, synthesis 2s) and the swaplift bank figures.

Nothing here reuses the audited code. Methods:
 G1  Gilbreath lone-defect chain rebuilt (BFS); floats by value iteration AND by a sparse direct solve (scipy);
     exact values verified by solving the linear system modulo two large primes (a claimed rational a/b is confirmed iff
     the modular mean equals a*b^-1 mod p for both primes); exact Fractions recomputed for small F by my own elimination.
 G2  Theorem A tested against the FULL automaton: all 2^(F-1+R) (left word, right sea) configurations are iterated with
     the actual Gilbreath rule for R rows; the exact count of configurations whose leading 1 is destroyed within R rows
     must equal the chain's finite-horizon absorption count (integer weights). Also the F = 3 hand computation.
 G3  per-z decomposition, closed forms, shift law, halving, ratios, growth claims (1.4 log_2 F, 0.7 log_2 F / F).
 C1  U^j(n) = (3^j n + S_j)/2^A and S_(j+1) = 3 S_j + 2^(A_j); the bound S_j <= 2^A ((3/2)^j - 1); Terras class structure
     (mod 2^A versus mod 2^(A+1)); representative formula.
 C2  2^A - 3^j > (3/2)^j - 1 at the minimal A for j <= 5000 (exact integers), worst ratio.
 C3  first-descent classes with j <= 14: rigorous exception search (stop when N(w) < 1 is forced) + replication of the
     audited enumeration's counting rule.
 C4  sigma versus sigma_inf, Syracuse coding, all odd 3 <= n <= NMAX (own vectorised code); max stopping time; PLUS the
     T-coded (Terras) tau_T versus sigma_T comparison on the same range.
 C5  D(k) by an exact integer DP (no floating barrier tests), W_k, rates to large k; the claimed limit 0.9659 versus 0.9466.
 C6  thresholds N(w)/2^A per (j, A); the swaplift bank (density, J+1 <= 41, K <= 65, containment in C_41);
     window-avoidance density (2^-L claim) exactly.
Usage: python3 gilbreath_size6_collatz_precision_20260926_audit.py [--quick]
"""
import sys, time, math, os, re
from fractions import Fraction
from math import comb
import numpy as np

QUICK = '--quick' in sys.argv
T0 = time.time()
RESULTS = []
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, (" -- " + detail) if detail else ""), flush=True)


def note(msg):
    print("  " + msg, flush=True)


def elapsed():
    return time.time() - T0


P1 = (1 << 61) - 1            # Mersenne prime
P2 = (1 << 62) - 57           # prime (largest below 2^62); primality re-checked below


def is_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d = n - 1; r = 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ----------------------------------------------------------------------------------------------------------------------
# G. GILBREATH
# ----------------------------------------------------------------------------------------------------------------------

def g_step(state, u):
    F = len(state)
    return tuple([abs(state[i] - state[i + 1]) for i in range(F - 1)] + [abs(state[F - 1] - u)])


def g_build(d, F):
    init = [tuple(2 * ((b >> i) & 1) for i in range(F - 1)) + (d,) for b in range(1 << (F - 1))]
    index = {}; states = []
    for s in init:
        if s not in index:
            index[s] = len(states); states.append(s)
    i = 0
    while i < len(states):
        s = states[i]; i += 1
        if s[0] >= 4 or max(s) <= 2:
            continue
        for u in (0, 2):
            t = g_step(s, u)
            if t not in index:
                index[t] = len(states); states.append(t)
    n = len(states)
    n0 = np.zeros(n, dtype=np.int64); n1 = np.zeros(n, dtype=np.int64); kind = np.full(n, -1, dtype=np.int64)
    for i, s in enumerate(states):
        if s[0] >= 4:
            kind[i] = 1
        elif max(s) <= 2:
            kind[i] = 0
        else:
            n0[i] = index[g_step(s, 0)]; n1[i] = index[g_step(s, 2)]
    return states, index, n0, n1, kind, [index[s] for s in init]


def g_iterate(n0, n1, kind, target, iters=200000, tol=1e-15):
    """value iteration for P(hit target-kind absorbing set); target = 1 (extinction) or 'any' (absorbed at all)"""
    trans = kind < 0
    if target == 'any':
        P = (kind >= 0).astype(np.float64)
    else:
        P = (kind == target).astype(np.float64)
    for it in range(iters):
        Pn = np.where(trans, 0.5 * (P[n0] + P[n1]), P)
        if np.max(np.abs(Pn - P)) < tol:
            return Pn, it + 1
        P = Pn
    return P, iters


def g_direct(n0, n1, kind):
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl
    trans = np.nonzero(kind < 0)[0]
    m = len(trans)
    col = -np.ones(len(kind), dtype=np.int64); col[trans] = np.arange(m)
    rhs = np.zeros(m)
    rows = [np.arange(m)]; cols = [np.arange(m)]; vals = [np.ones(m)]
    for nb in (n0, n1):
        succ = nb[trans]
        rhs += 0.5 * (kind[succ] == 1)
        mask = kind[succ] < 0
        rows.append(np.nonzero(mask)[0]); cols.append(col[succ[mask]]); vals.append(-0.5 * np.ones(int(mask.sum())))
    M = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(m, m))
    x = spl.spsolve(M, rhs)
    P = (kind == 1).astype(np.float64); P[trans] = x
    return P


def g_solve_mod(n0, n1, kind, p, order=None):
    """exact solution mod prime p of x_i = (x_n0 + x_n1)/2 (transient), x = 1 / 0 on the absorbing kinds"""
    n = len(kind)
    inv2 = pow(2, -1, p)
    trans = [i for i in range(n) if kind[i] < 0]
    if order == 'rev':
        trans = trans[::-1]
    pivots = {}
    porder = []
    for i in trans:
        row = {i: 1}; rhs = 0
        for j in (int(n0[i]), int(n1[i])):
            if kind[j] == 1:
                rhs = (rhs + inv2) % p
            elif kind[j] == 0:
                pass
            else:
                row[j] = (row.get(j, 0) - inv2) % p
        while True:
            ks = [k for k in row if k in pivots]
            if not ks:
                break
            for k in ks:
                c = row.pop(k)
                if c == 0:
                    continue
                prow, prhs = pivots[k]
                for kk, v in prow.items():
                    if kk == k:
                        continue
                    row[kk] = (row.get(kk, 0) - c * v) % p
                rhs = (rhs - c * prhs) % p
            row = {kk: v for kk, v in row.items() if v != 0}
        if not row:
            assert rhs == 0, "inconsistent system"
            raise RuntimeError("singular system (closed transient class?)")
        k = i if i in row else min(row)
        c = pow(row[k], -1, p)
        row = {kk: (v * c) % p for kk, v in row.items()}
        rhs = (rhs * c) % p
        pivots[k] = (row, rhs); porder.append(k)
    sol = {}
    for k in reversed(porder):
        row, rhs = pivots[k]
        val = rhs
        for kk, v in row.items():
            if kk != k:
                val = (val - v * sol[kk]) % p
        sol[k] = val
    x = {}
    for i in range(n):
        if kind[i] == 1:
            x[i] = 1
        elif kind[i] == 0:
            x[i] = 0
        else:
            x[i] = sol[i]
    return x


def g_solve_frac(n0, n1, kind):
    """exact rationals by my own elimination (BFS order, own-variable pivot)"""
    n = len(kind)
    half = Fraction(1, 2)
    trans = [i for i in range(n) if kind[i] < 0]
    pivots = {}; porder = []
    for i in trans:
        row = {i: Fraction(1)}; rhs = Fraction(0)
        for j in (int(n0[i]), int(n1[i])):
            if kind[j] == 1:
                rhs += half
            elif kind[j] == 0:
                pass
            else:
                row[j] = row.get(j, Fraction(0)) - half
        while True:
            ks = [k for k in row if k in pivots]
            if not ks:
                break
            for k in ks:
                c = row.pop(k)
                if c == 0:
                    continue
                prow, prhs = pivots[k]
                for kk, v in prow.items():
                    if kk == k:
                        continue
                    row[kk] = row.get(kk, Fraction(0)) - c * v
                rhs -= c * prhs
            row = {kk: v for kk, v in row.items() if v != 0}
        if not row:
            raise RuntimeError("singular")
        k = i if i in row else min(row)
        c = 1 / row[k]
        row = {kk: v * c for kk, v in row.items()}; rhs = rhs * c
        pivots[k] = (row, rhs); porder.append(k)
    sol = {}
    for k in reversed(porder):
        row, rhs = pivots[k]
        val = rhs
        for kk, v in row.items():
            if kk != k:
                val -= v * sol[kk]
        sol[k] = val
    x = {}
    for i in range(n):
        x[i] = Fraction(1) if kind[i] == 1 else (Fraction(0) if kind[i] == 0 else sol[i])
    return x


def front_path_weight(bits, F):
    """my derivation: the first front (speed 1) at row s+1 sits at column F-1-s and meets the sea cell (row s, column F-1-s)
    = XOR_(j subset s) b(F-1-s+j) (sea kernel), s = 0..F-2; bits[c-1] = bit of column c."""
    w = 0
    for s in range(F - 1):
        c = F - 1 - s
        v = 0
        for j in range(s + 1):
            if (j & s) == j:
                v ^= bits[c - 1 + j]
        w += v
    return w


def zeros_before_wall(bits, F):
    z = 0
    while z < F - 1 and bits[F - 2 - z] == 0:
        z += 1
    return z


def front_only(d, F):
    j = d // 2
    return Fraction(sum(comb(F - 1, t) for t in range(0, j - 1)), 2 ** (F - 1))


def modfrac(fr, p):
    return (fr.numerator % p) * pow(fr.denominator % p, -1, p) % p


CLAIMED_P6 = {3: Fraction(13, 16), 4: Fraction(17, 32), 5: Fraction(349, 1024), 6: Fraction(413, 2048),
              7: Fraction(493, 4096), 8: Fraction(557, 8192), 9: Fraction(159469, 4194304),
              10: Fraction(175853, 8388608), 11: Fraction(196333, 16777216)}
CLAIMED_P8 = {3: Fraction(1), 4: Fraction(29, 32), 5: Fraction(703, 960), 6: Fraction(16243, 30720),
              7: Fraction(388363, 1044480), 8: Fraction(509743, 2088960), 9: Fraction(2576519, 16711680),
              10: Fraction(202620131, 2139095040), 11: Fraction(64449019483, 1099494850560)}
CLAIMED_C6 = {3: Fraction(1, 2), 4: Fraction(1, 2), 5: Fraction(29, 32), 6: Fraction(29, 32), 7: Fraction(45, 32),
              8: Fraction(45, 32), 9: Fraction(12013, 8192), 10: Fraction(12013, 8192), 11: Fraction(16109, 8192)}
CLAIMED_C6_FLOAT = {12: 1.9664, 13: 2.3727, 14: 2.3727, 15: 2.8727, 16: 2.8727, 17: 2.8732}
CLAIMED_STATES6 = [30, 97, 253, 669, 1639, 3952, 8795, 20084, 44559, 101200, 220506, 484666, 1047839, 2271338, 4766630]
CLAIMED_RATIOS6 = [1.083, 1.063, 1.091, 1.076, 1.100, 1.088, 1.081, 1.073, 1.089, 1.082, 1.091, 1.085, 1.096, 1.090, 1.085]
# per-z claimed table (2^F * excess by z), exact for F <= 11
CZ = {}
for F in (3, 4):
    CZ[F] = {0: Fraction(1, 2)}
for F in (5, 6):
    CZ[F] = {0: Fraction(21, 32), 1: Fraction(1, 8), 2: Fraction(1, 8)}
for F in (7, 8):
    CZ[F] = {0: Fraction(21, 32), 1: Fraction(1, 8), 2: Fraction(1, 8), 4: Fraction(1, 2)}
for F in (9, 10):
    CZ[F] = {0: Fraction(5461, 8192), 1: Fraction(273, 2048), 2: Fraction(273, 2048), 3: Fraction(1, 128),
             4: Fraction(261, 512), 5: Fraction(1, 128), 6: Fraction(1, 128)}
CZ[11] = dict(CZ[9]); CZ[11][8] = Fraction(1, 2)


def gilbreath_part():
    print("=" * 110)
    print("G. GILBREATH SIZE-6 NOTE / HYP-9163 / gilbreath_size6_exact_20260926.py")
    print("=" * 110)
    check("P1, P2 are primes (Miller-Rabin, deterministic bases)", is_probable_prime(P1) and is_probable_prime(P2))
    Fmax = {4: 12, 6: (13 if QUICK else 15), 8: (10 if QUICK else 12)}
    Fexact = {6: (9 if QUICK else 11), 8: (8 if QUICK else 11)}
    have_scipy = True
    try:
        import scipy  # noqa
    except Exception:
        have_scipy = False
    data = {}
    for d in (4, 6, 8):
        print("-- d = %d --" % d)
        for F in range(3, Fmax[d] + 1):
            t0 = time.time()
            states, index, n0, n1, kind, init_idx = g_build(d, F)
            nst = len(states)
            P, its = g_iterate(n0, n1, kind, 1)
            Q, _ = g_iterate(n0, n1, kind, 'any')
            pf = float(np.mean(P[init_idx]))
            fo = front_only(d, F)
            ratio = pf / float(fo)
            line = " F=%2d states=%8d p=%.12f front-only=%.12f ratio=%.6f 2^F*excess=%.6f iters=%d absorb-min=%.3e (%.1fs)" % (
                F, nst, pf, float(fo), ratio, (pf - float(fo)) * 2 ** F, its, 1 - Q.min(), time.time() - t0)
            print(line, flush=True)
            rec = {'states': nst, 'p': pf, 'fo': fo, 'P': P, 'init': init_idx, 'Qmin': float(Q.min())}
            if have_scipy and nst <= 250000:
                Pd = g_direct(n0, n1, kind)
                rec['pdirect'] = float(np.mean(Pd[init_idx]))
                if abs(rec['pdirect'] - pf) > 1e-10:
                    check("d=%d F=%d value iteration vs sparse direct solve agree" % (d, F), False, "%.15g vs %.15g" % (pf, rec['pdirect']))
            if d in (6, 8) and F <= Fexact[d]:
                t1 = time.time()
                x1 = g_solve_mod(n0, n1, kind, P1)
                x2 = g_solve_mod(n0, n1, kind, P2)
                inv1 = pow(1 << (F - 1), -1, P1); inv2 = pow(1 << (F - 1), -1, P2)
                m1 = sum(x1[i] for i in init_idx) * inv1 % P1
                m2 = sum(x2[i] for i in init_idx) * inv2 % P2
                rec['mod'] = (x1, x2, m1, m2)
                claimed = (CLAIMED_P6 if d == 6 else CLAIMED_P8).get(F)
                ok = claimed is not None and modfrac(claimed, P1) == m1 and modfrac(claimed, P2) == m2
                check("d=%d F=%d claimed exact p = %s confirmed modulo two primes (%.1fs)" % (d, F, claimed, time.time() - t1), ok,
                      "" if ok else "modular means %d, %d" % (m1, m2))
                if F <= 8:
                    xf = g_solve_frac(n0, n1, kind)
                    pe = sum(xf[i] for i in init_idx) / (1 << (F - 1))
                    check("d=%d F=%d own exact Fraction solve = %s equals claimed" % (d, F, pe), pe == claimed)
                    rec['frac'] = xf
                # excess and its 2^F multiple
                if d == 6:
                    ex = claimed - fo
                    check("d=6 F=%d 2^F * (p - F 2^(1-F)) = %s equals claimed c_6(F) = %s" % (F, ex * 2 ** F, CLAIMED_C6[F]), ex * 2 ** F == CLAIMED_C6[F])
            data[(d, F)] = rec
    # ---- state counts, ratios, float c_6 ----
    got = [data[(6, F)]['states'] for F in range(3, Fmax[6] + 1)]
    check("d=6 reachable state counts F=3..%d match the note's list" % Fmax[6], got == CLAIMED_STATES6[:len(got)], "%s" % got)
    got_r = [data[(6, F)]['p'] / float(data[(6, F)]['fo']) for F in range(3, Fmax[6] + 1)]
    check("d=6 ratios p_6/(F 2^(1-F)) F=3..%d match the note's list to +-0.0006 (1.0625 is printed 1.063 there)" % Fmax[6],
          all(abs(a - b) <= 6e-4 for a, b in zip(got_r, CLAIMED_RATIOS6)), "%s" % ["%.4f" % r for r in got_r])
    for F, cl in CLAIMED_C6_FLOAT.items():
        if (6, F) in data:
            c6 = (data[(6, F)]['p'] - float(data[(6, F)]['fo'])) * 2 ** F
            check("d=6 F=%d floating c_6(F) = %.6f matches claimed %.4f" % (F, c6, cl), abs(c6 - cl) < 6e-5)
    # p_4 = 2^(1-F)
    ok4 = all(abs(data[(4, F)]['p'] - 2.0 ** (1 - F)) < 1e-12 for F in range(3, Fmax[4] + 1))
    check("d=4: p_4(F) = 2^(1-F) to 1e-12 for F = 3..%d (THM-4511 consistency)" % Fmax[4], ok4)
    # absorption certainty
    okQ = all(data[k]['Qmin'] > 1 - 1e-12 for k in data)
    check("every transient state is absorbed with probability 1 (no closed transient class; linear system nonsingular)", okQ)
    # size-8 ratios range and Fermat denominators
    r8 = [data[(8, F)]['p'] / float(data[(8, F)]['fo']) for F in range(4, Fmax[8] + 1)]
    check("d=8 excess ratios lie in [1.036, 1.082] for F=4..%d (note: 1.036-1.082 for F=4..17)" % Fmax[8], min(r8) >= 1.035 and max(r8) <= 1.0825, "min %.4f max %.4f" % (min(r8), max(r8)))
    fermat = {3: 1, 15: 3 * 5, 255: 3 * 5 * 17, 65535: 3 * 5 * 17 * 257}
    dens = []
    for F, fr in CLAIMED_P8.items():
        den = fr.denominator; m = 0
        while den % 2 == 0:
            den //= 2; m += 1
        dens.append((F, m, den))
    check("d=8 claimed denominators are 2^m * {1, 15, 255, 65535} (Fermat products)", all(o in (1, 15, 255, 65535) for _, _, o in dens), "%s" % dens)
    # S10 comparison: F=8: 0.067993 ; F=12: 0.006338 (finite-context R=9)
    if (6, 12) in data:
        p12 = data[(6, 12)]['p']
        note("S10 finite-context values: F=8 quoted 0.067993 (chain %.6f); F=12 quoted 0.006338 (chain %.7f: differs in the 4th significant digit; S10 used R = 9 there)" % (data[(6, 8)]['p'], p12))
        check("S10's F=12 finite-context value 0.006338 agrees with the chain to 'the digits it had'", abs(p12 - 0.006338) < 5e-7, "chain %.7f vs 0.006338: |diff| = %.1e" % (p12, abs(p12 - 0.006338)))
    # ---- G2: exhaustive full-automaton test of Theorem A ----
    print("-- G2: Theorem A against the full automaton (exact finite-horizon counts) --")
    R = 14
    for d, Fs in ((6, (3, 4, 5, 6, 7)), (8, (3, 4, 5, 6)), (4, (3, 4, 5, 6))):
        for F in Fs:
            states, index, n0, n1, kind, init_idx = g_build(d, F)
            # chain finite horizon with integer weights
            w = {}
            for i in init_idx:
                w[i] = w.get(i, 0) + 1
            ext_chain = 0
            for s in range(R):
                nw = {}
                for i, c in w.items():
                    if kind[i] == 1:
                        ext_chain += c << (R - s)
                    elif kind[i] == 0:
                        continue
                    else:
                        a = int(n0[i]); b = int(n1[i])
                        nw[a] = nw.get(a, 0) + c; nw[b] = nw.get(b, 0) + c
                w = nw
            # full automaton
            nL = F - 1; total = 1 << (nL + R)
            idx = np.arange(total, dtype=np.int64)
            width = F + R + 1
            row = np.zeros((total, width), dtype=np.int8)
            row[:, 0] = 1
            for i in range(nL):
                row[:, 1 + i] = 2 * ((idx >> i) & 1)
            row[:, F] = d
            for r in range(R):
                row[:, F + 1 + r] = 2 * ((idx >> (nL + r)) & 1)
            destroyed = np.zeros(total, dtype=bool)
            lead_bad = np.zeros(total, dtype=bool)
            firstrow = np.full(total, -1, dtype=np.int64)
            for s in range(R):
                hit = (row[:, 1] >= 4) & ~destroyed
                firstrow[hit] = s
                destroyed |= hit
                row = np.abs(row[:, :-1] - row[:, 1:])
                lead_bad |= (row[:, 0] != 1)
            ext_full = int(destroyed.sum())
            check("d=%d F=%d R=%d: full-automaton destruction count %d == chain finite-horizon count %d (of 2^%d)" % (d, F, R, ext_full, ext_chain, nL + R), ext_full == ext_chain)
            check("d=%d F=%d: 'a(1) >= 4 at row s' == 'leading element != 1 at row s+1' (event equivalence)" % (d, F), bool(np.array_equal(destroyed, lead_bad)))
            if d == 6 and F == 3:
                # word (columns 1,2) = (0,2): bit0 = 0, bit1 = 1
                sel = ((idx & 1) == 0) & (((idx >> 1) & 1) == 1)
                r0 = (idx >> 2) & 1; r1 = (idx >> 3) & 1
                pred = (r0 == 0) & (r1 == 0)
                check("F=3 by hand: word (0,2) destroyed iff the first two right bits are 0 (exhaustive, R=14)", bool(np.array_equal(destroyed[sel], pred[sel])))
                # the other three words are destroyed by the first front regardless of the sea
                check("F=3: the three weight<=1 words are destroyed for every sea", bool(destroyed[~sel].all()))
    # ---- G3: per-z decomposition, closed forms, laws ----
    print("-- G3: per-z anatomy, closed forms, shift law, halving, growth --")
    byz_float = {}
    for F in range(3, Fmax[6] + 1):
        rec = data[(6, F)]
        states, index, n0, n1, kind, init_idx = g_build(6, F)   # rebuild (cheap for these F) to get word->state map
        P = rec['P']
        bz = {}; bw = {}
        for b in range(1 << (F - 1)):
            bits = [(b >> i) & 1 for i in range(F - 1)]
            s = tuple(2 * x for x in bits) + (6,)
            z = zeros_before_wall(bits, F); wgt = front_path_weight(bits, F)
            pr = P[index[s]]
            bz[z] = bz.get(z, 0.0) + (pr - (1.0 if wgt <= 1 else 0.0)) * 2.0    # 2^F * excess_z = 2^F * sum/2^(F-1) = 2 * sum
        byz_float[F] = bz
        # front-only words: exactly F words have weight <= 1, and each is destroyed with probability 1
        cnt = sum(1 for b in range(1 << (F - 1)) if front_path_weight([(b >> i) & 1 for i in range(F - 1)], F) <= 1)
        surely = all(abs(P[index[tuple(2 * ((b >> i) & 1) for i in range(F - 1)) + (6,)]] - 1) < 1e-12
                     for b in range(1 << (F - 1)) if front_path_weight([(b >> i) & 1 for i in range(F - 1)], F) <= 1)
        if F <= 12:
            check("F=%d: exactly F=%d left words have first-front path weight <= 1, all destroyed with probability 1" % (F, F), cnt == F and surely)
        # exact per-z modulo primes for F <= Fexact
        if 'mod' in rec:
            x1, x2, _, _ = rec['mod']
            okz = True; bad = []
            for z in range(0, F - 2 + 1):
                s1 = 0; s2 = 0
                for b in range(1 << (F - 1)):
                    bits = [(b >> i) & 1 for i in range(F - 1)]
                    if zeros_before_wall(bits, F) != z:
                        continue
                    st = index[tuple(2 * x for x in bits) + (6,)]
                    fo_ind = 1 if front_path_weight(bits, F) <= 1 else 0
                    s1 = (s1 + x1[st] - fo_ind) % P1; s2 = (s2 + x2[st] - fo_ind) % P2
                v1 = s1 * 2 % P1; v2 = s2 * 2 % P2
                cl = CZ.get(F, {}).get(z, Fraction(0))
                if not (modfrac(cl, P1) == v1 and modfrac(cl, P2) == v2):
                    okz = False; bad.append(z)
            check("F=%d exact per-z table (2^F * excess by z) confirmed modulo two primes for all z" % F, okz, "" if okz else "mismatch at z=%s" % bad)
    # closed forms
    ok = True
    for F, cl in CZ.items():
        k = 1
        while not (2 ** k < F <= 2 ** (k + 1)):
            k += 1
        c0 = Fraction(2, 3) * (1 - Fraction(1, 4 ** (2 ** k - 1)))
        if cl.get(0) != c0:
            ok = False
        if F >= 5:
            c1 = Fraction(2, 15) * (1 - Fraction(1, 16 ** (2 ** (k - 1) - 1)))
            if cl.get(1) != c1 or cl.get(2) != c1:
                ok = False
    check("closed forms c_0 = (2/3)(1 - 4^-(2^k-1)), c_1 = c_2 = (2/15)(1 - 16^-(2^(k-1)-1)) reproduce the exact table F<=11", ok)
    check("identities 1/2 = (2/3)(1-1/4), 21/32 = (2/3)(1-4^-3), 5461/8192 = (2/3)(1-4^-7), 1/8 = (2/15)(1-1/16), 273/2048 = (2/15)(1-16^-3)",
          Fraction(1, 2) == Fraction(2, 3) * (1 - Fraction(1, 4)) and Fraction(21, 32) == Fraction(2, 3) * (1 - Fraction(1, 64)) and
          Fraction(5461, 8192) == Fraction(2, 3) * (1 - Fraction(1, 4 ** 7)) and Fraction(1, 8) == Fraction(2, 15) * (1 - Fraction(1, 16)) and
          Fraction(273, 2048) == Fraction(2, 15) * (1 - Fraction(1, 16 ** 3)))
    check("261/512 = 1/2 + 1/128 + 1/512", Fraction(261, 512) == Fraction(1, 2) + Fraction(1, 128) + Fraction(1, 512))
    # shift law in floating point for F <= Fmax: c_(z+8)(F) = c_z(F-8), z = 0,1,2,4, F <= 16; and what happens at F=17
    for F in range(11, Fmax[6] + 1):
        bz = byz_float[F]; bz8 = byz_float[F - 8]
        diffs = {z: bz.get(z + 8, 0.0) - bz8.get(z, 0.0) for z in (0, 1, 2, 4)}
        ok = all(abs(v) < 1e-6 for v in diffs.values())
        check("shift law c_(z+8)(F) = c_z(F-8), z in {0,1,2,4}, F=%d" % F, ok, "diffs %s" % {z: "%.2e" % v for z, v in diffs.items()})
        d4 = bz.get(4, 0.0) - byz_float[F - 4].get(0, 0.0)
        note("F=%d: c_4(F) - c_0(F-4) = %.6f (the note says c_4 = c_0(F-4) holds for F = 7, 8 only)" % (F, d4))
    for F in (7, 8):
        note("F=%d: c_4(F) = %.6f, c_0(F-4) = %.6f" % (F, byz_float[F].get(4, 0.0), byz_float[F - 4].get(0, 0.0)))
    check("c_4(F) = c_0(F-4) for F = 7, 8 (1/2) and fails at F = 9 (261/512 vs 21/32)",
          abs(byz_float[7][4] - 0.5) < 1e-9 and abs(byz_float[8][4] - 0.5) < 1e-9 and abs(byz_float[9][4] - 261 / 512) < 1e-9 and abs(byz_float[5][0] - 21 / 32) < 1e-9)
    ok356 = all(abs(byz_float[F].get(z, 0.0) - 1 / 128) < 1e-9 for F in range(9, Fmax[6] + 1) for z in (3, 5, 6))
    ok7 = all(abs(byz_float[F].get(7, 0.0)) < 1e-9 for F in range(9, Fmax[6] + 1))
    check("c_3 = c_5 = c_6 = 1/128 and c_7 = 0 for 9 <= F <= %d (float)" % Fmax[6], ok356 and ok7)
    okz0 = True
    for F in range(3, Fmax[6] + 1):
        k = 1
        while not (2 ** k < F <= 2 ** (k + 1)):
            k += 1
        if abs(byz_float[F][0] - float(Fraction(2, 3) * (1 - Fraction(1, 4 ** (2 ** k - 1))))) > 1e-9:
            okz0 = False
        if F >= 5 and (abs(byz_float[F][1] - float(Fraction(2, 15) * (1 - Fraction(1, 16 ** (2 ** (k - 1) - 1))))) > 1e-9 or abs(byz_float[F][2] - byz_float[F][1]) > 1e-9):
            okz0 = False
    check("closed forms for c_0, c_1 = c_2 hold in floating point for every computed F = 3..%d (including F = 17's new block if computed)" % Fmax[6], okz0)
    print("  per-z table (2^F * excess by z, float):")
    for F in range(3, Fmax[6] + 1):
        print("   F=%2d %s" % (F, {z: "%.6f" % v for z, v in sorted(byz_float[F].items()) if abs(v) > 1e-12}))
    # halving
    okh = True
    for F in range(3, Fmax[6], 2):
        e1 = data[(6, F)]['p'] - float(data[(6, F)]['fo']); e2 = data[(6, F + 1)]['p'] - float(data[(6, F + 1)]['fo'])
        if abs(e2 - e1 / 2) > 1e-13:
            okh = False
            note("halving fails at F=%d: %.6e vs %.6e" % (F, e1 / 2, e2))
    check("E(F+1) = E(F)/2 for odd F (float, to 1e-13) for F = 3..%d" % (Fmax[6] - 1), okh)
    # growth claims
    print("  growth: F, c_6(F), 1.4 log_2 F, ratio-1 = c_6/(2F), 0.7 log_2 F / F")
    for F in range(3, Fmax[6] + 1):
        c6 = (data[(6, F)]['p'] - float(data[(6, F)]['fo'])) * 2 ** F
        print("   %2d  %.4f  %.4f  %.4f  %.4f" % (F, c6, 1.4 * math.log2(F), c6 / (2 * F), 0.7 * math.log2(F) / F))
    c6_9 = (data[(6, 9)]['p'] - float(data[(6, 9)]['fo'])) * 2 ** 9
    c6_15 = (data[(6, Fmax[6])]['p'] - float(data[(6, Fmax[6])]['fo'])) * 2 ** Fmax[6]
    note("c_6(17) - c_6(9) from the note's floats = 2.8732 - 1.4664 = 1.4068: the '1.4 per doubling' rests on ONE doubling (9 -> 17); "
         "at F = 17, 1.4 log_2 17 = %.2f versus c_6(17) = 2.87, and 0.7 log_2 F / F = %.3f versus the observed ratio-1 = 0.085" % (1.4 * math.log2(17), 0.7 * math.log2(17) / 17))
    check("the quoted 'ratio - 1 ~ 0.7 log_2 F / F' is NOT matched by the quoted data points (0.083, 0.100, 0.096, 0.085 vs %.3f, %.3f, %.3f, %.3f)" % tuple(0.7 * math.log2(F) / F for F in (3, 7, 15, 17)),
          all(abs(0.7 * math.log2(F) / F - v) > 0.05 for F, v in ((3, 0.083), (7, 0.100), (15, 0.096), (17, 0.085))),
          "the formula is a conjectured asymptotic, not a fit of the F <= 17 data (ratios are flat 1.06-1.10)")


# ----------------------------------------------------------------------------------------------------------------------
# C. COLLATZ
# ----------------------------------------------------------------------------------------------------------------------
LOG23 = math.log2(3)


def syr_word(n, j):
    w = []
    for _ in range(j):
        m = 3 * n + 1; v = 0
        while m % 2 == 0:
            m //= 2; v += 1
        w.append(v); n = m
    return w, n


def S_of(word):
    """S_j = sum_(t<j) 3^(j-1-t) 2^(A_t), A_0 = 0"""
    j = len(word); A = 0; S = 0
    At = [0]
    for v in word:
        A += v; At.append(A)
    for t in range(j):
        S += 3 ** (j - 1 - t) * 2 ** At[t]
    return S, A


def collatz_part():
    print("=" * 110)
    print("C. COLLATZ PRECISION NOTE / THM-4512 / the two scripts")
    print("=" * 110)
    # ---- C1: formulas ----
    import random
    random.seed(20260926)
    ok_formula = True; ok_rec = True; ok_bound = True; tight = False
    for trial in range(3000):
        n = random.randrange(1, 1 << 40) | 1
        j = random.randint(1, 30)
        w, m = syr_word(n, j)
        S, A = S_of(w)
        if (3 ** j * n + S) % 2 ** A != 0 or (3 ** j * n + S) // 2 ** A != m:
            ok_formula = False
        # recursion S_(j+1) = 3 S_j + 2^(A_j)
        w1, _ = syr_word(n, j + 1); S1, A1 = S_of(w1)
        if S1 != 3 * S + 2 ** A:
            ok_rec = False
        # bound S_j <= 2^A ((3/2)^j - 1)  <=>  S_j 2^j <= 2^A (3^j - 2^j)
        if S * 2 ** j > 2 ** A * (3 ** j - 2 ** j):
            ok_bound = False
        if S * 2 ** j == 2 ** A * (3 ** j - 2 ** j):
            tight = True
    check("U^j(n) = (3^j n + S_j)/2^A with S_j = sum_(t<j) 3^(j-1-t) 2^(A_t) (3000 random n, j)", ok_formula)
    check("recursion S_(j+1) = 3 S_j + 2^(A_j)", ok_rec)
    check("bound S_j <= 2^A ((3/2)^j - 1) on random words (proof: A_t <= A-(j-t) since every v >= 1; equality iff all v = 1)", ok_bound)
    ok_eq = all(S_of([1] * j)[0] == 3 ** j - 2 ** j for j in range(1, 40))
    okb2 = True
    for trial in range(2000):
        n = random.randrange(1, 1 << 40) | 1; j = random.randint(1, 30)
        w, _ = syr_word(n, j); S, A = S_of(w); Ap = A - w[-1]
        if S * 2 ** j > 2 ** Ap * 2 * (3 ** j - 2 ** j):
            okb2 = False
    check("S_j <= 2^(A_(j-1)) * 2((3/2)^j - 1), i.e. S_j/2^A <= 2^-v_j 2((3/2)^j - 1) (the bound used in THM-4512's proof of part 4)", okb2)
    note("THM-4512's proof of part 4 compares S_j/2^A with 1 <= rho_w; the exception criterion is rho_w <= N(w) = S_j/(2^A - 3^j) > S_j/2^A, so the written inequality is the wrong quantity (repairable: 2^A - 3^j >= 2^A/2 once v >= vcrit+2)")
    check("equality case: all-ones word has S_j = 3^j - 2^j = 2^A((3/2)^j - 1)", ok_eq)
    # Terras class structure: the odd n with a given exact word form one class mod 2^(A+1), NOT mod 2^A
    okA1 = True; okA = False   # okA: whether "mod 2^A" ever suffices
    for word in ([1], [2], [1, 2], [1, 1, 3], [2, 1, 1, 2], [1, 2, 1, 2, 2], [3], [1, 1, 1, 5]):
        S, A = S_of(word); j = len(word)
        members = [n for n in range(1, 1 << (A + 4), 2) if syr_word(n, j)[0] == word]
        rho = (-S * pow(3, -j, 1 << A)) % (1 << A)
        classA1 = set(m % (1 << (A + 1)) for m in members)
        classA = set(m % (1 << A) for m in members)
        if len(classA1) != 1 or len(members) != (1 << (A + 4)) // (1 << (A + 1)):
            okA1 = False
        # is the class equal to {odd n = rho mod 2^A}?
        full = [n for n in range(1, 1 << (A + 4), 2) if n % (1 << A) == rho]
        if set(full) == set(members):
            okA = True
        note("word %s: A=%d, rho_w = -S 3^-j mod 2^A = %d; members < 2^(A+4): %s ...; class mod 2^(A+1) = %s; rho in class: %s" % (
            word, A, rho, members[:4], sorted(classA1), rho in members))
    check("odd n with a given exact Syracuse word form ONE residue class modulo 2^(A+1) (density 2^-A among odd n)", okA1)
    check("the class is NOT the full set of odd n = rho_w mod 2^A (THM-4512 'one residue class modulo 2^A' is off by one bit)", not okA)
    # ---- C2: gap inequality j <= 5000 ----
    t0 = time.time()
    worst = (Fraction(0), 0, 0); bad = []
    pow3 = 1;
    for j in range(1, 5001):
        pow3 *= 3
        A = pow3.bit_length()          # minimal A with 2^A > 3^j (3^j is never a power of 2)
        assert (1 << A) > pow3 and (1 << (A - 1)) < pow3
        gap = (1 << A) - pow3
        # (3/2)^j - 1 < gap  <=>  3^j - 2^j < gap 2^j
        if not (pow3 - (1 << j) < gap << j):
            bad.append((j, A))
        r = Fraction(pow3 - (1 << j), gap << j)
        if r > worst[0]:
            worst = (r, j, A)
    check("2^A - 3^j > (3/2)^j - 1 at the minimal admissible A = ceil(j log_2 3) for all j <= 5000 (exact integers, %.1fs)" % (time.time() - t0), not bad)
    check("worst ratio ((3/2)^j - 1)/(2^A - 3^j) = %.4f at (j, A) = (%d, %d) (claimed 0.507 at (5, 8))" % (float(worst[0]), worst[1], worst[2]), worst[1] == 5 and worst[2] == 8 and abs(float(worst[0]) - 0.5072) < 1e-3)
    note("Baker-type remark: with |A ln2 - j ln3| >= C j^-mu one gets 2^A - 3^j >= C 3^j j^-mu; the theorem states this without C and without a j_0; the 'for all j' clause is conditional as written")
    # ---- C3: rigorous exception search j <= 14 ----
    JMAX = 14
    t0 = time.time()
    pow3l = [3 ** i for i in range(JMAX + 2)]
    prefixes = 0; classes_rig = 0; exceptions = []; maxAseen = 0
    best = {}   # (j, A) -> (ratio, word)
    classes_theirs = 0

    def rec(j, A, S, word):
        nonlocal prefixes, classes_rig, maxAseen, classes_theirs
        prefixes += 1
        jn = j + 1
        Sn = 3 * S + (1 << A)
        p3 = pow3l[jn]
        # critical: largest v with 2^(A+v) < 3^jn
        vcrit = 0
        while (1 << (A + vcrit + 1)) < p3:
            vcrit += 1
        # first-descent valuations v >= vcrit+1: exception iff rho (2^An - 3^jn) <= Sn ; once 2^An - 3^jn > Sn no exception for any larger v
        v = vcrit + 1
        while True:
            An = A + v
            gap = (1 << An) - p3
            classes_rig += 1
            maxAseen = max(maxAseen, An)
            rho = (-Sn * pow(3, -jn, 1 << An)) % (1 << An)
            if rho == 0:
                rho = 1 << An
            ratio = Fraction(Sn, gap << An)
            key = (jn, An)
            if key not in best or ratio > best[key][0]:
                best[key] = (ratio, tuple(word + [v]))
            if rho * gap <= Sn:
                w, m = syr_word(rho, jn)
                exceptions.append((rho, tuple(word + [v]), w == word + [v], m))
            if gap > Sn:          # N(w) < 1 for this and every larger v: rigorous stop
                break
            v += 1
        # replicate the audited script's counting rule: v in range(vcrit+1, vcrit+40), break after An > 40
        for v2 in range(vcrit + 1, vcrit + 40):
            classes_theirs += 1
            if (1 << (A + v2)) > (1 << 40):
                break
        if jn < JMAX:
            for v2 in range(1, vcrit + 1):
                rec(jn, A + v2, Sn, word + [v2])

    rec(0, 0, 0, [])
    check("j <= 14: rigorous enumeration of all first-descent classes that could have an uncertified member (stop once 2^A - 3^j > S_j): %d classes over %d no-descent prefixes, max A = %d (%.1fs)" % (classes_rig, prefixes, maxAseen, time.time() - t0), True)
    check("j <= 14: the only class with rho_w <= N(w) is rho = 1, word (2) (found: %s)" % exceptions, exceptions == [(1, (2,), True, 1)])
    check("replication of the audited script's enumeration rule gives its quoted 606746 classes (rule = v up to vcrit+39 AND stop after A > 40, not 'all v up to 40 beyond critical')", classes_theirs == 606746, "%d" % classes_theirs)
    note("what the audited enumeration actually covers: for j=14 prefixes with A_13 = 13 it checks v = 10..28 (A <= 41), not 40 valuations; the conclusion survives because N(w) < 1 for every larger A (the rigorous stop above)")
    rows = sorted(best.items(), key=lambda kv: -kv[1][0])
    top = [(j, A, float(r), w) for (j, A), (r, w) in rows if r > 0.05]
    print("  thresholds N(w)/2^A > 0.05 by (j, A), max over ALL first-descent words j <= 14: %s" % top)
    check("threshold table: 0.25 at (1,2) word (2); 0.1437 at (3,5) word (1,2,2); 0.0959 at (5,8) word (1,2,1,2,2); nothing else above 0.05",
          len(top) == 3 and top[0][:2] == (1, 2) and abs(top[0][2] - 0.25) < 1e-9 and top[1][:2] == (3, 5) and abs(top[1][2] - 0.1437) < 1e-3 and top[1][3] == (1, 2, 2)
          and top[2][:2] == (5, 8) and abs(top[2][2] - 0.0959) < 1e-3 and top[2][3] == (1, 2, 1, 2, 2))
    note("N(w)/2^A at (12, 20) = %.2e, at (41, 65) = (checked below by the j<=5000 bound); the '(12, 19) [no descent]' remark in the script output is correct: 2^19 < 3^12" % float(max((r for (j, A), (r, w) in best.items() if j == 12 and A == 20), default=0)))
    # ---- C4: sigma vs sigma_inf (Syracuse) and tau_T vs sigma_T (T-coding), vectorised ----
    NMAX = 2 * 10 ** 6 if QUICK else 10 ** 7
    t0 = time.time()
    JT = 400
    descS = np.zeros((JT + 1, 3 * JT + 1), dtype=bool)   # descS[j, A] = 3^j < 2^A
    for j in range(JT + 1):
        p = 3 ** j
        for A in range(3 * JT + 1):
            descS[j, A] = (1 << A) > p
    nvals = np.arange(3, NMAX + 1, 2, dtype=np.int64)
    N = len(nvals)
    sig = np.zeros(N, dtype=np.int32); sig_inf = np.zeros(N, dtype=np.int32)
    maxpeak = 0
    CH = 1 << 20
    for c0 in range(0, N, CH):
        n = nvals[c0:c0 + CH]; x = n.copy(); A = np.zeros(len(n), dtype=np.int64)
        act = np.arange(len(n)); s = np.zeros(len(n), dtype=np.int32); si = np.zeros(len(n), dtype=np.int32)
        step = 0
        while len(act):
            step += 1
            xa = x[act]
            m = 3 * xa + 1
            low = m & (-m)
            v = (np.frexp(low.astype(np.float64))[1] - 1).astype(np.int64)
            xa = m >> v
            A[act] += v
            x[act] = xa
            maxpeak = max(maxpeak, int(xa.max()))
            coef = descS[step, np.minimum(A[act], 3 * JT)]
            newc = coef & (si[act] == 0)
            si[act[newc]] = step
            desc = xa < n[act]
            s[act[desc]] = step
            act = act[~desc]
            assert step < JT
        sig[c0:c0 + CH] = s; sig_inf[c0:c0 + CH] = si
    diff = np.nonzero(sig != sig_inf)[0]
    check("Syracuse: sigma(n) == sigma_inf(n) for every odd 3 <= n <= %d (%d numbers, %.0fs)" % (NMAX, N, time.time() - t0), len(diff) == 0, "differences: %s" % [(int(nvals[i]), int(sig[i]), int(sig_inf[i])) for i in diff[:10]])
    check("Syracuse: maximal stopping time for odd n <= %d is %d (claimed 155 for 10^7)" % (NMAX, int(sig.max())), int(sig.max()) == 155 or NMAX < 10 ** 7, "max at n = %d; max orbit value before descent %d (< 2^63 ok)" % (int(nvals[int(np.argmax(sig))]), maxpeak))
    check("sigma_inf <= sigma always (Terras inequality) on the range", bool((sig_inf <= sig).all()))
    # T-coded comparison
    t0 = time.time()
    KT = 700
    descT = np.zeros((KT + 1, KT + 1), dtype=bool)  # descT[k, o] = 3^o < 2^k
    for k in range(KT + 1):
        for o in range(min(k, KT) + 1):
            descT[k, o] = (3 ** o) < (1 << k)
    sigT = np.zeros(N, dtype=np.int32); tauT = np.zeros(N, dtype=np.int32)
    for c0 in range(0, N, CH):
        n = nvals[c0:c0 + CH]; x = n.copy(); o = np.zeros(len(n), dtype=np.int64)
        act = np.arange(len(n)); s = np.zeros(len(n), dtype=np.int32); tt = np.zeros(len(n), dtype=np.int32)
        k = 0
        while len(act):
            k += 1
            xa = x[act]
            odd = (xa & 1) == 1
            xa = np.where(odd, (3 * xa + 1) >> 1, xa >> 1)
            o[act] += odd
            x[act] = xa
            coef = descT[k, np.minimum(o[act], KT)]
            newc = coef & (tt[act] == 0)
            tt[act[newc]] = k
            desc = xa < n[act]
            s[act[desc]] = k
            act = act[~desc]
            assert k < KT
        sigT[c0:c0 + CH] = s; tauT[c0:c0 + CH] = tt
    diffT = np.nonzero(sigT != tauT)[0]
    check("T-coding (Terras's actual conjecture): tau_T(n) == sigma_T(n) for every odd 3 <= n <= %d (%.0fs); max sigma_T = %d" % (NMAX, time.time() - t0, int(sigT.max())), len(diffT) == 0,
          "differences (n, sigma_T, tau_T): %s (count %d)" % ([(int(nvals[i]), int(sigT[i]), int(tauT[i])) for i in diffT[:12]], len(diffT)))
    note("the note/THM check 'sigma = sigma_inf' is in Syracuse coding, which is coarser than Terras's T-coding: a T-coded exception with tau_T inside a Syracuse step (k* = ceil(j log_2 3) < A_j) would be invisible to it")
    # ---- C5: D(k), W_k exact DP and rates ----
    t0 = time.time()
    KMAX = 60
    KBIG = 400 if QUICK else 1500
    # Syracuse no-descent words: counts N[j][A]
    def D_series(K):
        cur = {0: 1}
        out = []
        p3 = 1
        for j in range(1, K + 1):
            p3 *= 3
            Amax = p3.bit_length() - 1          # largest A with 2^A < 3^j
            nxt = {}
            # nxt[A'] = sum_{A < A'} cur[A]
            keys = sorted(cur)
            pref = 0; ki = 0
            for Ap in range(1, Amax + 1):
                while ki < len(keys) and keys[ki] < Ap:
                    pref += cur[keys[ki]]; ki += 1
                if pref:
                    nxt[Ap] = pref
            cur = nxt
            out.append(cur)
        return out
    levels = D_series(KBIG)
    Dk = [sum(Fraction(c, 1 << A) for A, c in lv.items()) for lv in levels[:KMAX]]
    claimed_small = {1: Fraction(1, 2), 2: Fraction(3, 8), 3: Fraction(1, 4), 4: Fraction(13, 64), 5: Fraction(19, 128), 6: Fraction(1, 8), 7: Fraction(113, 1024), 8: Fraction(367, 4096)}
    check("D(k) exact for k = 1..8 = 1/2, 3/8, 1/4, 13/64, 19/128, 1/8, 113/1024, 367/4096", all(Dk[k - 1] == v for k, v in claimed_small.items()), "%s" % [str(Dk[k - 1]) for k in range(1, 9)])
    claimed_float = {12: 0.05212, 16: 0.03092, 20: 0.01925, 24: 0.01314, 28: 0.00875, 32: 0.00593, 36: 0.00427, 40: 0.00299, 41: 0.00265, 48: 0.00156, 52: 0.00112, 56: 0.00082, 60: 0.00062}
    okD = all(abs(float(Dk[k - 1]) - v) < 6e-6 for k, v in claimed_float.items())
    check("D(k) floats for k = 12..60 match the note's table (5 decimals)", okD, "D(41) = %.6e, D(60) = %.6e" % (float(Dk[40]), float(Dk[59])))
    check("D(41) = 0.00265..., 1 - D(41) = 0.9973 (certified density among odd n within 41 Syracuse steps)", abs(float(Dk[40]) - 2.654723e-3) < 1e-8 and abs(1 - float(Dk[40]) - 0.997345) < 1e-6)
    check("D(k)^(1/k): 0.865 at k=40, 0.884 at k=60 as quoted", abs(float(Dk[39]) ** (1 / 40) - 0.865) < 1e-3 and abs(float(Dk[59]) ** (1 / 60) - 0.884) < 1e-3, "%.4f, %.4f" % (float(Dk[39]) ** (1 / 40), float(Dk[59]) ** (1 / 60)))
    # rate at large k (log domain)
    def logD(lv):
        # log2 of sum c 2^-A : use max-scaling
        Amin = min(lv)
        tot = sum(Fraction(c, 1 << (A - Amin)) for A, c in lv.items())
        return -Amin + math.log2(float(tot.numerator)) - math.log2(float(tot.denominator)) if tot.numerator < (1 << 1000) else -Amin + (tot.numerator.bit_length() - tot.denominator.bit_length())
    l1 = logD(levels[KBIG - 1]); l2 = logD(levels[KBIG // 2 - 1])
    half = KBIG - KBIG // 2
    rate_syr = 2 ** ((l1 - l2) / half)                       # geometric mean of D(k+1)/D(k) over the upper half
    rate_syr_c = rate_syr * (KBIG / (KBIG // 2)) ** (1.5 / half)   # remove the k^-3/2 factor
    h = -(math.log(2, 3) * math.log2(math.log(2, 3)) + (1 - math.log(2, 3)) * math.log2(1 - math.log(2, 3)))
    pred_T = 2 ** (-(1 - h)); pred_syr = 3 ** (-(1 - h))
    note("h* = h(log_3 2) = %.7f; per-T-step 2^-(1-h*) = %.5f; per-Syracuse-step 2^-((1-h*) log_2 3) = 3^-(1-h*) = %.5f; (1-h*)/rho* = %.4f" % (h, pred_T, pred_syr, (1 - h) / math.log(2, 3)))
    note("D(k): geometric-mean rate of D(k+1)/D(k) over k in [%d, %d]: %.5f raw, %.5f after removing k^-3/2; D(k)^(1/k) at k = %d: %.5f" % (KBIG // 2, KBIG, rate_syr, rate_syr_c, KBIG, 2 ** (l1 / KBIG)))
    check("the note's claimed limit of D(k)^(1/k), 2^-(1-h*) = 0.9659, is WRONG for the Syracuse-coded D(k): the rate is 3^-(1-h*) = 0.9466 (the INDEX/synthesis wiring says exponent (1-h*)/rho* = 0.0793 per Syracuse step, i.e. 2^-0.0793 = 0.9466, contradicting the note)",
          abs(rate_syr_c - pred_syr) < 0.002 and abs(rate_syr_c - pred_T) > 0.01, "measured %.4f" % rate_syr_c)
    check("wiring: (1-h*)/rho* = 0.0793 (INDEX, synthesis 2s)", abs((1 - h) / math.log(2, 3) - 0.0793) < 5e-5)
    # W_k (T-coding)
    def W_series(K):
        p3 = [3 ** i for i in range(K + 2)]
        cur = {0: 1}; out = []
        for j in range(1, K + 1):
            nxt = {}
            p2 = 1 << j
            for o, c in cur.items():
                for letter in (0, 1):
                    oo = o + letter
                    if p3[oo] > p2:
                        nxt[oo] = nxt.get(oo, 0) + c
            cur = nxt; out.append(sum(cur.values()))
        return out
    W = W_series(KBIG)
    check("W_1..W_20 = 1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328 (nodescent note / THM-4495)",
          W[:20] == [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328])
    check("W_41 = 12805670000 and W_60 = 2216134944775156", W[40] == 12805670000 and W[59] == 2216134944775156, "%d, %d" % (W[40], W[59]))
    lw1 = math.log2(W[KBIG - 1]) - KBIG; lw2 = math.log2(W[KBIG // 2 - 1]) - KBIG // 2
    rW = 2 ** ((lw1 - lw2) / half); rWc = rW * (KBIG / (KBIG // 2)) ** (1.5 / half)
    note("W_k/2^k: geometric-mean rate over k in [%d, %d]: %.5f raw, %.5f after removing k^-3/2 (THM-4495 per-T-step rate %.5f)" % (KBIG // 2, KBIG, rW, rWc, pred_T))
    check("T-coded W_k/2^k rate agrees with 2^-(1-h*) = 0.9659 (THM-4495), while the Syracuse D(k) rate is 3^-(1-h*)", abs(rWc - pred_T) < 0.002 and abs(rate_syr_c - pred_T) > 0.01)
    note("D(k) and W_k index precision differently (k Syracuse steps ~ k log_2 3 T-steps at the no-descent barrier), as the note itself says; the exponent must be converted")
    print("  (%.1fs)" % (time.time() - t0))
    # ---- C6: swaplift bank ----
    t0 = time.time()
    path = os.path.join(ROOT, '05-knowledge', 'results', 'reset_20260926_swaplift.out')
    rows = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        lines = fh.read().splitlines()
    in_tab = False
    for ln in lines:
        if ln.startswith('certificate table:'):
            in_tab = True; continue
        if in_tab:
            parts = ln.split()
            if len(parts) == 9 and all(p.lstrip('-').isdigit() for p in parts):
                rows.append(tuple(int(p) for p in parts))
            else:
                if rows:
                    in_tab = False
    check("swaplift certificate table parsed: 171 rows", len(rows) == 171, "%d" % len(rows))
    maxJ = max(r[1] for r in rows); maxK = max(r[8] for r in rows); maxA = max(r[2] for r in rows); maxR = max(r[6] for r in rows)
    check("bank: max J = 40 (descent within J+1 <= 41 U steps), max A = 66, max R_* = 64, max K = 65 bits", (maxJ, maxA, maxR, maxK) == (40, 66, 64, 65), "%s" % ((maxJ, maxA, maxR, maxK),))
    # minimal disjoint cylinders
    cls = [(r[7], r[8], r[1], r[0]) for r in rows]   # residue, K, J, q
    keep = []
    for (res, K, J, q) in cls:
        contained = False
        for (res2, K2, J2, q2) in cls:
            if (res2, K2) == (res, K):
                continue
            if K2 <= K and (res - res2) % (1 << K2) == 0:
                contained = True; break
        if not contained:
            keep.append((res, K, J, q))
    keep = sorted(set((r, K) for r, K, _, _ in keep))
    dens = sum(Fraction(1, 1 << K) for _, K in keep)
    check("bank: 65 pairwise disjoint minimal cylinders", len(keep) == 65, "%d" % len(keep))
    check("bank density = 6985206796614369409/36893488147419103232 = 0.189334 among all integers, 0.378669 among odd", dens == Fraction(6985206796614369409, 36893488147419103232), "%.6f, x2 = %.6f" % (float(dens), 2 * float(dens)))
    check("all bank residues are odd", all(r % 2 == 1 for r, _ in keep))
    # containment in C_41: every member has an actual descent within J+1 steps (check a sample per class) => coefficient descent within 41 steps
    okd = True; okc = True; worstJ = 0
    Jof = {}
    for r in rows:
        Jof[(r[7], r[8])] = min(Jof.get((r[7], r[8]), 10 ** 9), r[1])
    for (res, K) in keep:
        J = Jof[(res, K)]
        for mlt in range(0, 12):
            n = res + (mlt << K)
            if n < 3:
                continue
            x = n; A = 0; dj = None; cj = None
            for j in range(1, 42):
                m = 3 * x + 1; v = 0
                while m % 2 == 0:
                    m //= 2; v += 1
                x = m; A += v
                if cj is None and (1 << A) > 3 ** j:
                    cj = j
                if x < n:
                    dj = j; break
            if dj is None or dj > J + 1:
                okd = False
            if cj is None or cj > 41 or (dj is not None and cj > dj):
                okc = False
            worstJ = max(worstJ, dj or 0)
    check("bank sample (12 members per class): actual descent within J+1 steps for every member tested (max observed %d)" % worstJ, okd)
    check("bank sample: each tested member has a coefficient descent within 41 Syracuse steps, at or before its actual descent (bank subset of C_41, THM-4512 part 1)", okc)
    # window-avoidance density: exact enumeration over valuation words
    def window_density(k, L):
        # words v_1..v_(L+k-1) with every window of length k starting at i < L having no coefficient descent; density sum 2^-sum v
        # enumerate with pruning: valuations bounded by the barrier (no-descent windows force v small); exact
        tot = Fraction(0)
        maxv = 2 * k + 2
        def ok_window(word):
            A = 0
            for j, v in enumerate(word, 1):
                A += v
                if (1 << A) > 3 ** j:
                    return False
            return True
        def rec(word):
            nonlocal tot
            m = len(word)
            # check the windows that are complete
            if m >= k:
                i = m - k
                if i < L and not ok_window(word[i:i + k]):
                    return
            if m == L + k - 1:
                tot += Fraction(1, 1 << sum(word))
                return
            for v in range(1, maxv + 1):
                rec(word + [v])
        rec([])
        return tot
    okw = True; detail = []
    for k in (1, 2, 3, 4):
        for L in (1, 2, 3):
            dd = window_density(k, L)
            pred = Fraction(2, 1 << L) * Dk[k - 1]
            detail.append((k, L, str(dd), str(pred)))
            if dd != pred or dd > Fraction(1, 1 << L) or (k == 1 and dd != Fraction(1, 1 << L)):
                okw = False
    check("window-avoidance density (first L Syracuse points all outside the k-step certified region) is exactly 2^(1-L) D(k): = 2^-L for k = 1 and < 2^-L for k >= 2 (the note's bound is correct but not sharp)", okw, "%s" % detail)
    print("  (%.1fs)" % (time.time() - t0))


def main():
    print("gilbreath_size6_collatz_precision_20260926_audit.py  (QUICK=%s)  numpy %s" % (QUICK, np.__version__))
    gilbreath_part()
    collatz_part()
    print("=" * 110)
    nfail = sum(1 for _, ok, _ in RESULTS if not ok)
    print("SUMMARY: %d checks, %d failed, total time %.0fs" % (len(RESULTS), nfail, elapsed()))
    for name, ok, detail in RESULTS:
        if not ok:
            print("  FAILED: %s -- %s" % (name, detail))


if __name__ == '__main__':
    main()
