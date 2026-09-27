#!/usr/bin/env python3
"""gilbreath_size6_exact_20260926.py -- exact extinction probabilities of a lone defect of size d in Gilbreath's
automaton via an absorbing Markov chain (session gilbreath6-collatz-precision-20260926, opus, 2026-09-26).

Columns 1..F-1 hold the left sea word (uniform in {0,2}), column F holds d, and the right sea enters only through
column F+1, whose successive values are i.i.d. fair {0,2} bits (each new row's value at F+1 involves a fresh sea bit
with coefficient 1 in the Lucas kernel). So the state (values at columns 1..F) is a Markov chain with two equally
likely transitions new[i] = |a[i]-a[i+1]| (i < F), new[F] = |a[F]-u|, u in {0,2}. Absorbing: a[1] >= 4 (the leading
1 is destroyed at the next row; probability 1) and max(a) <= 2 (no defect left; probability 0). The extinction
probability p_d(F) is the average over left words of the absorption probability, an exact rational.
Exact rationals (sparse Gaussian elimination over Fractions) for small F, floating value iteration for larger F,
and for d = 6 the excess over the front-only tail decomposed by z = number of zeros between the defect and the
first 2 on its left (the wall), and by the weight of the front path.
Usage: python3 gilbreath_size6_exact_20260926.py [Fmax_exact=11] [Fmax_float=17]
"""
import sys, time
from fractions import Fraction
from math import comb
import numpy as np


def step(state, u):
    F = len(state)
    new = [abs(state[i] - state[i + 1]) for i in range(F - 1)]
    new.append(abs(state[F - 1] - u))
    return tuple(new)


def reachable(initials):
    seen = {}
    order = []
    for s in initials:
        if s not in seen:
            seen[s] = len(order); order.append(s)
    i = 0
    while i < len(order):
        s = order[i]; i += 1
        if s[0] >= 4 or max(s) <= 2:
            continue
        for u in (0, 2):
            t = step(s, u)
            if t not in seen:
                seen[t] = len(order); order.append(t)
    return seen, order


def solve_float(seen, order, iters=200000, tol=1e-15):
    n = len(order)
    absorb = np.zeros(n); mask = np.zeros(n, dtype=bool)
    n0 = np.zeros(n, dtype=np.int64); n1 = np.zeros(n, dtype=np.int64)
    for idx, s in enumerate(order):
        if s[0] >= 4:
            absorb[idx] = 1.0; mask[idx] = True
        elif max(s) <= 2:
            mask[idx] = True
        else:
            n0[idx] = seen[step(s, 0)]; n1[idx] = seen[step(s, 2)]
    P = absorb.copy()
    for it in range(iters):
        Pn = np.where(mask, absorb, 0.5 * (P[n0] + P[n1]))
        if np.max(np.abs(Pn - P)) < tol:
            P = Pn; break
        P = Pn
    return P


def solve_exact(seen, order, seconds=600):
    t0 = time.time()
    trans = [i for i, s in enumerate(order) if not (s[0] >= 4 or max(s) <= 2)]
    col = {i: k for k, i in enumerate(trans)}
    rows = []
    for i in trans:
        s = order[i]
        eq = {col[i]: Fraction(1)}; rhs = Fraction(0)
        for u in (0, 2):
            j = seen[step(s, u)]
            t = order[j]
            if t[0] >= 4:
                rhs += Fraction(1, 2)
            elif max(t) <= 2:
                pass
            else:
                eq[col[j]] = eq.get(col[j], Fraction(0)) - Fraction(1, 2)
        rows.append((eq, rhs))
    pivots = {}
    for eq, rhs in rows:
        while True:
            if time.time() - t0 > seconds:
                return None
            ks = [k for k in eq if k in pivots]
            if not ks:
                break
            k = ks[0]
            peq, prhs = pivots[k]
            f = eq[k] / peq[k]
            for kk, v in peq.items():
                eq[kk] = eq.get(kk, Fraction(0)) - f * v
            rhs -= f * prhs
            eq = {kk: v for kk, v in eq.items() if v != 0}
        if not eq:
            assert rhs == 0
            continue
        k = min(eq)
        pivots[k] = (eq, rhs)
    sol = {}
    for k in sorted(pivots, reverse=True):
        eq, rhs = pivots[k]
        val = rhs
        for kk, v in eq.items():
            if kk != k:
                val -= v * sol[kk]
        sol[k] = val / eq[k]
    P = {}
    for i, s in enumerate(order):
        if s[0] >= 4:
            P[i] = Fraction(1)
        elif max(s) <= 2:
            P[i] = Fraction(0)
        else:
            P[i] = sol[col[i]]
    return P


def front_only(d, F):
    j = d // 2
    return Fraction(sum(comb(F - 1, t) for t in range(0, j - 1)), 2 ** (F - 1))


def path_weight(bits, F):
    """number of twos met by the first front: L_s(b) = XOR_(j subset s) b(F-1-s+j), s = 0..F-2; bits[i] = column i+1"""
    w = 0
    for s in range(F - 1):
        col = F - 1 - s
        v = 0; sub = s
        while True:
            v ^= bits[col - 1 + sub]
            if sub == 0:
                break
            sub = (sub - 1) & s
        w += v
    return w


def main():
    Fmax_exact = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    Fmax_float = int(sys.argv[2]) if len(sys.argv) > 2 else 17
    for d in (4, 6, 8):
        print("== d = %d ==" % d, flush=True)
        for F in range(3, Fmax_float + 1):
            initials = []
            for bits in range(1 << (F - 1)):
                s = tuple(2 * ((bits >> i) & 1) for i in range(F - 1)) + (d,)
                initials.append(s)
            t0 = time.time()
            seen, order = reachable(initials)
            P = solve_float(seen, order)
            pf = float(np.mean([P[seen[s]] for s in initials]))
            fo = front_only(d, F)
            line = " F=%2d states=%8d  p=%.12f  front-only=%.12f  ratio=%.6f  excess=%.6e  (%.1fs)" % (F, len(order), pf, float(fo), pf / float(fo) if fo else float('nan'), pf - float(fo), time.time() - t0)
            Pe = None
            if F <= Fmax_exact and d in (6, 8):
                t0 = time.time()
                Pe = solve_exact(seen, order, seconds=900)
                if Pe is not None:
                    pe = sum(Pe[seen[s]] for s in initials) / len(initials)
                    assert abs(float(pe) - pf) < 1e-9, (pe, pf)
                    ex = pe - fo
                    line += "\n        exact p = %s ; excess = %s = %s / 2^%d  (%.1fs)" % (pe, ex, ex * 2 ** (2 * F), 2 * F, time.time() - t0)
                else:
                    line += "\n        exact solve timed out"
            print(line, flush=True)
            if d == 6 and F <= 16:
                # decomposition of the excess by z (zeros between the wall and the defect) and by path weight
                byz = {}; byw = {}
                for bits in range(1 << (F - 1)):
                    b = [(bits >> i) & 1 for i in range(F - 1)]
                    s = tuple(2 * x for x in b) + (d,)
                    z = 0
                    while z < F - 1 and b[F - 2 - z] == 0:
                        z += 1
                    w = path_weight(b, F)
                    pr = P[seen[s]]
                    fo_ind = 1.0 if w <= 1 else 0.0
                    byz[z] = byz.get(z, 0.0) + (pr - fo_ind) / (1 << (F - 1))
                    byw[w] = byw.get(w, 0.0) + (pr - fo_ind) / (1 << (F - 1))
                print("        excess by z (zeros before the wall): %s" % {z: "%.4e" % v for z, v in sorted(byz.items())})
                print("        2^F * excess by z:                   %s" % {z: "%.6f" % (v * 2 ** F) for z, v in sorted(byz.items())})
                if Pe is not None:
                    byz_e = {}
                    for bits in range(1 << (F - 1)):
                        b = [(bits >> i) & 1 for i in range(F - 1)]
                        s_ = tuple(2 * x for x in b) + (d,)
                        z = 0
                        while z < F - 1 and b[F - 2 - z] == 0:
                            z += 1
                        w = path_weight(b, F)
                        byz_e[z] = byz_e.get(z, Fraction(0)) + (Pe[seen[s_]] - (1 if w <= 1 else 0)) / 2 ** (F - 1)
                    print("        exact 2^F * excess by z:             %s" % {z: str(v * 2 ** F) for z, v in sorted(byz_e.items())})
                print("        excess by first-front path weight:   %s" % {w: "%.3e" % v for w, v in sorted(byw.items())}, flush=True)


if __name__ == '__main__':
    main()
