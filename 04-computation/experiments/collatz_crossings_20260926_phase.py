#!/usr/bin/env python3
"""collatz_crossings_20260926_phase.py -- the real-place rigidity of an orbit: mantissa rotation, the two-place
potential Q_l, and what they say about crossings (session collatz-crossings-20260926, opus, 2026-09-26).

Odd iterates m_0 = n, m_(l+1) = (3 m_l + 1)/2^(v_l), d_l = v_0 + ... + v_(l-1), alpha = log_2 3, Delta_l = d_l - l alpha.
Exact identity (Prop. 6 / Bernstein):  2^(d_l) m_l = 3^l n + S_(l-1),  S_(l-1) = sum_(j<l) 3^(l-1-j) 2^(d_j),
so  m_l = n 2^(-Delta_l) C_l  with  C_l = prod_(j<l) (1 + 1/(3 m_j)) = 1 + R_l/n,  R_l = S_(l-1)/3^l.
(P) Phase: {log_2 m_l} = {log_2 n + l alpha + log_2 C_l}. Along any orbit with sum 1/m_j < infinity (every
    non-eventually-periodic orbit, THM-4476 Cor. 6) C_l converges, so the mantissa phases are an irrational
    rotation by alpha up to a convergent shift: equidistributed (Benford's law along the orbit). We print the
    discrepancy of the phases against uniform and against the pure rotation for long segments.
(Q) Potential: Q_l = n C_l = m_l 2^(Delta_l) is increasing, and 2-adically |Q_l|_2 = 2^(-d_l) exactly
    (Q_l = 2^(d_l) m_l / 3^l). We verify the identity exactly with Fractions and print C_l, R_l/n on segments.
(X) Crossings: for a level Y, the odd iterates' visits below Y: we print, for the record-holder orbits and
    5n+1 (a proxy for divergence), the number of downcrossings of dyadic levels versus the number of odd
    iterates below them, and the segment Diophantine pairs (M, D) between successive visits to a band.
Usage: python3 collatz_crossings_20260926_phase.py
"""
import math
from fractions import Fraction

ALPHA = math.log2(3)


def odd_iterates(n, steps, q=3):
    ms = [n]; ds = [0]; d = 0
    for _ in range(steps):
        x = q * n + 1; v = 0
        while x % 2 == 0:
            x //= 2; v += 1
        d += v; n = x
        ms.append(n); ds.append(d)
        if n == 1: break
    return ms, ds


def discrepancy(xs):
    xs = sorted(x % 1 for x in xs); N = len(xs)
    return max(max(abs((i + 1) / N - x), abs(i / N - x)) for i, x in enumerate(xs))


def main():
    print("(Q) exact two-place identity 2^(d_l) m_l = 3^l n + S_(l-1), and Q_l = n C_l = m_l 2^(Delta_l) increasing, |Q_l|_2 = 2^(-d_l):")
    for n in (27, 703, 26623):
        ms, ds = odd_iterates(n, 200)
        S = 0; ok = True; Cl = Fraction(1); prevQ = None; incr = True
        for l in range(1, len(ms)):
            S = 3 * S + 2 ** ds[l - 1]
            if 2 ** ds[l] * ms[l] != 3 ** l * n + S: ok = False
            Cl *= (1 + Fraction(1, 3 * ms[l - 1]))
            Q = Fraction(2 ** ds[l] * ms[l], 3 ** l)
            if Q != n * Cl: ok = False
            if prevQ is not None and Q <= prevQ: incr = False
            prevQ = Q
        print("   n=%6d: %d odd iterates; identity exact: %s; Q_l increasing: %s; C_final = %.6f, R_final/n = %.6f" % (n, len(ms) - 1, ok, incr, float(Cl), float(Cl - 1)))
    print("(P) phase equidistribution: discrepancy D_N of {log_2 m_l} (uniform) and of {log_2 m_l - l alpha} (should converge, not equidistribute):")
    cases = [("3n+1 from 2^60-1", odd_iterates(2 ** 60 - 1, 3000), 3), ("3n+1 from 2^200-1", odd_iterates(2 ** 200 - 1, 6000), 3),
             ("5n+1 from 7 (alpha = log2 5)", odd_iterates(7, 3000, q=5), 5), ("3n+1 record 63728127", odd_iterates(63728127, 3000), 3)]
    for name, (ms, ds), q in cases:
        alpha_q = math.log2(q)
        ph = [math.log2(m) for m in ms]
        N = len(ph)
        # the shift log2(C_l) = log2(m_l) - log2(n) + d_l - l log2 q should converge (C_l = prod(1 + 1/(q m_j)))
        shift = [ph[l] - ph[0] + ds[l] - l * alpha_q for l in range(N)]
        print("   %-24s N=%5d  D_N(phase)=%.4f (uniform ~ %.4f)  shift log2 C_l: first %.5f, last %.5f, max-min over last half %.2e" % (
            name, N, discrepancy(ph), 1 / math.sqrt(N), shift[1] if N > 1 else 0, shift[-1], max(shift[N // 2:]) - min(shift[N // 2:])))
    print("(X) crossings of dyadic levels by the odd iterates (record-holder 3n+1 orbits and 5n+1): level 2^t, #odd iterates below, #downcrossings (m_(l-1) > 2^t >= m_l)")
    for name, (ms, ds) in [("3n+1 from 63728127", odd_iterates(63728127, 3000)), ("3n+1 from 2^60-1", odd_iterates(2 ** 60 - 1, 3000)), ("5n+1 from 7", odd_iterates(7, 2000, q=5))]:
        row = []
        for t in (10, 14, 18, 22, 26, 30):
            Y = 2 ** t
            below = sum(1 for m in ms if m <= Y)
            down = sum(1 for l in range(1, len(ms)) if ms[l - 1] > Y >= ms[l])
            row.append("2^%d: %d/%d" % (t, below, down))
        print("   %-22s " % name + "  ".join(row))
        # Diophantine pairs between successive visits of the band (Y/2, Y] for Y = 2^18
        Y = 2 ** 18; vis = [l for l, m in enumerate(ms) if Y // 2 < m <= Y]
        pairs = [(vis[i + 1] - vis[i], ds[vis[i + 1]] - ds[vis[i]]) for i in range(min(8, len(vis) - 1))]
        print("      band (2^17, 2^18]: visits at l = %s; (M, D) between successive visits: %s; |M alpha - D| = %s" % (
            vis[:9], pairs, ["%.2f" % abs(M * ALPHA - D) for M, D in pairs]))


if __name__ == '__main__':
    main()
