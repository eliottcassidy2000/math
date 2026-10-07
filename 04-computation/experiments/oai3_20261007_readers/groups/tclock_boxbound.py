"""Certified lower bounds for merge probabilities of the T-clock relation chain (lane 'groups').

h_B(state) = P(reach (0,0) before leaving the finite box B(J, C) = {|j| <= J, |N| <= C 3^|j|}) is a lower bound for
the merge probability h(state).  Value iteration from 0 is monotone: after n sweeps it equals
P(merge within n steps without leaving B), a certified lower bound (float64; rounding ~1e-12 relative).
Targets:  M1  x = y + 1                       h(0, 1)
          S19 lag-1 Mersenne pair (p = 3q + 2) h(2, 1)  [state after the forced prefix 1,0,0,0]
          and min over the small box B(1,10).
Usage: python3 tclock_boxbound.py J C SWEEPS
"""
import sys, time
import numpy as np


def build(J, C):
    offs = {}
    tot = 0
    for j in range(-J, J + 1):
        R = C * 3 ** abs(j)
        offs[j] = (tot, R)
        tot += 2 * R + 1
    n = tot
    succ = np.full((2, n), -1, dtype=np.int64)   # -1 = exit; -2 = merge
    for j in range(-J, J + 1):
        base, R = offs[j]
        N = np.arange(-R, R + 1, dtype=object)  # exact python ints
        idx = base + np.arange(2 * R + 1)
        for p in (0, 1):
            jn = np.empty(len(N), dtype=np.int64)
            Nn = np.empty(len(N), dtype=object)
            for k, v in enumerate(N):
                e = v & 1
                if p == 0 and e == 0:
                    a, b = j, v >> 1
                elif p == 0:
                    a, b = (j + 1, (3 * v + 1) >> 1) if j >= 0 else (j + 1, (v + 3 ** (-j - 1)) >> 1)
                elif e == 0:
                    a, b = (j, (3 * v + 1 - 3 ** j) >> 1) if j >= 0 else (j, (3 * v + 3 ** (-j) - 1) >> 1)
                else:
                    a, b = (j - 1, (v - 3 ** (j - 1)) >> 1) if j >= 1 else (j - 1, (3 * v - 1) >> 1)
                jn[k] = a; Nn[k] = b
            for k in range(len(N)):
                a, b = int(jn[k]), Nn[k]
                if a == 0 and b == 0:
                    succ[p, idx[k]] = -2
                elif abs(a) <= J and abs(b) <= C * 3 ** abs(a):
                    ob, oR = offs[a]
                    succ[p, idx[k]] = ob + b + oR
    return offs, succ


def main():
    J, C, SW = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    offs, succ = build(J, C)
    n = succ.shape[1]
    zero = offs[0][0] + offs[0][1]          # index of (0,0)
    print('box B(%d,%d): %d states (built in %.0f s)' % (J, C, n, time.time() - t0))
    h = np.zeros(n)
    s0, s1 = succ[0], succ[1]
    m0 = s0 >= 0; m1 = s1 >= 0
    g0 = (s0 == -2).astype(float); g1 = (s1 == -2).astype(float)
    i0 = np.where(m0, s0, 0); i1 = np.where(m1, s1, 0)
    def at(j, N):
        b, R = offs[j]
        return b + N + R
    report = {'M1 (0,1)': at(0, 1), 'S19 lag-1 after prefix (2,1)': at(2, 1), '(1,1)': at(1, 1), '(0,-1)': at(0, -1)}
    checkpoints = sorted(set([SW // 16, SW // 4, SW]))
    small = [at(j, N) for j in (-1, 0, 1) for N in range(-10 * 3 ** abs(j), 10 * 3 ** abs(j) + 1) if (j, N) != (0, 0)]
    for it in range(1, SW + 1):
        hn = 0.5 * (np.where(m0, h[i0], 0.0) + g0) + 0.5 * (np.where(m1, h[i1], 0.0) + g1)
        hn[zero] = 1.0
        h = hn
        if it in checkpoints:
            print('  after %6d sweeps: ' % it + ';  '.join('%s >= %.5f' % (k, h[v]) for k, v in report.items())
                  + ';  min over B(1,10) >= %.5f' % min(h[small]))
    print('  (each value is P(merge within %d steps without leaving B(%d,%d)), a certified lower bound)' % (SW, J, C))


if __name__ == '__main__':
    main()
