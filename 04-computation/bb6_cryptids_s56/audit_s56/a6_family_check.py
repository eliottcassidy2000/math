#!/usr/bin/env python3
"""Referee: T_d halting tops, N_d entry, capture family P(v0,4v0+29)/P(v0,4v0+33) on the raw TM-verified rules."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_phase_check import fast, my_depth, raw_until
from tm_rules_check import HALT


def run_raw(s, max_rules):
    """Raw collapsed rules from s; stop at HALT or a return J(0,0,.) (after >= 1 rule). Return (status, n, s)."""
    for n in range(1, max_rules + 1):
        s = fast(s)
        if s == HALT:
            return 'HALT', n, s
        if s[0] == 0 and s[1] == 0:
            return 'RET', n, s
    return 'RUNNING', max_rules, s


def first_phase_out(v0, w0):
    """Raw: from J(0,v0,w0), run the cascade (my_depth) and then raw rules to the next u == 0 state."""
    j, v, w = my_depth(v0, w0)
    s = (0, v, w)
    for n in range(100):
        s = fast(s)
        if s == HALT or s[0] == 0:
            return j, s
    return j, None


if __name__ == '__main__':
    # (1) T_d
    ok = True
    for d in range(0, 70):
        T = 24 * 2 ** d - 3 * d - 16
        jd, _, _ = my_depth(0, T)
        j1, _, _ = my_depth(0, T + 1)
        st, n, s = run_raw((0, 0, T), 10 ** 6)
        if not (jd == d and j1 == d + 1 and st == 'HALT'):
            ok = False
            print('T_d FAIL', d, jd, j1, st)
    print('T_d (d = 0..69): depth d, T_d + 1 has depth d+1, raw rules halt in the first sweep:', ok)
    # largest W of depth d: check a few W just below T_d too (depth d) -- covered by monotonicity
    # (2) N_d
    ok = True
    for d in range(1, 70):
        N = 18 * 2 ** d - 3 * d - 4
        j, s = first_phase_out(0, N)
        exp = (0, 3 * 2 ** d - 5, 12 * 2 ** d + 13)
        if d == 1:
            if s != HALT:
                ok = False
                print('N_1 not halting?', s)
            continue
        if not (j == d and s == exp):
            ok = False
            print('N_d FAIL', d, j, s, exp)
    print('N_d (d = 2..69): first phase at depth d ends in P(3*2^d-5, 12*2^d+13); N_1 = T_1 = 29 halts:', ok)
    # (3) capture family on raw rules
    for v0 in list(range(0, 14)) + [40, 41, 1000, 1001, 2 ** 50, 2 ** 50 + 1]:
        w0 = 4 * v0 + (29 if v0 % 2 == 0 else 33)
        st, n, s = run_raw((0, v0, w0), 200000)
        # follow phases explicitly through the family for 30 phases
        chain = []
        v, w = v0, w0
        good = True
        for k in range(30):
            j, s2 = first_phase_out(v, w)
            if s2 == HALT or s2 is None or (s2[0] == 0 and s2[1] == 0):
                chain.append((j, s2))
                good = False
                break
            vv, ww = s2[1], s2[2]
            if v % 2 == 0:
                expo = (v + 1, 4 * (v + 1) + 33)
            else:
                expo = (v + 3, 4 * (v + 3) + 29)
            if (vv, ww) != expo or j != 1:
                chain.append(('off-family', j, (vv, ww), expo))
                good = False
                break
            v, w = vv, ww
        print(f'v0={v0} P({v0},{w0}): raw 200000 rules status {st}; 30-phase family chain {"OK" if good else chain}')
