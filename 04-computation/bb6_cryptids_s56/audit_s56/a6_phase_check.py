#!/usr/bin/env python3
"""Referee: A6 phase map vs the raw (TM-verified) rule system, on random (v0, w0) of many sizes, plus
halting tops T_d, capture family, N_d entry, and the cascade/end-case derivation."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), 'session_copy'))
from tm_rules_check import a6_rule, HALT
from a6_phase_map_20261025 import phase as s_phase


def fast(s):
    u, v, w = s
    if u >= 1 and w >= 2:
        t = min(u, w // 2)
        return (u - t, v + 2 * t, w - 2 * t)
    return a6_rule(s)


def raw_until(s0, target, budget=10 ** 6):
    """Follow the raw collapsed rules from s0 until `target` (a state or HALT) appears; report states with u == 0
    passed on the way (excluding s0), and whether a return J(0,0,.) or HALT occurred before target."""
    s = s0
    passed0 = []
    for n in range(budget):
        s = fast(s)
        if s == target:
            return True, passed0, n + 1
        if s == HALT:
            return False, passed0 + ['HALT'], n + 1
        if s[0] == 0:
            if s[1] == 0:
                return False, passed0 + [('ret', s[2])], n + 1
            passed0.append(s)
    return False, passed0, budget


def my_depth(v0, w0):
    """Count cascade levels by direct simulation of the raw rules: (0,v,w) with w >= 2v+9 maps to (0,2v+6,w-2v-9)."""
    v, w = v0, w0
    j = 0
    while w >= 2 * v + 9:
        s = fast((0, v, w))          # (v+3, 0, w-3)
        s = fast(s)                  # collapsed loop
        assert s[0] == 0, s
        v, w = s[1], s[2]
        j += 1
    return j, v, w


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    bad = []
    count = {}
    for trial in range(ntr):
        sz = rng.choice([8, 16, 64, 1000, 10 ** 6, 2 ** 40, 2 ** 100])
        v0 = rng.randrange(0, sz) if rng.random() < 0.7 else rng.choice([0, 8, rng.randrange(0, 20)])
        w0 = rng.randrange(0, sz * rng.choice([1, 4, 64, 2 ** 20]))
        res, j, lab = s_phase(v0, w0)
        jj, vj, wj = my_depth(v0, w0)
        if jj != j:
            bad.append(('depth', v0, w0, j, jj))
            continue
        if res[0] == 'HALT':
            target = HALT
        elif res[0] == 'ret':
            target = (0, 0, res[1])
        else:
            target = (0, res[1], res[2])
        ok, passed, n = raw_until((0, v0, w0), target)
        # passed must consist only of the cascade levels (and, for w2vo, J(0, v+4, 1))
        count[lab] = count.get(lab, 0) + 1
        if not ok:
            bad.append(('reach', v0, w0, res, lab, passed[-3:]))
            continue
        # the cascade levels are u==0 states; check that every passed u==0 state before the end is a cascade level
        exp_levels = jj
        extra = len(passed) - exp_levels
        if extra not in (0, 1) or (extra == 1 and lab != 'w2vo'):
            bad.append(('extra-u0-states', v0, w0, lab, passed[-3:]))
        # halting criterion: HALT iff w_j = 2 v_j + 8
        if (res[0] == 'HALT') != (wj == 2 * vj + 8):
            bad.append(('halt-crit', v0, w0, vj, wj))
    print('phase map vs raw rules:', ntr, 'trials; labels seen', sorted(count.items()), '; bad', len(bad), bad[:8])
