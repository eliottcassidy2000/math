#!/usr/bin/env python3
"""A6: exact phase map and halting criterion, verified against the paper's rule system (S56).

Phase start P(v0, w0) = configuration J(0, v0, w0).  With A = v0 + 6 the cascade
    J(0,v,w) -> J(v+3,0,w-3) -> J(0, 2v+6, w-2v-9)        (while w >= 2v + 9)
gives v_j = A 2^j - 6, w_j = w0 - 2A(2^j - 1) + 3j; the phase ends at the first j with w_j < 2A 2^j - 3.
End cases at J(0, v, w) (w < 2v + 9):
  w >= 3: t = (w-3)//2, u = v + 3 - t (>= 1), eps = (w-3) % 2.
     eps=1: u=1 HALT | u=2 return 4t+10 | u even>=4 -> P(u+4, 4t+1) | u odd>=3 -> P(u, 4t+5)
     eps=0: u=1 return 4t+12 | u even -> P(8, 2u+4t+1) | u odd>=3 return 2u+4t+10
  w <= 2: v even: w=0 return 2v+12 | w=1 -> P(8, 2v+5) | w=2 return 2v+18
          v odd : w=0 -> P(8, 2v+3) | w=1 return 2v+14 | w=2 return 2v+22
("return W'" = J(0,0,W') = P(0, W').)
Halting criterion (THM candidate): a phase P(v0, w0) halts iff w0 = 4A 2^j - 2A - 3j - 4 for its end depth j
(equivalently w_j = 2 v_j + 8 at the end).  For v0 = 0: W = 24 2^j - 3j - 16 (the main chain of the halting set).
Usage: python a6_phase_map_20261025.py [n_returns]
"""
import sys, json, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, HALT


def phase(v0, w0):
    """Return ('HALT',) | ('ret', W') | ('P', v0', w0'), plus the depth j and end-case label."""
    A = v0 + 6
    # smallest j >= 0 with w0 - 2A(2^j - 1) + 3j < 2A 2^j - 3  <=>  w0 + 2A + 3j + 3 < 4A 2^j
    j = 0
    if w0 + 2 * A + 3 >= 4 * A:
        j = max(0, (w0 + 2 * A + 3) // (4 * A)).bit_length() - 1
        while j > 0 and w0 + 2 * A + 3 * (j - 1) + 3 < 4 * A * (1 << (j - 1)):
            j -= 1
        while not (w0 + 2 * A + 3 * j + 3 < 4 * A * (1 << j)):
            j += 1
    v = A * (1 << j) - 6
    w = w0 - 2 * A * ((1 << j) - 1) + 3 * j
    if w >= 3:
        t = (w - 3) // 2
        u = v + 3 - t
        eps = (w - 3) % 2
        if eps == 1:
            if u == 1:
                return ('HALT',), j, 'e1u1'
            if u == 2:
                return ('ret', 4 * t + 10), j, 'e1u2'
            if u % 2 == 0:
                return ('P', u + 4, 4 * t + 1), j, 'e1ueven'
            return ('P', u, 4 * t + 5), j, 'e1uodd'
        if u == 1:
            return ('ret', 4 * t + 12), j, 'e0u1'
        if u % 2 == 0:
            return ('P', 8, 2 * u + 4 * t + 1), j, 'e0ueven'
        return ('ret', 2 * u + 4 * t + 10), j, 'e0uodd'
    if v % 2 == 0:
        if w == 0:
            return ('ret', 2 * v + 12), j, 'w0ve'
        if w == 1:
            return ('P', 8, 2 * v + 5), j, 'w1ve'
        return ('ret', 2 * v + 18), j, 'w2ve'
    if w == 0:
        return ('P', 8, 2 * v + 3), j, 'w0vo'
    if w == 1:
        return ('ret', 2 * v + 14), j, 'w1vo'
    return ('ret', 2 * v + 22), j, 'w2vo'


def sweep_phase(W):
    """From J(0,0,W) to the next return J(0,0,W') via phases; returns (W' or 'HALT', signature)."""
    v0, w0 = 0, W
    sig = []
    for _ in range(10 ** 6):
        res, j, lab = phase(v0, w0)
        sig.append((v0 if v0 in (0, 8) else 'x', j, lab))
        if res[0] == 'HALT':
            return 'HALT', tuple(sig)
        if res[0] == 'ret':
            return res[1], tuple(sig)
        v0, w0 = res[1], res[2]
    raise RuntimeError('too many phases')


def sweep_rules(W, cap=10 ** 8):
    p = step('A6', (0, 0, W))
    n = 1
    while p != HALT and not (p[0] == 0 and p[1] == 0) and n < cap:
        p = step('A6', p)
        n += 1
    return 'HALT' if p == HALT else p[2]


if __name__ == '__main__':
    nret = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    # (1) verify phase map against the rule system on all W < 20000
    bad = [W for W in range(0, 20000) if sweep_phase(W)[0] != sweep_rules(W)]
    print(f'phase map == rule system for all W < 20000: {not bad} {bad[:5]}')
    # (2) halting criterion for v0 = 0 single phase: W = 24 2^j - 3j - 16
    chain = [24 * (1 << j) - 3 * j - 16 for j in range(0, 14)]
    print('main chain halts:', all(sweep_phase(W)[0] == 'HALT' for W in chain), chain[:8])
    # (3) the orbit
    W = 0
    t0 = time.time()
    n = 0
    halted = False
    labs = {}
    while n < nret:
        r, sig = sweep_phase(W)
        if r == 'HALT':
            halted = True
            break
        W = r
        n += 1
    print(f'orbit: {n} returns, halted {halted}, W has {W.bit_length()} bits ({time.time() - t0:.1f}s)')
    json.dump({'verified_below_20000': not bad, 'returns': n, 'halted': halted, 'W_bits': W.bit_length()},
              open(os.path.join(HERE, 'a6_phase_map.json'), 'w'), indent=1)
