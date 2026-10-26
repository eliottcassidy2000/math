#!/usr/bin/env python3
"""Referee: independent implementation of the A6 simplified rule system (Dos Santos et al. A.6, second system)
transcribed directly from the paper text, plus a direct Turing-machine check of every rule branch on random
parameters.  Written independently of the session code.

J(u, v, w) = 0^inf 1^5 10^u 1 A> 11 01^v 1^w 01 0^inf
"""
import random, sys

HALT = 'HALT'


def rule(u, v, w):
    """One rule application of the paper's simplified A.6 system (no collapsing)."""
    if u == 0:
        if w == 0:
            return (0, 0, 2 * v + 12) if v % 2 == 0 else (2, 4, 2 * v + 7)
        if w == 1:
            return (2, 4, 2 * v + 9) if v % 2 == 0 else (0, 0, 2 * v + 14)
        if w == 2:
            return (2, v + 4, 1) if v % 2 == 0 else (0, v + 4, 1)
        return (v + 3, 0, w - 3)
    # u >= 1
    if w >= 2:
        return (u - 1, v + 2, w - 2)
    if w == 0:
        if u == 1:
            return (0, 0, 2 * v + 12)
        if u % 2 == 0:
            return (2, 4, 2 * u + 2 * v + 5)
        return (0, 0, 2 * u + 2 * v + 10)
    # w == 1
    if u == 1:
        return HALT
    if u == 2:
        return (0, 0, 2 * v + 10)
    if u % 2 == 0:
        return (2, u, 2 * v + 5)
    return (0, u, 2 * v + 5)


def rule_fast(u, v, w):
    """Same, with the loop J(u,v,w)->J(u-1,v+2,w-2) collapsed (valid while u>=1 and w>=2)."""
    if u >= 1 and w >= 2:
        t = min(u, w // 2)
        return (u - t, v + 2 * t, w - 2 * t)
    return rule(u, v, w)


# ------------------------------------------------------------------ TM
TM = '1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA'


def parse(code):
    tab = {}
    for si, row in enumerate(code.split('_')):
        for sym in (0, 1):
            t = row[3 * sym:3 * sym + 3]
            tab[(si, sym)] = None if t == '---' else (int(t[0]), 1 if t[1] == 'R' else -1, 'ABCDEF'.index(t[2]))
    return tab


TAB = parse(TM)


def tape_of(u, v, w):
    left = '1' * 5 + '10' * u + '1'
    right = '11' + '01' * v + '1' * w + '01'
    return left, right


def tm_check(u, v, w, max_steps=10 ** 6):
    """Build J(u,v,w) on a tape, run the TM until the configuration predicted by `rule` appears (state A,
    head on the first symbol of the right part) or HALT.  Return True if the prediction is reached first
    among all J-shaped configurations, i.e. the TM reaches exactly the predicted configuration."""
    pred = rule(u, v, w)
    left, right = tape_of(u, v, w)
    N = 4 * (len(left) + len(right)) + 200
    tape = bytearray(N)
    off = N // 4
    for i, ch in enumerate(left):
        tape[off + i] = 1 if ch == '1' else 0
    head = off + len(left)
    for i, ch in enumerate(right):
        tape[head + i] = 1 if ch == '1' else 0
    state = 0
    if pred == HALT:
        target = None
    else:
        L2, R2 = tape_of(*pred)
        target = (bytes(1 if c == '1' else 0 for c in L2), bytes(1 if c == '1' else 0 for c in R2))
    for step in range(1, max_steps + 1):
        sym = tape[head]
        tr = TAB[(state, sym)]
        if tr is None:
            return pred == HALT, step
        wsym, mv, ns = tr
        tape[head] = wsym
        head += mv
        state = ns
        if head <= 1 or head >= N - 2:
            return None, step
        if target is not None and state == 0:
            L2, R2 = target
            if tape[head - len(L2):head] == L2 and tape[head:head + len(R2)] == R2:
                # rest of the tape blank?
                if not any(tape[:head - len(L2)]) and not any(tape[head + len(R2):]):
                    return True, step
    return False, max_steps


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    bad = []
    tested = 0
    cases = {}
    # every branch: u in 0..9, v in 0..9, w in 0..9 exhaustively, plus random larger
    pts = [(u, v, w) for u in range(10) for v in range(10) for w in range(10)]
    pts += [(rng.randrange(0, 60), rng.randrange(0, 60), rng.randrange(0, 60)) for _ in range(600)]
    for (u, v, w) in pts:
        ok, st = tm_check(u, v, w)
        tested += 1
        if ok is not True:
            bad.append(((u, v, w), rule(u, v, w), ok, st))
    print('A6 TM check of single rules: tested', tested, 'bad', len(bad), bad[:10])
