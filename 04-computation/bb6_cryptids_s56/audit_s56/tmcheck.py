#!/usr/bin/env python3
"""Referee: generic direct-TM checker for a configuration family.

Given a TM, a family encoder enc(params) -> (state_index, left_string, right_string) [head on right[0]],
a parser parse(state, left, right) -> params or None, and a rule function rule(params) -> params | 'HALT',
build the tape for random params, run the TM and check that the FIRST family-shaped configuration reached
(after >= 1 TM step) is exactly rule(params), or that the TM halts when rule says HALT (and reaches no
family configuration before halting).
"""
import re

HALT = 'HALT'


def parse_tm(code):
    tab = {}
    for si, row in enumerate(code.split('_')):
        for sym in (0, 1):
            t = row[3 * sym:3 * sym + 3]
            tab[(si, sym)] = None if t == '---' else (int(t[0]), 1 if t[1] == 'R' else -1, 'ABCDEF'.index(t[2]))
    return tab


def check(code, enc, parse, rule, params, max_steps=2 * 10 ** 6, state_filter=None, norm=lambda x: x, max_pass=4):
    passed = []
    tab = parse_tm(code)
    st0, L, R = enc(params)
    pred = rule(params)
    N = 6 * (len(L) + len(R)) + 400
    tape = bytearray(N)
    off = N // 3
    for i, ch in enumerate(L):
        tape[off + i] = 1 if ch == '1' else 0
    head = off + len(L)
    for i, ch in enumerate(R):
        tape[head + i] = 1 if ch == '1' else 0
    state = st0
    lo, hi = off, head + len(R)
    for step in range(1, max_steps + 1):
        tr = tab[(state, tape[head])]
        if tr is None:
            return ('ok' if pred == HALT else 'halted-unexpectedly'), step, None
        w, mv, ns = tr
        tape[head] = w
        head += mv
        state = ns
        if head < lo:
            lo = head
        if head > hi:
            hi = head
        if head <= 2 or head >= N - 3:
            return 'tape-overflow', step, None
        if state_filter is not None and state != state_filter:
            continue
        # extract nonblank window
        a = lo
        while a < head and tape[a] == 0:
            a += 1
        b = hi
        while b >= head and tape[b] == 0:
            b -= 1
        left = ''.join('1' if x else '0' for x in tape[a:head])
        right = ''.join('1' if x else '0' for x in tape[head:b + 1])
        got = parse(state, left, right)
        if got is not None:
            if pred == HALT:
                return 'family-config-before-halt', step, got
            if norm(tuple(got)) == norm(tuple(pred)):
                return ('ok' if not passed else 'ok-via-%d' % len(passed)), step, passed[:3]
            passed.append(tuple(got))
            if len(passed) > max_pass:
                return 'mismatch', step, passed[:3]
    return 'timeout', max_steps, passed[:3]
