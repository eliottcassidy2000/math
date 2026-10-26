#!/usr/bin/env python3
"""Referee: compare the session's step() (copied, read-only) with the referee's own transcription, jump=False,
on random states including large ones (valid domains only)."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), 'session_copy'))
import tm_rules_check as M
from deepmind_cryptids_20261025 import step, HALT as SH

rng = random.Random(7)
mine = {'A6': M.a6_rule, 'A11': M.a11_rule, 'A4': M.a4_rule, 'A2': M.a2_rule, 'A5': M.a5_rule, 'A10': M.a10_rule}
for name, f in mine.items():
    bad = 0
    n = 0
    for trial in range(200000):
        big = rng.choice([10, 100, 10 ** 6, 2 ** 80])
        dim = M.FAM[name][4]
        p = tuple(rng.randrange(0, big) for _ in range(dim))
        if name == 'A11':
            p = p[:3] + (rng.choice((3, 6)),)
        if name == 'A2':
            if p[1] == 0 or (p[0] >= 1 and p[1] % 2 == 0):
                continue
        if name == 'A5':
            if p[1] == 0 or (p[0] >= 1 and p[1] % 2 == 0):
                continue
        q1 = f(p)
        if q1 is None:
            continue
        try:
            q2 = step(name, p, jump=False)
        except AssertionError:
            continue
        if q2 == SH:
            q2 = 'HALT'
        if name == 'A11':
            q1, q2 = M.a11_norm(q1), M.a11_norm(q2)
        n += 1
        if q1 != q2:
            bad += 1
            if bad <= 5:
                print(name, 'DIFF', p, 'mine', q1, 'session', q2)
    print(name, 'compared', n, 'differences', bad, flush=True)
