#!/usr/bin/env python3
"""Referee: A1 rule transcription (paper A.1) checked against the TM; then the session's generic engine
(copied, read-only) vs plain stepping with the referee's rules on the blank-tape orbit."""
import sys, os, re, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), 'session_copy'))
from tmcheck import check, HALT
A1 = '1RB1LD_1RC0RA_1LA1RE_0LA0LB_0RB1RF_---1RD'
def enc(p):
    a, b, c, d = p
    return 0, '1' * a + '10' * b, '01' * c + '1' * d
RL = re.compile(r'^(1*)((?:10)*)$'); RR = re.compile(r'^((?:01)*)(1*)$')
def parse(state, left, right):
    if state != 0: return None
    m1 = RL.match(left); m2 = RR.match(right)
    if not m1 or not m2: return None
    return (len(m1.group(1)), len(m1.group(2)) // 2, len(m2.group(1)) // 2, len(m2.group(2)))
def rule(p):
    a, b, c, d = p
    if c > 0: return (a, b + 1, c - 1, d)
    if d == 0: return (a + 1, 0, 0, 2 * b + 2)
    if a >= 3: return (a - 3, b + 2, 0, d - 1)
    if a == 2: return (0, 0, b + 1, d - 1)
    if a == 1: return (1, 1, b + 1, d - 1) if b > 0 else (0, 0, 0, d + 2)
    if d == 1: return (0, 0, 2, 5) if b == 0 else (0, 0, 2, 2 * b + 1)
    if d >= 2 and b >= d - 1: return (3 * d - 4, 1, 2, 2 * b - 2 * d + 3)
    if d == b + 2: return HALT
    if d == b + 3: return (3 * b + 2, 0, 0, 6)
    if d == b + 4: return (3 * b + 4, 0, 0, 4)
    return (3 * b + 3, 2, 0, d - b - 5)
rng = random.Random(9)
pts = [(a, b, c, d) for a in range(6) for b in range(6) for c in range(4) for d in range(8)]
pts += [tuple(rng.randrange(0, 30) for _ in range(4)) for _ in range(300)]
res = {}
bad = []
for p in pts:
    r, st, got = check(A1, enc, parse, rule, p, max_steps=10 ** 6)
    res[r] = res.get(r, 0) + 1
    if not r.startswith('ok'): bad.append((p, rule(p), r, got))
print('A1 TM check:', res, bad[:5], flush=True)
# engine vs plain
import rs_accel_20261025 as R
eng = R.Engine('A1')
s = (0, 0, 0, 0)
while not R.ANCHOR['A1'](s) or s == (0, 0, 0, 0):
    s = rule(s)
    if s == (0, 0, 0, 0) or s == HALT: break
# plain anchors from the same start
def plain_anchors(s0, n_anchor, max_steps):
    out = []; s = s0
    for i in range(max_steps):
        s = rule(s)
        if s == HALT: out.append('HALT'); break
        if s[2] == 0 and s[3] == 0:
            out.append(s)
            if len(out) >= n_anchor: break
    return out, i
t0 = time.time()
pa, steps = plain_anchors(s, 10 ** 9, 3 * 10 ** 7)
print('plain: anchors', len(pa), 'in', steps + 1, 'rule steps', f'({time.time() - t0:.0f}s)', flush=True)
# engine segments
seg_states = []
x = s
for k in range(len(pa)):
    x2, toks = eng.segment(x)
    seg_states.append(x2)
    x = x2
    if len(seg_states) >= len(pa): break
pset = {q: i for i, q in enumerate(pa)}
idx = [pset.get(q) for q in seg_states]
mono = all(i is not None for i in idx) and all(idx[i] < idx[i + 1] for i in range(len(idx) - 1))
print('engine segment ends are plain anchors in increasing order:', mono, '| engine segments', len(seg_states),
      '| anchors skipped by jumps:', (idx[-1] + 1 - len(seg_states)) if mono else None, '| engine l1 jumps', eng.stats)
