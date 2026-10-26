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
    r, st, got = check(A1, enc, parse, rule, p, max_steps=10 ** 6, max_pass=10 ** 6)
    res[r] = res.get(r, 0) + 1
    if not r.startswith('ok'): bad.append((p, rule(p), r, got))
print('A1 TM check:', res, bad[:5], flush=True)
# engine vs plain (collapsed loops for plain stepping: c-loop, a-loop; both trivially valid)
import rs_accel_20261025 as R
def rule_fast(p):
    a, b, c, d = p
    if c > 0: return (a, b + c, 0, d)
    if d > 0 and a >= 3:
        t = min(a // 3, d); return (a - 3 * t, b + 2 * t, 0, d - t)
    return rule(p)
s = (0, 0, 0, 0)
anch = []
t0 = time.time()
x = s
for i in range(10 ** 8):
    x = rule_fast(x)
    if x == HALT: anch.append('HALT'); break
    if x[2] == 0 and x[3] == 0:
        anch.append(x)
        if len(anch) >= 80: break
print('plain (collapsed) anchors:', len(anch), 'max bits', max(max(q).bit_length() for q in anch), f'({time.time()-t0:.0f}s)', flush=True)
eng = R.Engine('A1')
x = anch[0]
seg = []
for k in range(len(anch) - 1):
    x, toks = eng.segment(x)
    seg.append(x)
pset = {q: i for i, q in enumerate(anch)}
idx = [pset.get(q) for q in seg]
ok = all(i is not None for i in idx) and all(idx[i] < idx[i + 1] for i in range(len(idx) - 1))
print('engine segment ends coincide with plain anchors in order:', ok, '| skipped anchors:', (idx[-1] - len(seg)) if ok else idx[:10], eng.stats)
