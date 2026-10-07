#!/usr/bin/env python3
"""Audit A, item 2: Lemma 1 (a)-(d) and the one-step drift of V = |f|^theta s^h, exhaustively on a box and on
random large states. Exact rationals for sizes; theta = 1/2 for the drift."""
import random, math
from fractions import Fraction
from a1_chain_vs_direct import chain_step

def v2(x): return (x & -x).bit_length() - 1 if x else 10**9
def m_h(h): return v2(3 ** h - 1)
def fsize(k, N):  # |f| = |e| 3^-max(k,0), e = N / 3^max(0,-k)
    return abs(Fraction(N, 3 ** max(0, -k))) / 3 ** max(k, 0)
def coin(k, N, beta): return beta if k >= 0 else beta ^ (N & 1)

th = 0.5
s = 2 ** th * (1 - math.sqrt(1 - 0.75 ** th)); rho = 2 ** -th * s
def V(k, N): return float(fsize(k, N)) ** th * s ** abs(k)
def eps(h): return 2 ** -th * s ** (h - 1)

fails = {'size': 0, 'dep': 0, 'dir': 0, 'k0run': 0, 'runval': 0, 'flipdrift': 0, 'rundrift': 0, 'k0drift': 0,
         'away_add_1/18': 0, 'fair': 0}
cnt = 0; worst = {'flip': -9, 'run': -9}
def check_state(k, N):
    global cnt
    if (k, N) == (0, 0): return
    h = abs(k); f = fsize(k, N); sig = N & 1
    outs = {}
    for beta in (0, 1):
        c = coin(k, N, beta)
        k2, N2 = chain_step(k, N, beta)
        outs[c] = (k2, N2)
        f2 = fsize(k2, N2); cnt += 1
        if k == 0 and sig == 1:
            if f2 > f / 2 + Fraction(1, 6): fails['dep'] += 1
        else:
            A = Fraction(1, 2) if c == 0 else Fraction(3, 2)
            if f2 > A * f + Fraction(1, 2): fails['size'] += 1
        if sig == 1 and h >= 1:
            toward = abs(k2) == h - 1
            if toward != (c == 1): fails['dir'] += 1
            if c == 0 and f2 > f / 2 + Fraction(1, 18): fails['away_add_1/18'] += 1
        if k == 0 and sig == 0:
            if Fraction(N2) != (Fraction(N, 2) if c == 0 else Fraction(3 * N, 2)): fails['k0run'] += 1
        if sig == 0 and h >= 1:
            w = v2(N); w2 = v2(N2); m = m_h(h)
            if N == 0:
                ok = (w2 >= 10**9) if c == 0 else (w2 == m - 1)
            elif w < m: ok = (w2 == w - 1)
            elif w > m: ok = (w2 == w - 1) if c == 0 else (w2 == m - 1)
            else: ok = (w2 == m - 1) if c == 0 else (w2 >= m)
            if not ok: fails['runval'] += 1
    if set(outs) != {0, 1}: fails['fair'] += 1
    EV = 0.5 * (V(*outs[0]) + V(*outs[1]))
    if h >= 1:
        slack = EV - V(k, N) - eps(h)
        if sig == 1:
            worst['flip'] = max(worst['flip'], slack / (1 + V(k, N)))
            if slack > 1e-9 * (1 + V(k, N)): fails['flipdrift'] += 1
        else:
            slack2 = EV - 0.5 * (0.5 ** th + 1.5 ** th) * V(k, N) - 0.5 * eps(h)
            worst['run'] = max(worst['run'], slack2 / (1 + V(k, N)))
            if slack2 > 1e-9 * (1 + V(k, N)): fails['rundrift'] += 1
    elif sig == 0:
        if EV > V(k, N) * 0.5 * (0.5 ** th + 1.5 ** th) * (1 + 1e-9): fails['k0drift'] += 1

for k in range(-12, 13):
    for N in range(-3000, 3001):
        check_state(k, N)
rng = random.Random(7)
for _ in range(200000):
    k = rng.randint(-60, 60); N = rng.getrandbits(rng.choice([8, 40, 200])) * rng.choice([1, -1])
    if rng.random() < 0.3: N <<= rng.randint(1, 12)
    if rng.random() < 0.2 and k != 0:  # states at the edge of the run regime
        h = abs(k); N = (3 ** h - 1) * rng.getrandbits(30) * rng.choice([1, -1])
    check_state(k, N)
print(f"[lemmas] (state, bit) pairs checked: {cnt}; failures: {fails}")
print(f"[lemmas] worst slack of E[V'] - V - eps_h at flips: {worst['flip']:.3e}; at runs vs 0.966V+eps/2: {worst['run']:.3e}")
print(f"[lemmas] s = {s:.6f}, rho = {rho:.6f}, run bracket = {0.5*(0.5**th+1.5**th):.6f}")
