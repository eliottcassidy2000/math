#!/usr/bin/env python3
"""Audit (b2) of S35: the 2-adic claims about b(q), from the independently built newform table (audit_b_newform_bq.json),
compared with the session's table paley_k6p6_cq_table_20261014.json (c_q := sqrt(M_q + (q+1)^2), from the K6-P6 sum).

Checks:
  1. |b(q)| = c_q at every prime of the session's table (this upgrades c_q from 'perfect-square root' to 'newform coefficient').
  2. depth table: max(v2(c - q - 1), v2(c + q + 1)) minima per class mod 16, and v2(Y_q - 1) minima.
  3. 'the value 7 is attained at every q = 11 (mod 16) tested'.
  4. valuation law v2(b(q)) >= min(v2(q+1), 7), sharpness of 7, minimal examples, zeros of b.
  5. q = 1 (mod 4): b(q) = 2 (mod 4) and v2(Y_q - 1) = 2.
  6. Mersenne primes, q = 191, 383.
Also extends 2-4 to all primes q = 3 (mod 4) below 524,288 (newform side only).
"""
import json, math
from collections import defaultdict

SES = '../../../05-knowledge/results/paley_k6p6_cq_table_20261014.json'


def v2(n):
    n = abs(n)
    return 99 if n == 0 else (n & -n).bit_length() - 1


b = {int(k): tuple(v) for k, v in json.load(open('audit_b_newform_bq.json')).items()}
ses = {int(k): v for k, v in json.load(open(SES)).items()}


def absb(q):
    re, im = b[q]
    s = re * re + im * im
    r = math.isqrt(s)
    assert r * r == s
    return r


mism = [q for q in ses if absb(q) != ses[q]]
print(f'1. session table: {len(ses)} primes q = 3 mod 4 in [{min(ses)}, {max(ses)}]; |b(q)| (newform) = c_q (K6-P6 sum) at all: {not mism}; mismatches: {mism[:10]}')

for label, QMAX in (('session range (q < 200,000)', 200000), ('extended (q < 524,288)', 524288)):
    mind = defaultdict(lambda: 99); minv = defaultdict(lambda: 99); argm = {}
    exact7_11 = []; count = defaultdict(int)
    byS = defaultdict(list)
    for q in sorted(b):
        if q > QMAX or q % 4 != 3 or q < 7:
            continue
        c = absb(q)
        r = q % 16
        d = max(v2(c - q - 1), v2(c + q + 1))
        M = c * c - (q + 1) ** 2            # = Y_q - 1 at Paley primes, under the square law
        count[r] += 1
        if d < mind[r]:
            mind[r] = d; argm[r] = q
        minv[r] = min(minv[r], v2(M))
        if r == 11:
            exact7_11.append(v2(M) == 7)
        byS[v2(q + 1)].append((v2(c), q, c))
    print(f'2. {label}: per class mod 16: ' + '; '.join(
        f'q = {r} mod 16: {count[r]} primes, min depth {mind[r]} (first at q = {argm[r]}), min v2(Y-1) {minv[r]}' for r in (3, 7, 11, 15)))
    print(f'3. {label}: v2(Y_q - 1) = 7 exactly at every q = 11 mod 16: {all(exact7_11)} ({sum(exact7_11)}/{len(exact7_11)})')
    K = 99
    line = []
    for s in sorted(byS):
        lst = byS[s]
        mn = min(v for v, _, _ in lst)
        ex = min(lst)
        line.append(f's={s}:{len(lst)}p,min v2(b)={mn if mn < 99 else "inf"}(q={ex[1]})')
        if mn < s:
            K = min(K, mn)
        assert mn >= min(s, 7) or True
    viol = [(q, c) for s in byS for v, q, c in byS[s] if v < min(s, 7)]
    print(f'4. {label}: valuation law v2(b) >= min(v2(q+1), 7) violations: {viol[:5]}; largest working K = {K}')
    print('   ' + '  '.join(line))
    zeros = sorted(q for q in b if q <= QMAX and b[q] == (0, 0))
    print(f'   zeros of b(q) at primes q <= {QMAX}: {zeros}')

# minimal examples of v2(b) = v2(q+1) for s = 2..7 and v2(b) = 7 for s = 8..14
ex1 = {}
ex2 = {}
for q in sorted(b):
    if q % 4 != 3 or q < 7 or q > 200000:
        continue
    s, vb = v2(q + 1), v2(absb(q))
    if 2 <= s <= 7 and vb == s and s not in ex1:
        ex1[s] = q
    if 8 <= s <= 14 and vb == 7 and s not in ex2:
        ex2[s] = q
print('   minimal q with v2(b) = v2(q+1), s = 2..7:', [ex1.get(s) for s in range(2, 8)], ' (S35: 11, 7, 47, 31, 191, 383)')
print('   minimal q with v2(b) = 7, s = 8..14:', [ex2.get(s) for s in range(8, 15)], ' (S35: 6911, 23039, 21503, 6143, 20479, 73727, 81919)')

# 5. q = 1 mod 4
bad = [q for q in b if q % 4 == 1 and (b[q][0] % 4 != 2 or b[q][1] != 0)]
badY = []
for q in b:
    if q % 4 == 1:
        c = absb(q)
        Y = c * c - q * q + 2 * q
        if v2(Y - 1) != 2:
            badY.append(q)
print(f'5. q = 1 mod 4 (all primes < 524,288): b(q) = 2 mod 4 fails at {bad[:5]}; v2(Y_q - 1) = 2 fails at {badY[:5]}')

# 6. Mersenne and special values
for p in (3, 5, 7, 13, 17, 19):
    q = 2 ** p - 1
    c = absb(q)
    M = c * c - (q + 1) ** 2
    print(f'6. Mersenne p = {p}: q = {q}: b(q) = {b[q]}, |b(q)| = {c}, M_q = Y_q - 1 = {M}, v2 = {v2(M)}')
for q in (3, 7, 31, 127, 191, 199, 223, 383):
    c = absb(q)
    print(f'   q = {q}: |b(q)| = {c}, Y_q - 1 (Paley) = {c * c - (q + 1) ** 2}')
Y1 = [q for q in sorted(b) if q % 4 == 3 and q < 700 and absb(q) == q + 1]
print('   primes q = 3 mod 4 below 700 with |b(q)| = q + 1:', Y1)
# Sato-Tate moment for the 124 primes 3..691 (both classes)
import statistics
vals = [(absb(q) / (2 * q)) ** 2 for q in sorted(b) if q <= 691]
print(f'   mean of (|b(q)|/2q)^2 over {len(vals)} primes 3..691: {statistics.mean(vals):.4f} (se {statistics.stdev(vals) / len(vals) ** 0.5:.4f})')
