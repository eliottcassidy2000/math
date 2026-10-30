#!/usr/bin/env python3
# Verification gates are explicit and remain active under python -O.
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

def v(n):
    if not n:
        raise AssertionError()
    k = 0
    while n % 5 == 0:
        k += 1
        n //= 5
    return k

def offsets(P, Q, h):
    ds = []
    d = 0
    T = 1
    for e in range(1, h + 1):
        q = 5 ** e
        m = P * pow(Q, -1, q) % q
        c = 3 * pow(2, -1, q) % q
        candidates = range(2) if e == 1 else range(5)
        good = [d + k * T for k in candidates if m * pow(c, d + k * T, q) % q == 1]
        if not len(good) == 1:
            raise AssertionError()
        d = good[0]
        ds.append(d)
        T = 2 * 5 ** (e - 1)
    return ds

def audit(w, name):
    P = Q = 1
    B = 0
    actual = 0
    ep = []
    for a in w:
        P, Q, B = (3 * P, 2 ** a * Q, 3 * B + Q)
        actual += v(gcd(P - Q, B))
        if a > 1:
            ep.append([P, Q, v(B + Q - P), 0])
        else:
            ep[-1][3] += 1
    if not all((p * pow(q, -1, 5) % 5 in (1, 4) for p, q, h, L in ep)):
        raise AssertionError()
    excess = sum((h + F(5, 8) * (1 - F(1, 5 ** h)) * L for p, q, h, L in ep), F()) - F(5 * sum(w), 8)
    trunc = {}
    for E in [8, 9, 10]:
        total = sum((min(h, E) + F(5, 8) * (1 - F(1, 5 ** min(h, E))) * L for p, q, h, L in ep), F())
        trunc[E] = str(total - F(5, 8) * (1 - F(1, 5 ** E)) * sum(w))
    refined = {}
    dsall = []
    for p, q, h, L in ep:
        ds = offsets(p, q, max(h, 2))
        dsall.append(ds)
        if not all((ds[e] >= ds[e - 1] and (ds[e] - ds[e - 1]) % (2 * 5 ** (e - 1)) == 0 for e in range(1, len(ds)))):
            raise AssertionError()
        exact = sum((max(0, 1 + (L - d) // (2 * 5 ** e)) for e, d in enumerate(ds[:h])))
        if not exact == sum((min(h, v(p * 3 ** j - q * 2 ** j)) for j in range(L + 1))):
            raise AssertionError()
    for r in [1, 2]:
        total = F()
        for (p, q, h, L), ds in zip(ep, dsall):
            low = sum((max(0, 1 + (L - ds[e - 1]) // (2 * 5 ** (e - 1))) for e in range(1, min(h, r) + 1)))
            upper = F(low)
            if h > r and L >= ds[r - 1]:
                upper += h - r + F(5, 8) * (F(1, 5 ** r) - F(1, 5 ** h)) * (L - ds[r - 1])
            exact = sum((max(0, 1 + (L - d) // (2 * 5 ** e)) for e, d in enumerate(ds[:h])))
            if not upper >= exact:
                raise AssertionError()
            pending = True
            online = F()
            for j in range(L + 1):
                mm = p * 3 ** j * pow(q * 2 ** j, -1, 5 ** r) % 5 ** r
                online += sum((mm % 5 ** e == 1 for e in range(1, min(h, r) + 1)))
                if h > r:
                    if pending and mm == 1:
                        online += h - r
                        pending = False
                    elif not pending:
                        online += F(5, 8) * (F(1, 5 ** r) - F(1, 5 ** h))
            if not online == upper:
                raise AssertionError()
            total += upper
        refined[r] = str(total)
        if not actual <= total < F(5 * sum(w), 8):
            raise AssertionError()
    y = B * pow(Q, -1, P) % P
    if y % 2 == 0:
        y += P
    while Q * y <= B:
        y += 2 * P
    x = (Q * y - B) // P
    if not (x > 0 and x % 2 == 1):
        raise AssertionError()
    n = x
    for a in w:
        nn = 3 * n + 1
        if not (nn % 2 ** a == 0 and nn // 2 ** a % 2 == 1):
            raise AssertionError()
        n = nn // 2 ** a
    if not n == y:
        raise AssertionError()
    return dict(name=name, length=len(w), clock=sum(w), heights=[e[2] for e in ep], runs=[e[3] for e in ep], actual=actual, excess=str(excess), truncated=trunc, refined=refined, first_offsets=dsall, positive_source=x, positive_endpoint=y)
R = [5] + [1] * 8 + [5]
W = R + [1] * 233 + R
U5 = W + [1] * 8 + W
V = [3, 1, 3] + [1] * 6 + [3, 1, 3]
U3 = V + [1] * 231 + V
rs = [audit(U3, 'U3'), audit(U5, 'U5')]
if not (rs[0]['actual'] == 154 and rs[0]['excess'] == '769/5000' and (rs[0]['truncated'][9] == '60112/390625')):
    raise AssertionError()
if not (rs[1]['actual'] == 317 and rs[1]['excess'] == '479371/312500' and (rs[1]['truncated'][10] == '1043012/1953125')):
    raise AssertionError()
if not (rs[1]['refined'][1] == '26448046/78125' and rs[1]['refined'][2] == '5094453/15625'):
    raise AssertionError()
if not rs[0]['refined'][1] == '198661/1250':
    raise AssertionError()
Path(__file__).with_name('new_boundaries_independent.json').write_text(json.dumps(rs, indent=2) + '\n')
print('PASS independent U3/U5, exact content and budgets, lifted discrete logs, epoch counts, pending-bit machine, reverse-CRT positive-source witnesses')
checks = 0
for m in range(1, 125):
    if m % 5 == 0:
        continue
    active = m % 5 in (1, 4)
    ds = offsets(m, 1, 3) if active else []
    for h in range(4):
        for r in [1, 2]:
            for L in range(31):
                actual = 0
                mm = m
                for j in range(L + 1):
                    actual += sum((mm % 5 ** e == 1 for e in range(1, h + 1)))
                    mm = mm * 3 * pow(2, -1, 125) % 125
                upper = F()
                if active:
                    upper = sum((max(0, 1 + (L - ds[e - 1]) // (2 * 5 ** (e - 1))) for e in range(1, min(h, r) + 1)))
                    if h > r and L >= ds[r - 1]:
                        upper += h - r + F(5, 8) * (F(1, 5 ** r) - F(1, 5 ** h)) * (L - ds[r - 1])
                if not actual <= upper:
                    raise AssertionError((m, h, r, L, actual, upper))
                checks += 1
print('PASS', checks, 'general multiplier/height/precision/length epoch bounds')
