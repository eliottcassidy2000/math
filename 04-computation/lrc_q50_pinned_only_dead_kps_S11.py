#!/usr/bin/env python3
"""
lrc_q50_pinned_only_dead_kps_S11.py  --  kind-pasteur-2026-07-05-S11, HYP-4137.

Companion to lrc_q50_refutation_kps_S11.py.  Shows the "pinned-only repair" of the
Q50 conjecture is ALSO dead: Q0 = infinity.  A single runner == 0 (mod L) with
L = lcm(2..25) is == 0 (mod q) for every PINNED-ONLY modulus q | L, hence always in
the danger band -> blocks every pinned-only witness.  Since an L-runner makes the
profile's pinning (F5) and covering vacuous, the other runners are unconstrained, so
the loose-branch witness modulus is UNBOUNDED.  Verdict: no fixed-modulus template
(bound-50, pinned-only q|L, or any fixed set) closes the loose branch; only the
real-valued TightLooseDichotomy survives.
"""
from math import gcd
from functools import reduce
import random
random.seed(2026)

def lcm(a, b): return a * b // gcd(a, b)
L = reduce(lcm, range(2, 26))
def mu(s): return (2 * s + 24) // 25
def distZ(x, q):
    r = x % q; return min(r, q - r)
def is_pinned_only(q):
    n = q; f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return all(p**e <= 25 for p, e in f.items())
def has_template_witness(base, s):
    m = mu(s)
    for k in range(1, s):
        if all(m <= (v * k) % s <= s - m for v in base): return k
    return None
def witness_at(vs, q):
    if any(v % q == 0 for v in vs): return None
    m = mu(q)
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        if all(distZ(v * a, q) >= m for v in vs): return a
    return None
def covering_2_14(v): return all(any(x % q == 0 for x in v) for q in range(2, 15))
def compressed(v):
    return all(any(j != i and abs(v[i]) <= 13 * abs(v[j]) for j in range(len(v)))
               for i in range(len(v)))
def primitive(v): return reduce(gcd, [abs(x) for x in v]) == 1
def is_tight(base):
    a = [abs(x) for x in base]; g = reduce(gcd, a)
    for c in range(2, g + 1):
        if g % c: continue
        if all(x % c == 0 and 1 <= x // c <= 12 for x in a): return c
    return None

print("=" * 66)
print("PART A -- the composite-blocking family: Q0 = infinity")
print("=" * 66)
smalls = [1, 3, 4, 5, 7, 8, 9, 11, 13, 17, 19]
base = [L] + smalls
istar = 2 * L
V = base + [istar]
print(f"L = lcm(2..25) = {L}")
print(f"base = {{L, {smalls}}},  istar = 2L")
for name, ok in [("CoveringFamily(2..14)", covering_2_14(V)),
                 ("compressed", compressed(V)),
                 ("primitive", primitive(V)),
                 ("istar=2L argmax", all(abs(x) <= abs(istar) for x in V)),
                 ("base NOT tight", is_tight(base) is None)]:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
po = [q for q in range(26, 3001) if is_pinned_only(q) and has_template_witness(base, q) is not None]
print(f"  pinned-only witnesses q|L (scan <=3000): {po}   => Q0 = INFINITY")
fw = next((q for q in range(26, 200)
           if not is_pinned_only(q) and has_template_witness(base, q) is not None), None)
print(f"  smallest FREE witness: q = {fw}  (loose only via a height-dependent witness)")

print("\n" + "=" * 66)
print("PART B -- high-height shapes push the pinned-only bound past the census 69")
print("=" * 66)
PIN_Q = [q for q in range(26, 601) if is_pinned_only(q)]
def Q0(vs):
    for q in PIN_Q:
        if witness_at(vs, q) is not None: return q
    return 9999
def pin_le25(vs):
    for q in range(2, 26):
        if any(v % q == 0 for v in vs): continue
        for b in range(1, q // 2 + 1):
            if gcd(b, q) != 1: continue
            if not any(v % q in (b, q - b) for v in vs): return False
    return True
def ok(vs): return len(set(vs)) == 12 and covering_2_14(vs) and pin_le25(vs) and reduce(gcd, vs) == 1
best = 69; bf = None
cur = [1, 3, 4, 5, 7, 8, 9, 11, 12, 19, 23, 30]
for _ in range(120000):
    cand = list(cur); cand[random.randrange(12)] = random.randint(1, 5000)
    cand = sorted(set(cand))
    if len(cand) != 12 or not ok(cand): continue
    if Q0(cand) >= (Q0(cur) if ok(cur) else -1):
        cur = cand
        if Q0(cand) < 9999 and Q0(cand) > best:
            best = Q0(cand); bf = list(cand)
print(f"census-bounded max pinned-only bound = 69; high-height hill-climb reaches Q0 = {best}")
if bf:
    comp = [v for v in bf if sum(v % p == 0 for p in [2,3,5,7,11,13,17,19,23]) >= 3]
    print(f"  attained by shape with highly-composite (blocking) runners: {comp}")

print("\nVERDICT: Q0 = inf (Part A) and grows with height (Part B). No fixed-modulus")
print("template closes the loose branch; use the real-valued TightLooseDichotomy.")
