#!/usr/bin/env python3
"""An exact counterexample to the epoch budget, not to the content bound."""
from fractions import Fraction
from pathlib import Path
import json, math
V = (3, 1, 3) + (1,) * 6 + (3, 1, 3)
w = V + (1,) * 231 + V

def val(n):
    if not n != 0:
        raise RuntimeError('verification failed at original line 10')
    e = 0
    while n % 5 == 0:
        e += 1
        n //= 5
    return e
P = Q = 1
B = 0
content = 0
epochs = []
for a in w:
    P, Q, B = (3 * P, 2 ** a * Q, 3 * B + Q)
    h = val(B + Q - P)
    active = P * pow(Q, -1, 5) % 5 in (1, 4)
    content += val(math.gcd(abs(P - Q), B))
    if not active:
        raise RuntimeError('verification failed at original line 20')
    if a > 1:
        epochs.append([a, h, 0])
    else:
        epochs[-1][2] += 1
if not epochs == [[3, 0, 1], [3, 2, 6], [3, 0, 1], [3, 5, 231], [3, 0, 1], [3, 2, 6], [3, 0, 1], [3, 9, 0]]:
    raise RuntimeError('verification failed at original line 23')
excess = sum((Fraction(h) - Fraction(5 * a, 8) - Fraction(5 * L, 8 * 5 ** h) for a, h, L in epochs), Fraction())
E = 9
q = 5 ** E
D = 8 * q
N = 5 * q - 5
scaled = sum((D * min(h, E) - N * a - L * (5 * q // 5 ** min(h, E) - 5) for a, h, L in epochs))
if not (scaled == 2404480 and excess > 0 and (content == 154)):
    raise RuntimeError('verification failed at original line 28')
if not content < Fraction(N * sum(w), D) < Fraction(5 * sum(w), 8):
    raise RuntimeError('verification failed at original line 29')
modulus = 2 ** (sum(w) + 1)
source = (2 ** sum(w) - B) * pow(P, -1, modulus) % modulus
if not (source > 0 and source % 2 == 1):
    raise RuntimeError('verification failed at original line 32')
y = source
for a in w:
    t = 3 * y + 1
    if not (t % 2 ** a == 0 and t // 2 ** a % 2 == 1):
        raise RuntimeError('verification failed at original line 36')
    y = t // 2 ** a
report = dict(positive_source=source, positive_endpoint=y, word='V 1^231 V, V=(3,1,3,1^6,3,1,3)', length=len(w), clock=sum(w), epochs_a_h_L=epochs, actual_untruncated_content=content, original_budget_excess=str(excess), E9_envelope_excess=str(Fraction(scaled, D)), E9_content_bound=str(Fraction(N * sum(w), D)), status='REFUTES universal epoch budget, including tighter h-only version; NOT a counterexample to sharp content bound or Collatz')
p = Path(__file__).with_name('envelope_obstruction.json')
p.write_text(json.dumps(report, indent=2) + '\n')
print(p.read_text())
