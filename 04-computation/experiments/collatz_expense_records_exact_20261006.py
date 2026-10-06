#!/usr/bin/env python3
"""Exact recomputation of the S9 expense records by isolating log2 3 in a rational bracket (the packing proof's
root-isolation idea), as the input to the Lean certificate
04-computation/lean/standalone/collatz_expense_records_certificate_20261006.lean.

S9 (collatz_expense_diophantine_20261005.md): for a first-descent segment type (l, A) with 2^A > 3^l,
q(l, A) = ceil(l / (A - l theta)), theta = log2 3.  Worst type at length l: A = A_min(l) (least A with 2^A > 3^l);
q_max(l) = q(l, A_min(l)); records of q_max sit at the upper best approximations of theta.

Certificate shape (all integer arithmetic):
  (B) bracket: 3^QLO < 2^PLO is false and ... precisely  2^PLO < 3^QLO  (so PLO/QLO < theta)  and  3^QHI < 2^PHI
      (so theta < PHI/QHI), with PLO/QLO, PHI/QHI consecutive convergents of theta;
  (A) for each l: 2^(A-1) < 3^l < 2^A;
  (Q) f(theta) = l/(A - l theta) is increasing in theta, so f(PLO/QLO) < f(theta) < f(PHI/QHI); with
      f(p/q') = l q' / (A q' - l p), the value q = ceil f(theta) is certified by
      (q - 1) (A QLO - l PLO) <= l QLO   and   l QHI <= q (A QHI - l PHI)   (strict since theta is irrational).
"""
import sys, json, time
from fractions import Fraction as Fr
import mpmath as mp
import flint

LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 320
mp.mp.dps = 80
theta = mp.log(3) / mp.log(2)

# continued fraction convergents of theta
cf = []
x = theta
for _ in range(25):
    a = int(mp.floor(x))
    cf.append(a)
    x = 1 / (x - a)
convs = []
h0, h1, k0, k1 = 0, 1, 1, 0
for a in cf:
    h0, h1 = h1, a * h1 + h0
    k0, k1 = k1, a * k1 + k0
    convs.append((h1, k1))
print('continued fraction of log2 3:', cf[:16])
print('convergents:', convs[:16])
# pick consecutive convergents bracketing theta with denominators around 1e5..1e7
for i in range(len(convs) - 1):
    (p1, q1), (p2, q2) = convs[i], convs[i + 1]
    if q2 > 10 ** 6:
        break
lo, hi = sorted([(p1, q1), (p2, q2)], key=lambda pq: Fr(pq[0], pq[1]))
PLO, QLO = lo
PHI, QHI = hi
t0 = time.time()
two, three = flint.fmpz(2), flint.fmpz(3)
okB = (two ** PLO < three ** QLO) and (three ** QHI < two ** PHI)
print(f'bracket {PLO}/{QLO} < log2 3 < {PHI}/{QHI}: certified by exact powers = {okB} ({time.time() - t0:.1f}s); width {float(Fr(PHI, QHI) - Fr(PLO, QLO)):.3e}')
assert okB


def a_min(l):
    t = 3 ** l
    A = t.bit_length()
    assert 2 ** (A - 1) < t < 2 ** A
    return A


rows, records, best, failures = [], [], 0, []
for l in range(1, LMAX + 1):
    A = a_min(l)
    dlo = A * QLO - l * PLO      # > 0
    dhi = A * QHI - l * PHI      # > 0 required (A/l > PHI/QHI)
    if dhi <= 0:
        failures.append((l, A, 'A/l not above the upper bracket end'))
        continue
    flo = Fr(l * QLO, dlo)
    fhi = Fr(l * QHI, dhi)
    qlo = -(-flo.numerator // flo.denominator)   # ceil
    qhi = -(-fhi.numerator // fhi.denominator)
    if qlo != qhi:
        failures.append((l, A, f'bracket too wide: ceil range {qlo}..{qhi}'))
        continue
    q = qlo
    # the two certificate inequalities
    assert (q - 1) * dlo <= l * QLO and l * QHI <= q * dhi
    rows.append((l, A, q))
    if q > best:
        best = q
        records.append((l, A, q))
print('failures (bracket insufficient):', failures if failures else 'none')
print(f'records of q_max(l), l <= {LMAX}:')
for l, A, q in records:
    print(f'   l = {l:4d}  A = {A:4d}  ({A}/{l})  q_max = {q}')
json.dump(dict(bracket=[PLO, QLO, PHI, QHI], records=records, rows=rows, lmax=LMAX),
          open('collatz_expense_records_exact_20261006.json', 'w'))
