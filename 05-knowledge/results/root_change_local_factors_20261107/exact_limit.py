"""exact_limit.py (kind-pasteur-2026-11-07-S15): the n -> infinity limit of the residue-automaton local factors.

In the letter-count model (resdp.c), prefixes with exactly two J dominate as n -> infinity, and almost all of them have
three long O/Z blocks. A long block from a start state converges, in the Cesaro sense, to the projector Pi onto the
fixed space of the block chain P (letters O and Z with probability 1/2). So the limiting class distribution of
F = D - E + U mod q is that of Pi J Pi J Pi delta_t, with t = (1,0,2,1).

Pi is applied to vectors by iterating P many times and averaging over 12 consecutive steps. That averaging removes only
periods dividing 12; the d = 0 classes have period ord_p(2) (10, 8, 18, 11 at p = 11, 17, 19, 23). Their U-phase affects
neither the value law nor the confinement event, so the class law is still exact, as the audit's lazy-chain projector
lim((I+P)/2)^k confirms to 1.2e-12 (audit m14). Convergence is reported by comparing two iteration depths. Fractions are
rationalised with denominators up to 10^7; the proved closed form is L_p(0) = 1 + (3p-2)(p-1)/(2p^3) (THM-4916).

usage: python3 exact_limit.py q [steps]
"""
import sys
import numpy as np
from fractions import Fraction

q = int(sys.argv[1]); steps = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
S = q ** 4
s = np.arange(S)
d = s % q; U = (s // q) % q; E = (s // q ** 2) % q; D = s // q ** 3
enc = lambda D, E, U, d: ((D * q + E) * q + U) * q + d
tO = enc(D, 0 * E, (2 * U) % q, d)
tZ = enc((D - E + U) % q, d, (2 * U - d) % q, d)
tJ = enc((D - E + U) % q, (D - E) % q, (D + U) % q, D)


def P(v):
    return 0.5 * (np.bincount(tO, weights=v, minlength=S) + np.bincount(tZ, weights=v, minlength=S))


def Jop(v):
    return np.bincount(tJ, weights=v, minlength=S)


def Pi(v, n):
    for _ in range(n):
        v = P(v)
    acc = np.zeros(S)
    for _ in range(12):
        acc += v; v = P(v)
    return acc / 12


def limit(n):
    v = np.zeros(S); v[enc(1, 0, 2 % q, 1)] = 1.0
    v = Pi(v, n); v = Jop(v); v = Pi(v, n); v = Jop(v); v = Pi(v, n)
    F = (D - E + U) % q
    cls = np.bincount(F, weights=v, minlength=q)
    return q * cls / cls.sum()


a = limit(steps // 2); b = limit(steps)
print('q=%d limiting local factor by class (steps %d vs %d, max diff %.2e):' % (q, steps // 2, steps, np.abs(a - b).max()))
p = min(f for f in range(2, q + 1) if q % f == 0)
for r in range(q):
    print('  class %d: %.10f ~ %s' % (r, b[r], Fraction(float(b[r])).limit_denominator(10**7)))
kq = 0; t = q
while t % p == 0:
    t //= p; kq += 1
byv = {}
for r in range(q):
    v_ = kq if r % p ** kq == 0 else next(k for k in range(kq) if r % p ** (k + 1))
    byv.setdefault(v_, []).append(b[r])
print('  by v_%d(N): %s' % (p, ', '.join('%d: %.6f' % (k, np.mean(byv[k])) for k in sorted(byv))))
