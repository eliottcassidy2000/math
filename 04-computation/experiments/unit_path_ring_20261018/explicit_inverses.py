#!/usr/bin/env python3
"""Explicit unit walks to 1/19, 1/3, 1/4 (and 1/37, 1/5) in K = Q(zeta15, sqrt-11) (mac-mini-2026-10-18).
Mechanism (building on clock_duality's B2, cubic_adelic_20261016e.md): for gamma in K^*, delta = conj(gamma)/gamma
is a unit vector, and delta + 1/delta = t^2/n - 2 with t = gamma + conj(gamma), n = gamma conj(gamma).
Unit vectors generate an additive group Gamma that is a RING (products of unit vectors are unit vectors).
So t^2/n is in Gamma; with t, n rational and coprime, Bezout gives 1/n in Gamma.
Checks everything in exact arithmetic with the repo's Field class (escaped_ports_20261017d.make('sqrt3'))."""
import os, sys
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from escaped_ports_20261017d import make
f, z, j, alpha, ports, units, eta, conj, escape, reduce_ = make('sqrt3')
bar = lambda x: f.apply(x, conj)
isunit = lambda x: f.mul(x, bar(x)) == f.one
h3 = f.sub(f.scale(eta, 2), f.one)          # eta = 1 + zeta5? check below; we build sqrt(-3) from zeta15 directly
z15 = z
zeta3 = f.power(z15, 5); s3 = f.add(f.scale(zeta3, 2), f.one)        # sqrt(-3) = 2 zeta3 + 1
assert f.mul(s3, s3) == f.scalar(-3)
zeta5 = f.power(z15, 3); s5 = f.add(f.scale(f.add(zeta5, f.power(zeta5, 4)), 2), f.one)   # sqrt5 = 2(z5 + z5^4) + 1
assert f.mul(s5, s5) == f.scalar(5)
s15 = f.mul(s3, s5); assert f.mul(s15, s15) == f.scalar(-15)
assert f.mul(j, j) == f.scalar(-11)
def walk_inverse(gamma, name):
    gb = bar(gamma); t = f.add(gamma, gb); n = f.mul(gamma, gb)
    delta = f.mul(gb, f.inv(gamma)); assert isunit(delta)
    s = f.add(delta, bar(delta))
    print('%-28s delta unit: %s ; t = %s ; n = %s ; delta + 1/delta = t^2/n - 2: %s' % (name, isunit(delta), t[0], n[0], s == f.sub(f.mul(f.mul(t, t), f.inv(n)), f.scalar(2))))
    return delta, t, n
# 1/19 from the Eisenstein integer (7 + 3 sqrt-3)/2 of norm 19
g19 = f.scale(f.add(f.scalar(7), f.scale(s3, 3)), Q(1, 2)); d19, t19, n19 = walk_inverse(g19, '(7 + 3 sqrt-3)/2')
w19 = f.sub(f.scale(f.add(d19, bar(d19)), 7), f.scalar(4))
print('   7 delta + 7 delta^-1 - 4 == 1/19:', w19 == f.scalar(Q(1, 19)), '(18 unit steps)')
# 1/3 from (1 + sqrt-11)/2 of norm 3
g3 = f.scale(f.add(f.one, j), Q(1, 2)); d3, _, _ = walk_inverse(g3, '(1 + sqrt-11)/2')
print('   delta + delta^-1 + 2 == 1/3:', f.add(f.add(d3, bar(d3)), f.scalar(2)) == f.scalar(Q(1, 3)), '(4 unit steps)')
# 1/4 from (1 + sqrt-15)/2 of norm 4
g2 = f.scale(f.add(f.one, s15), Q(1, 2)); d2, _, _ = walk_inverse(g2, '(1 + sqrt-15)/2')
print('   delta + delta^-1 + 2 == 1/4:', f.add(f.add(d2, bar(d2)), f.scalar(2)) == f.scalar(Q(1, 4)), '(4 unit steps)')
# 1/37 from (11 + 3 sqrt-3)/2 of norm 37: 121/37 - 2 = delta + 1/delta ; 1/37 = -11 (121/37) + 36
g37 = f.scale(f.add(f.scalar(11), f.scale(s3, 3)), Q(1, 2)); d37, _, _ = walk_inverse(g37, '(11 + 3 sqrt-3)/2')
w37 = f.add(f.scale(f.add(f.add(d37, bar(d37)), f.scalar(2)), -11), f.scalar(36))
print('   -11 (delta + delta^-1 + 2) + 36 == 1/37:', w37 == f.scalar(Q(1, 37)))
# 1/5 from (3 + sqrt-11)/2 of norm 5: 9/5 - 2 ; 1/5 = 4 (9/5) - 7
g5 = f.scale(f.add(f.scalar(3), j), Q(1, 2)); d5, _, _ = walk_inverse(g5, '(3 + sqrt-11)/2')
w5 = f.sub(f.scale(f.add(f.add(d5, bar(d5)), f.scalar(2)), 4), f.scalar(7))
print('   4 (delta + delta^-1 + 2) - 7 == 1/5:', w5 == f.scalar(Q(1, 5)))
# sqrt-11 from rho = (5 + sqrt-11)/6 (a spindle point): j = 6 rho - 5 (11 unit steps)
rho = f.scale(f.add(f.scalar(5), j), Q(1, 6)); print('rho = (5 + sqrt-11)/6 unit:', isunit(rho), '; 6 rho - 5 == sqrt-11:', f.sub(f.scale(rho, 6), f.scalar(5)) == j)
