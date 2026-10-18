#!/usr/bin/env python3
"""Referee check B6 (second S41 audit): the two contact-family identities, symbolically.

Only the relations used: conj(u) = 1/u, Z = d + u, conj(Z) = dbar + 1/u, s^2 = 1 - 4/(Z conj Z), conj(s) = -s,
v = Z(-1 + s)/2, conj(v) = conj(Z)(-1 - s)/2, e^2 = -Z/conj(Z).
  (a) u + v + dbar u v = u v conj(d + u + v), hence |u + v + dbar uv| = |d + u + v| = 1;
  (b) |u + v + dbar uv|^2 = 1 directly;
  (c) d + u + dbar e^2 = Z conj(u)/conj(Z), a unit since |u| = 1 and |Z/conj Z| = 1;
  (d) e^2 = (Z/R)^2 s^2 h^2 = Z^2 D/R equals -Z/conj(Z) (with D = 1 - 4/N, R = 4 - N, N = Z conj Z).
Usage: python b6_identities.py
"""
import sympy as sp

u, d, db, s = sp.symbols('u d dbar s')
Z = d + u; Zb = db + 1 / u
N = Z * Zb
v = Z * (-1 + s) / 2; vb = Zb * (-1 - s) / 2
red = lambda ex: sp.simplify(sp.expand(sp.numer(sp.together(ex))).subs(s ** 2, 1 - 4 / N))
print('|v|^2 - 1 =', sp.simplify(red(v * vb - 1)))
print('|d + u + v|^2 - 1 =', sp.simplify(red((Z + v) * (Zb + vb) - 1)))
lhs = u + v + db * u * v
print('(a) u + v + dbar uv - uv conj(d+u+v), using v conj(v) = 1:',
      sp.simplify(red(sp.expand(lhs - u * v * (db + 1 / u + vb)) - u * (v * vb - 1))))
lhsb = 1 / u + vb + d * (1 / u) * vb
print('(b) |u + v + dbar uv|^2 - 1 =', sp.simplify(red(lhs * lhsb - 1)))
e2 = -Z / Zb
print('(c) d + u + dbar e^2 - Z conj(u)/conj(Z) =', sp.simplify(d + u + db * e2 - Z * (1 / u) / Zb))
Dv = 1 - 4 / N; Rv = 4 - N
print('(d) Z^2 D/R + Z/conj(Z) =', sp.simplify(Z ** 2 * Dv / Rv + Z / Zb))
