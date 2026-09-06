#!/usr/bin/env python3
"""Exact bounded controls for the degree-five signed duplication margin.

No repository imports. Arithmetic is rational or in Q(sqrt(6)).
"""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from hashlib import sha256


@dataclass(frozen=True)
class Q6:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Q6) else Q6(other)

    def __add__(self, other):
        other = self.coerce(other)
        return Q6(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q6(-self.a, -self.b)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        return Q6(self.a*other.a+6*self.b*other.b,
                  self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        den = other.a**2-6*other.b**2
        if den == 0:
            raise ZeroDivisionError
        return self*Q6(other.a/den, -other.b/den)

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def __pow__(self, exponent):
        if exponent < 0:
            return (1/self)**(-exponent)
        out = Q6(1)
        for _ in range(exponent):
            out *= self
        return out

    def sign(self):
        def sg(x):
            return (x > 0)-(x < 0)
        a, b = self.a, self.b
        if b == 0:
            return sg(a)
        if a == 0 or sg(a) == sg(b):
            return sg(b)
        return sg(a*a-6*b*b) if a > 0 else sg(6*b*b-a*a)

    def __str__(self):
        return str(self.a) if self.b == 0 else f'({self.a})+({self.b})sqrt(6)'


checks = 0
trace = []


def gate(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)


def mul(a, b):
    out = [a[0]*0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def elementary(roots):
    one = roots[0]*0+1
    out = [one]
    for r in roots:
        out = mul(out, [one, r])
    return out


def data(roots):
    coeff = elementary(roots)
    energy = sum(x*x*y*y for x, y in combinations(roots, 2))
    doubled = mul(coeff, coeff)[4]
    return coeff, energy, -doubled


sharp = Q6(F(27, 29), -F(8, 87))
gate((sharp-F(1, 3)).sign() > 0, 'improve inherited terminal margin')
gate((F(7, 9)-sharp).sign() > 0, 'quartic-to-quintic strict decrease')

# Complete two-value orbits on sum r=sum r^2=1, up to permutation.
two_level = [
    (1, Q6(1), Q6(0), None),
    (1, Q6(-F(3, 5)), Q6(F(2, 5)), Q6(F(7, 3))),
    (2, Q6(F(1, 5), F(1, 5)), Q6(F(1, 5), -F(2, 15)), sharp),
    (2, Q6(F(1, 5), -F(1, 5)), Q6(F(1, 5), F(2, 15)), Q6(F(27, 29), F(8, 87))),
]
for multiplicity, x, y, expected in two_level:
    roots = [x]*multiplicity+[y]*(5-multiplicity)
    coeff, energy, d = data(roots)
    gate(sum(roots) == Q6(1), 'two-level first moment')
    gate(sum(r*r for r in roots) == Q6(1), 'two-level second moment')
    gate(coeff[2] == Q6(0), 'two-level e2 cancellation')
    p3, p4 = sum(r**3 for r in roots), sum(r**4 for r in roots)
    gate(energy == (1-p4)/2, 'two-level energy/power identity')
    gate(d == (5-8*p3+3*p4)/6, 'two-level doubled/power identity')
    gate((d-sharp*energy).sign() >= 0, 'two-level polynomial certificate')
    if expected is None:
        gate(energy == Q6(0) and d == Q6(0), 'zero-energy boundary')
    else:
        gate(d/energy == expected, 'two-level exact margin')
    trace.append(('two-level', multiplicity, str(x), str(y), str(energy), str(d)))

# The exact negative tangent curvature for every three-value multiplicity
# with singleton middle value in dimension five.
for left, right in [(1, 3), (2, 2), (3, 1)]:
    for u in map(F, [1, 2, 3]):
        for v in map(F, [1, 2, 3]):
            for leading in [F(1, 2), F(1), F(3)]:
                x, y, z = F(0), u, u+v
                dx, dy, dz = right*v, -left*right*(u+v), left*u
                gate(left*dx+dy+right*dz == 0, 'group tangent first moment')
                gate(left*x*dx+y*dy+right*z*dz == 0, 'group tangent second moment')
                qx, qy, qz = 4*leading*u*(u+v), -4*leading*u*v, 4*leading*v*(u+v)
                hessian = left*qx*dx*dx+qy*dy*dy+right*qz*dz*dz
                factored = -4*leading*left*right*u*v*(u+v)*(left*(right-1)*u+right*(left-1)*v)
                gate(hessian == factored < 0, 'strict negative group curvature')


def rational_audit(roots, label):
    coeff, energy, d = data(roots)
    gate(coeff[2] == 0, label+': cancellation')
    if energy == 0:
        gate(d == 0 and sum(r != 0 for r in roots) <= 1, label+': zero-energy boundary')
        trace.append((label, tuple(map(str, roots)), 'zero-energy'))
        return None
    gate((Q6(d)-sharp*energy).sign() > 0, label+': strict rational margin')
    p1 = sum(roots)
    gate(p1 != 0 and sum(r*r for r in roots) == p1*p1, label+': normalization eligible')
    normalized = [r/p1 for r in roots]
    p3, p4 = sum(r**3 for r in normalized), sum(r**4 for r in normalized)
    gate(d/energy == (5-8*p3+3*p4)/(3*(1-p4)), label+': moment quotient')
    gate(energy == -2*coeff[1]*coeff[3]+2*coeff[4], label+': coefficient energy')
    for scalar, variable in [(F(-2), F(3, 2)), (F(5, 3), F(-1, 2))]:
        changed = [scalar*variable**i*a for i, a in enumerate(coeff)]
        changed_energy = scalar*scalar*variable**4*energy
        gate(-mul(changed, changed)[4]/changed_energy == d/energy, label+': gauge')
    for shift in [1, 3]:
        shifted = [F(0)]*shift+coeff
        gate(-mul(shifted, shifted)[4+2*shift] == d, label+': monomial index')
    trace.append((label, tuple(map(str, roots)), str(d/energy)))
    return d/energy


full_rows = set()
full_rejected = 0
for first in combinations_with_replacement(map(F, [-3, -2, -1, 1, 2, 3]), 4):
    coeff = elementary(first)
    if coeff[1] == 0 or coeff[2] == 0:
        full_rejected += 1
        continue
    full_rows.add(tuple(sorted(first+(-coeff[2]/coeff[1],))))
full_ratios = [rational_audit(r, 'full') for r in sorted(full_rows)]

boundary_rows = set()
boundary_seeds = 0
boundary_rejected = 0
for first in combinations_with_replacement(map(F, [-2, -1, 0, 1, 2]), 4):
    if 0 not in first:
        continue
    boundary_seeds += 1
    coeff = elementary(first)
    if coeff[1] == 0:
        if all(r == 0 for r in first):
            boundary_rows.add(first+(F(0),))
        else:
            boundary_rejected += 1
        continue
    boundary_rows.add(tuple(sorted(first+(-coeff[2]/coeff[1],))))
for roots in sorted(boundary_rows):
    rational_audit(roots, 'zero-root boundary')

hostile = (F(1), F(1), -F(1, 6), -F(1, 6), -F(13, 60))
hostile_ratio = rational_audit(hostile, 'inherited degree-five hostile')
gate(hostile_ratio == F(18501, 26101), 'inherited hostile ratio')
gate(hostile_ratio < F(7, 9), 'inherited quartic bound refutation')

# Equality shape is tested by literal factors, independently of normalized
# two-level roots above. Its ratio of root magnitudes is 3+sqrt(6).
t = Q6(1, -F(1, 3))
eq_roots = [Q6(1), Q6(1), -t, -t, -t]
eq_coeff, eq_energy, eq_d = data(eq_roots)
gate(eq_coeff[2] == Q6(0), 'sharp unnormalized carrier cancellation')
gate(eq_d/eq_energy == sharp, 'sharp unnormalized carrier attainment')
gate(1/t == Q6(3, 1), 'sharp multiplicity/sign root geometry')

print('PASS sharp degree-five signed duplication margin')
print('Sharp constant: (81-8sqrt(6))/87')
print('Equality roots: scalar*(1,1,-t,-t,-t), t=1-sqrt(6)/3')
print('Two-level normalized orbits: 4; zero-energy boundary: 1; negative tangent controls: 81')
print('Full rational universe: 126 seeds;', full_rejected, 'rejections;', len(full_rows), 'unique rows')
print('Full rational ratio range:', min(full_ratios), max(full_ratios))
print('Zero-root boundary universe:', boundary_seeds, 'seeds;', boundary_rejected, 'infeasible seeds;', len(boundary_rows), 'unique rows')
print('Inherited hostile: 18501/26101; quartic bound 7/9 fails; sharp degree-five bound passes')
print('Exact gates:', checks)
print('Trace SHA256:', sha256(repr(trace).encode()).hexdigest())
