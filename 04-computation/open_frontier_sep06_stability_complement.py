#!/usr/bin/env python3
"""Exact complementary-region certificates for the sharp global constant.

Standard library only. Formal polynomial identities over Q(sqrt2,sqrt3),
radical sign comparisons, and a declared 49-prefix actual-row control bank.
The proof over the continuous domain is in the companion note.
"""
from fractions import Fraction as R
from itertools import product
import hashlib
import json

gates = 0


def require(test, name):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(name)


def rational_sign(x):
    return (x > 0) - (x < 0)


def quadratic_sign(a, b):
    """Sign of a+b sqrt2 without floating arithmetic."""
    sa, sb = rational_sign(a), rational_sign(b)
    if not sa:
        return sb
    if not sb or sa == sb:
        return sa
    return sa*rational_sign(a*a-2*b*b)


class K:
    """Coefficient order 1,sqrt2,sqrt3,sqrt6."""
    def __init__(self, x=0):
        if isinstance(x, K):
            self.c = x.c
        elif isinstance(x, (tuple, list)):
            self.c = tuple(map(R, x))
        else:
            self.c = (R(x), R(0), R(0), R(0))

    def __add__(self, other):
        other = K(other)
        return K(tuple(a+b for a,b in zip(self.c, other.c)))

    __radd__ = __add__

    def __neg__(self):
        return K(tuple(-a for a in self.c))

    def __sub__(self, other):
        return self+-K(other)

    def __rsub__(self, other):
        return K(other)+-self

    def __mul__(self, other):
        other = K(other)
        out = [R(0)]*4
        for i,a in enumerate(self.c):
            for j,b in enumerate(other.c):
                common = i & j
                factor = (2 if common & 1 else 1)*(3 if common & 2 else 1)
                out[i ^ j] += factor*a*b
        return K(out)

    __rmul__ = __mul__

    def inverse(self):
        a,b,c,d = self.c
        conjugate = K((a,b,-c,-d))
        norm = self*conjugate
        if norm.c[2:] != (0,0):
            raise RuntimeError('field norm implementation')
        aa,bb = norm.c[:2]
        denominator = aa*aa-2*bb*bb
        if denominator == 0:
            raise ZeroDivisionError
        return conjugate*K((aa/denominator,-bb/denominator,0,0))

    def __truediv__(self, other):
        return self*K(other).inverse()

    def __rtruediv__(self, other):
        return K(other)*self.inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse()**(-n)
        out = K(1)
        for _ in range(n):
            out *= self
        return out

    def __eq__(self, other):
        return self.c == K(other).c

    def sign(self):
        a,b,c,d = self.c
        sa,sb = quadratic_sign(a,b),quadratic_sign(c,d)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        difference = K((a,b,0,0))**2-3*K((c,d,0,0))**2
        return sa*quadratic_sign(*difference.c[:2])

    def serial(self):
        return list(map(str,self.c))


def clean(p):
    return {ij: K(c) for ij,c in p.items() if K(c) != 0}


def const(c):
    return clean({(0,0): K(c)})


def add(*ps):
    out = {}
    for p in ps:
        for ij,c in p.items():
            out[ij] = out.get(ij,K(0))+c
    return clean(out)


def scale(p,c):
    return clean({ij: x*c for ij,x in p.items()})


def mul(*ps):
    out = const(1)
    for p in ps:
        next_p = {}
        for (i,j),c in out.items():
            for (k,l),d in p.items():
                ij = i+k,j+l
                next_p[ij] = next_p.get(ij,K(0))+c*d
        out = clean(next_p)
    return out


def power(p,n):
    return mul(*([p]*n))


def derivative(p, coordinate):
    out = {}
    for ij,c in p.items():
        degree = ij[coordinate]
        if degree:
            key = list(ij)
            key[coordinate] -= 1
            out[tuple(key)] = degree*c
    return clean(out)


def fixed_b(p,b):
    return add(*(const(c*b**j) if i == 0 else {(i,0): c*b**j}
                 for (i,j),c in p.items()))


def diagonal(p):
    return add(*({(i+j,0): c} for (i,j),c in p.items()))


u,v = K((0,1,0,0)),K((0,0,1,0))
r,h = 1/v,1/u
A,B = 2-u,u-1
sharp = 4*v/(3*(1+u)*(1+v))
gamma = 3*sharp/(4*u)
alpha = 3*sharp/8
C0 = A+u*gamma
require(u*u == 2 and v*v == 3 and (u*v)**2 == 6, 'field basis')
require(gamma == (2-u)*(3-v)/4, 'gamma exact simplification')
require(C0 == (1+u+v-u*v)/2, 'C0 exact simplification')
require(alpha*u == gamma, 'objective coefficient identity')
require((C0-4*gamma*r).sign() > 0, 'secant maximum increases through r')
require((2-v).sign() > 0, 'strict secant derivative bound')
require((2*C0-12*gamma*r).sign() > 0, 'U derivative positive on entire region')
require((2*C0*r-1).sign() < 0, 'U(r) strictly negative')
require((R(8,3)-2*(1+u)/v).sign() < 0, 'boundary cubic derivative upper bound')
require((2*u*(2*v-3)/9).sign() > 0, 'boundary cubic endpoint lower bound')
require((4*(13*u-18)/3-R(1,2)).sign() > 0, 'two-atom coarse constant exceeds one half')
require((K(R(1,2))-sharp).sign() > 0, 'sharp constant below one half')

a,b = {(1,0):K(1)},{(0,1):K(1)}
one = const(1)
a2,b2 = power(a,2),power(b,2)
remaining = add(one,scale(a2,-1))
t = add(a,b)
d2 = add(const(2),scale(t,-u))
C = add(const(C0),scale(t,-gamma))
F = add(const(B),scale(d2,-alpha),scale(power(a,3),-1),
        scale(mul(remaining,b),-1),
        mul(C,add(power(a,4),mul(remaining,b2))))
T = add(scale(add(one,a2,scale(mul(a,b),-2),scale(b2,-3)),gamma),
        scale(b,2*C0),const(-1))
U = add(scale(add(const(2),scale(b2,-6)),gamma),scale(b,2*C0),const(-1))
require(derivative(F,1) == mul(remaining,T), 'formal envelope derivative')
require(add(U,scale(T,-1)) == scale(add(one,scale(a2,-1),scale(b2,-1),
        scale(mul(b,add(a,scale(b,-1))),2)),gamma), 'formal domain remainder')
diag_factor = scale(mul(add(one,scale(a,-1)),add(a,const(-r)),add(a,const(-h))),2*gamma)
require(diagonal(F) == diag_factor, 'formal diagonal factorization')
H = add(power(a,3),scale(a2,-(1+u)),scale(a,R(2,3)),const(R(2,3)))
boundary_factor = scale(mul(add(one,scale(a,-1)),add(a,const(-r)),H),gamma)
require(fixed_b(F,r) == boundary_factor, 'formal b=r boundary factorization')
endpoint_H = sum(c*(u/v)**i for (i,j),c in H.items())
require(endpoint_H == 2*u*(2*v-3)/9, 'boundary cubic endpoint identity')

# Independent univariate two-atom moment elimination (variable a means t).
t2minus1 = add(a2,const(-1))
E2 = scale(power(t2minus1,2),R(1,4))
p3 = scale(add(scale(a,3),scale(power(a,3),-1)),R(1,2))
p4 = add(one,scale(power(t2minus1,2),R(-1,2)))
g2 = add(const(B),scale(p3,-1),scale(p4,A))
twod = add(const(2),scale(a,-u))
require(mul(g2,power(add(a,one),2)) ==
        scale(mul(add(scale(a,A),one),E2,twod),u), 'exact two-atom ratio identity')

# Signed actual controls: a rational two-root prefix, then solve e2=0.
grid = list(map(R,[-2,-1]))+[R(-1,2),R(0),R(1,2),R(1),R(2)]
trace,exceptions = [],[]
for x,y in product(grid,repeat=2):
    if x+y == 0:
        exceptions.append([str(x),str(y),'zero prefix sum'])
        continue
    raw = (x,y,-x*y/(x+y))
    S = sum(raw)
    if not S:
        exceptions.append([str(x),str(y),'zero normalization'])
        continue
    row = tuple(z/S for z in raw)
    require(sum(row) == sum(z*z for z in row) == 1, 'actual first two moments')
    energy = (1-sum(z**4 for z in row))/2
    if not energy:
        exceptions.append([str(x),str(y),'zero energy'])
        continue
    positives = sorted((z for z in row if z>0),reverse=True)+[R(0),R(0)]
    aa,bb = positives[:2]
    distance = 2-u*(aa+bb)
    # Literal product/square coefficient, independently of the moment formula.
    coeffs = [R(1)]
    for z in row:
        out = [R(0)]*(len(coeffs)+1)
        for i,c in enumerate(coeffs):
            out[i] += c
            out[i+1] += c*z
        coeffs = out
    doubled4 = sum(coeffs[i]*coeffs[4-i] for i in range(len(coeffs))
                   if 0<=4-i<len(coeffs))
    actual_gap = -doubled4/energy-(13-8*u)/3-sharp*distance
    require(actual_gap.sign()>0, 'actual sharp global strict inequality')
    region = 'complement' if (K(bb)-r).sign()<=0 else 'packing'
    trace.append([str(x),str(y),region,actual_gap.serial()])

print('PASS complementary small-second-root region; combined with packing theorem gives exact global K3')
print('K3=4sqrt3/[3(1+sqrt2)(1+sqrt3)]; strict on every finite eligible normalized list')
print('Formal gates: signed secant bounds; derivative/domain remainder; diagonal; fixed b=1/sqrt3; two-atom ratio')
print('ACTUAL_PREFIX_UNIVERSE',len(grid)**2,'ELIGIBLE',len(trace),'EXCEPTIONS',len(exceptions))
print('EXCEPTIONS',json.dumps(exceptions,separators=(',',':')))
print('REGIONS',json.dumps({s:sum(row[2]==s for row in trace) for s in ('complement','packing')},sort_keys=True))
print('EXACT_GATES',gates)
print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(trace,separators=(',',':')).encode()).hexdigest())
