#!/usr/bin/env python3
"""Exact controls for regular-root global-modulus obstructions; no optimizer census."""
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
import sympy as s

gates = 0


def check(value, label):
    global gates
    gates += 1
    if not value:
        raise RuntimeError(label)


class I:
    def __init__(self, a, b=None):
        self.a, self.b = F(a), F(a if b is None else b)
        if self.a > self.b:
            raise ValueError('unordered interval')
    def __add__(self, y):
        y = iv(y)
        return I(self.a+y.a, self.b+y.b)
    __radd__ = __add__
    def __neg__(self): return I(-self.b, -self.a)
    def __sub__(self, y): return self + -iv(y)
    def __rsub__(self, y): return iv(y) + -self
    def __mul__(self, y):
        y = iv(y)
        p = [self.a*y.a, self.a*y.b, self.b*y.a, self.b*y.b]
        return I(min(p), max(p))
    __rmul__ = __mul__
    def __truediv__(self, y):
        y = iv(y)
        if y.a <= 0 <= y.b: raise ValueError('zero divisor')
        return self * I(1/y.b, 1/y.a)
    def __rtruediv__(self, y): return iv(y) / self
    def display(self):
        q = 10**12
        return [str(F((self.a*q).__floor__(), q)),
                str(F((self.b*q).__ceil__(), q))]


def iv(x): return x if isinstance(x, I) else I(x)


def root(q):
    q = F(q)
    scale = 10**40
    k = isqrt(q.numerator*scale**2//q.denominator)
    lo, hi = F(k, scale), F(k+1, scale)
    check(lo**2 <= q <= hi**2, 'exact radical enclosure')
    return I(lo, hi)


u, v, x, n = s.symbols('u v x n')


def identity(expr, label):
    num = s.cancel(expr).as_numer_denom()[0]
    reduced = s.rem(s.rem(s.expand(num), u*u-2, u), v*v-3, v)
    check(s.expand(reduced) == 0, label)


def main():
    B = u-1
    K = 4*B*v/(3*(1+v))
    cs = (13-8*u)/3
    k4 = 2*B*(v-1)/9
    J = (5-8*x+3*x*x)/(3*(1-x*x))
    R = (J-cs)/(2-2*u*x)
    regular = 2*B/(3*(1+v)*(1+x))
    identity(R-4*B/(3*(1+x)), 'regular-root cancellation identity')
    identity((R-K)/(2-2*v*x)-regular, 'regular distance quotient')
    identity(regular.subs(x, s.Rational(1,2))-k4, 'four-root barrier')
    identity(regular.subs(x, s.Rational(4,7))-s.Rational(21,22)*k4,
             'fractional multiplicity loses exactly one twenty-second')
    identity(s.diff(regular,x)+2*B/(3*(1+v)*(1+x)**2),
             'strict regular-family monotonicity derivative')
    identity((R-K)/(x*x-2*x/v+s.Rational(1,3))
             -4*B/(3*(1+1/v)*(1+x)*(1/v-x)), 'regular moment quotient')

    den = 2*n+1
    p3 = (4*n**3+1-8*n)/den**3
    p4 = (4*n**4+1+16*n)/den**4
    Jn = (7*n**3+32*n**2+62*n+34)/(3*(3*n**3+8*n**2+6*n-2))
    identity((4*n+1-2*n)/den-1, 'rational family first moment')
    identity((4*n*n+1+4*n)/den**2-1, 'rational family second moment')
    identity((5-8*p3+3*p4)/(3*(1-p4))-Jn, 'compressed rational family duplication')
    # General normalization lift: alpha*S+n*beta=1 and alpha^2+n*beta^2=1.
    S, alpha = s.symbols('S alpha')
    beta = (1-alpha*S)/n
    check(s.expand(n*(alpha**2+n*beta**2-1)
                   -((n+S*S)*alpha**2-2*S*alpha+1-n)) == 0,
          'normalization lift quadratic identity')

    U, V = root(2), root(3)
    KI = 4*(U-1)*V/(3*(1+V))
    CSI = (13-8*U)/3
    K4 = 2*(U-1)*(V-1)/9
    KM = 4*(U-1)/(1+V)
    check(F(67383,10**6)<K4.a<K4.b<F(67384,10**6), 'simple four-root bracket')
    check(K4.b<F(1,10), 'four-root barrier below one tenth')

    rows = []
    for k in (1, 10, 100, 1000):
        d = 2*k+1
        r = [F(k,d)]*4+[F(1,d)]+[F(-2,d)]*k
        check(len(r)==k+5 and all(r), 'actual exact nonzero length')
        check(sum(r)==1 and sum(t*t for t in r)==1, 'literal actual normalization')
        pos = sorted((t for t in r if t>0), reverse=True)
        check(pos[:3]==[F(k,d)]*3, 'actual ordered top-three labels')
        p3q, p4q = sum(t**3 for t in r), sum(t**4 for t in r)
        jq = (5-8*p3q+3*p4q)/(3*(1-p4q))
        check(jq==F(7*k**3+32*k*k+62*k+34,3*(3*k**3+8*k*k+6*k-2)),
              'literal and compressed duplication agree')
        dist2 = 2-U*sum(pos[:2])
        delta = 2-2*V*F(k,d)
        moment = p4q-2*p3q/V+F(1,3)
        check(p4q<1 and dist2.a>0 and delta.a>0 and moment.a>0,
              'every physical quotient has positive denominator')
        rr = (jq-CSI)/dist2
        qd, qm = (rr-KI)/delta, (rr-KI)/moment
        if k==1000:
            check(F(74921,10**6)<qd.a<qd.b<F(74922,10**6), 'named actual hostile bracket')
            check(qd.b<F(1,10), 'actual hostile to global coefficient one tenth')
        # Independent degree-four truncation suffices to recover [s^4]G^2.
        coeff=[F(1),F(0),F(0),F(0),F(0)]
        for t in r:
            for j in range(4,0,-1): coeff[j]+=t*coeff[j-1]
        D=-sum(coeff[j]*coeff[4-j] for j in range(5))
        E=(1-p4q)/2
        check(D/E==jq, 'ordinary polynomial square independently gives J')
        rows.append(dict(n=k,length=k+5,J=str(jq),distance_quotient=qd.display(),
                         moment_quotient=qm.display()))

    manifest=dict(status='PROVED regular-family barriers; global sharp lower bounds OPEN',
                  universe='symbolic identities plus four literal rational rows n=1,10,100,1000',
                  kappa4=K4.display(),moment_upper_barrier=KM.display(),
                  rational_rows=rows,
                  fractional_model='formal multiplicity49/16 at x4/7; quotient21*kappa4/22; NOT an actual list')
    data=json.dumps(manifest,sort_keys=True,separators=(',',':'))
    print('EXACT REGULAR-FAMILY GLOBAL-MODULUS OBSTRUCTIONS')
    print(data)
    print('EXPLICIT_GATES',gates)
    print('SEMANTIC_SHA256',hashlib.sha256(data.encode()).hexdigest())


if __name__ == '__main__': main()
