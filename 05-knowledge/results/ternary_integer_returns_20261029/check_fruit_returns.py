#!/usr/bin/env python3
"""Exact real-circle cover, over Q(sqrt(3)), on both sides of a 3-isogeny.

This checks the entire circle by its ordered algebraic endpoints, not by
sampling rational points or estimating elliptic logarithms. No floats.
"""
from collections import Counter
from fractions import Fraction as F
from functools import total_ordering
from hashlib import sha256
from itertools import groupby
from pathlib import Path
import json

CHECKS = 0


def require(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(msg)


@total_ordering
class Quad:
    """The real embedding of Q(sqrt(3)); comparisons are rational signs."""
    def __init__(self, a=0, b=0):
        if isinstance(a, Quad):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, v):
        v = Quad(v)
        return Quad(self.a+v.a, self.b+v.b)

    __radd__ = __add__

    def __neg__(self):
        return Quad(-self.a, -self.b)

    def __sub__(self, v):
        return self + -Quad(v)

    def __rsub__(self, v):
        return Quad(v) + -self

    def __mul__(self, v):
        v = Quad(v)
        return Quad(self.a*v.a+3*self.b*v.b, self.a*v.b+self.b*v.a)

    __rmul__ = __mul__

    def __truediv__(self, v):
        v = Quad(v)
        norm = v.a*v.a-3*v.b*v.b
        return self*Quad(v.a/norm, -v.b/norm)

    def __rtruediv__(self, v):
        return Quad(v)/self

    def __pow__(self, n):
        if n < 0:
            return 1/(self**-n)
        result, x = Quad(1), self
        while n:
            if n & 1:
                result = result*x
            x, n = x*x, n//2
        return result

    def sign(self):
        a, b = self.a, self.b
        if not a:
            return (b > 0)-(b < 0)
        if not b or a*b > 0:
            return (a > 0)-(a < 0)
        norm = a*a-3*b*b
        return ((a > 0)-(a < 0))*((norm > 0)-(norm < 0))

    def __eq__(self, v):
        v = Quad(v)
        return (self.a, self.b) == (v.a, v.b)

    def __lt__(self, v):
        return (self-v).sign() < 0

    def data(self):
        return [str(self.a), str(self.b)]


def point(x, y):
    return Quad(x), Quad(y)


def neg(p):
    return (p[0], -p[1]) if p is not None else None


def add(p, q, a2=109, a4=224):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    z, w = q
    if x == z:
        if y+w == 0:
            return None
        slope = (3*x*x+2*a2*x+a4)/(2*y)
    else:
        slope = (w-y)/(z-x)
    X = slope*slope-a2-x-z
    return X, slope*(x-X)-y


def mul(n, p, a2=109, a4=224):
    if n < 0:
        return mul(-n, neg(p), a2, a4)
    result = None
    while n:
        if n & 1:
            result = add(result, p, a2, a4)
        p, n = add(p, p, a2, a4), n//2
    return result


def on_curve(p, a2, a4):
    return p is None or p[1]**2 == p[0]**3+a2*p[0]**2+a4*p[0]


def circle_key(p):
    # All these points lie on the bounded component and avoid y=0.
    x, y = p
    require(y != 0, 'endpoint/anchor at a turning point')
    return (0, x) if y > 0 else (1, -x)


def positive(p):
    if p is None:
        return False
    x, y = p
    return x < F(-14, 3) and x*x+112*x+784 > 0


def chamber(p):
    return p is not None and p[0] < F(-686, 3)


def rational_image_guard(p):
    require(p is not None and p[0].b == 0 and p[1].b == 0, 'rational affine input')
    q = -p[0].a/196-1
    return q.denominator % 3 == 0 or (q.numerator+q.denominator) % 3 == 0


def steer_quotient(p, H):
    require(on_curve(p, 421, 43904) and p[0] < 0, 'bounded rational quotient input')
    candidates = [(e, add(p, mul(e, H, 421, 43904), 421, 43904)) for e in (0, 2, 4)]
    accepted = [(e, q) for e, q in candidates if rational_image_guard(q)]
    require(len(accepted) == 1, 'unique rational lifting phase')
    e, q = accepted[0]
    step = mul(6, H, 421, 43904)
    for j in range(13):
        if chamber(q):
            break
        q = add(q, step, 421, 43904)
    require(chamber(q) and rational_image_guard(q), 'positive rational image')
    require(add(q, neg(mul(e+6*j, H, 421, 43904)), 421, 43904) == p,
            'recover exact input without a subgroup index')
    return e, j, q


def fruit(a, b, c):
    z = 6*(a+b)-c
    return -28*(a+b+2*c)/z, 364*(a-b)/z


def isogeny(p):
    if p is None or p[0] == 4:
        return None
    x, y = p
    return (x*(x-56)**2/(x-4)**2,
            (x-56)*(x*x+44*x+224)*y/(x-4)**3)


def exact_cover(boundary, anchor, step, member, a2, a4):
    """Track exact membership on every cell AND every boundary point."""
    boundary = sorted(boundary, key=circle_key)
    events, active, translated, at = [], set(), boundary[:], anchor
    powers = [None]
    for j in range(14):
        if member(at):
            active.add(j)
        for i, p in enumerate(translated):
            require(on_curve(p, a2, a4) and p[0] < 0, 'bounded endpoint')
            # The leftmost real point is inside the target. Orientation by
            # dx/(2y) starts there: first crossing exits, next enters, etc.
            events.append((p, -1 if i % 2 == 0 else 1, j))
        at = add(at, step, a2, a4)
        translated = [add(p, neg(step), a2, a4) for p in translated]
        powers.append(add(powers[-1], step, a2, a4))

    events.sort(key=lambda e: circle_key(e[0]))
    split = sum(circle_key(p) < circle_key(anchor) for p, d, j in events)
    events = events[split:]+events[:split]
    initial = active.copy()
    minimum = {m: 100 for m in range(1, 14)}
    return_cells = Counter()
    records = []

    def inspect(state, kind):
        for m in minimum:
            minimum[m] = min(minimum[m], len(state & set(range(1, m+1))))
        if 0 in state:
            future = state-{0}
            require(bool(future), 'positive source with no return by 13')
            return_cells[str(min(future))+':'+kind] += 1

    inspect(active, 'cell')
    for _, group in groupby(events, key=lambda e: circle_key(e[0])):
        vals = list(group)
        exits = {j for p, d, j in vals if d < 0}
        enters = {j for p, d, j in vals if d > 0}
        require(not exits & enters, 'same interval enters and exits together')
        at_boundary = active-exits
        p = vals[0][0]
        # A separate exact evaluation verifies every boundary membership bit,
        # including simultaneous entry/exit events. Open intervals stay open.
        literal = {j for j in range(14)
                   if member(add(p, powers[j], a2, a4))}
        require(at_boundary == literal, 'event sweep vs direct membership')
        inspect(at_boundary, 'boundary')
        active = at_boundary | enters
        inspect(active, 'cell')
        records.append(dict(point=[q.data() for q in p],
                            exits=sorted(exits), enters=sorted(enters)))
    require(active == initial, 'closed circle sweep')
    require(minimum[13] >= 1 and minimum[12] == 0, 'sharp cover size')
    require({int(k.split(':')[0]) for k in return_cells} == {4, 9, 13},
            'first-return alphabet')
    return dict(minimum_coverage=minimum, return_cells=dict(return_cells),
                distinct_boundaries=len(records),
                boundary_digest=sha256(json.dumps(records, sort_keys=True).encode()).hexdigest())


def main():
    G, T, Hprime = point(-4, 28), point(56, 728), point(-196, 196)
    step = add(G, G)
    require(step == point(F(676, 49), F(55796, 343)), '2G')
    boundary = []
    for sign in (-1, 1):
        a, b, c = Quad(0), Quad(2, sign), Quad(1)
        boundary.extend((fruit(a, b, c), fruit(c, a, b), fruit(b, c, a)))
    for p in boundary:
        require(on_curve(p, 109, 224) and not positive(p), 'original boundary')
    source = exact_cover(boundary, G, step, positive, 109, 224)

    phiG, phistep = isogeny(G), isogeny(step)
    require(phiG == mul(3, Hprime, 421, 43904), 'phi(G)=3H')
    require(phistep == mul(6, Hprime, 421, 43904), 'step descends to 6H')
    quotient_boundary = [(Quad(F(-686, 3)), Quad(0, sign*F(686, 9)))
                         for sign in (-1, 1)]
    for p in boundary:
        require(isogeny(p) in quotient_boundary, 'six endpoints become two')
    quotient = exact_cover(quotient_boundary, phiG, phistep, chamber, 421, 43904)

    witnesses = []
    for n, torsion, expected in ((9, 0, 4), (17, 0, 13), (21, 1, 9)):
        p = add(mul(n, G), mul(torsion, T))
        require(positive(p), 'positive witness')
        q = p
        for j in range(1, 14):
            q = add(q, step)
            if positive(q):
                break
        require(j == expected, 'first-return witness')
        require(isogeny(q) == add(isogeny(p), mul(j, phistep, 421, 43904), 421, 43904),
                'source/clock transport')
        witnesses.append(dict(source=[n, torsion], target=[n+2*j, torsion], clock=j))

    steering = []
    for n in (1, 3, 5, 35, 53):
        e = 2*((-2*n) % 3)
        require(e in (0, 2, 4) and (n+e) % 3 == 0, 'component/lifting correction')
        m = (n+e)//3
        p = mul(m, G)
        for j in range(13):
            if positive(p):
                break
            p = add(p, step)
        require(positive(p) and j <= 12, 'entry with zero-time option')
        target = n+e+6*j
        require(target-n <= 76, 'bounded source coefficient displacement')
        require(isogeny(p) == mul(target, Hprime, 421, 43904), 'steered lift')
        steering.append(dict(source=n, correction=e, clock=j, target=target,
                             lift_coefficient=target//3))
    require(steering[-1] == dict(source=53, correction=4, clock=12,
                                target=129, lift_coefficient=43), 'sharp 76 witness')

    full_steering = []
    for n in (-5, 1, 3, 5, 35, 53):
        for torsion in (0, 1):
            p = add(mul(n, Hprime, 421, 43904),
                    point(0, 0) if torsion else None, 421, 43904)
            # The controller receives p, not n or its torsion label.
            e, j, q = steer_quotient(p, Hprime)
            require((n+e) % 3 == 0 and e+6*j <= 76, 'full-group coefficient audit')
            lift = add(mul((n+e+6*j)//3, G), mul(torsion, T))
            require(positive(lift) and isogeny(lift) == q, 'authenticated named lift')
            full_steering.append(dict(source_index=n, torsion=torsion,
                                      correction=e, clock=j))

    result = dict(status='PASS; all real bounded-component phases covered',
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  explicit_checks=CHECKS, source_circle=source,
                  quotient_circle=quotient, rational_witnesses=witnesses,
                  steering_witnesses=steering,
                  index_free_steering_controls=full_steering,
                  scope='sharp uniform clock 13; first positive return clocks 4,9,13; no Collatz implication')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
