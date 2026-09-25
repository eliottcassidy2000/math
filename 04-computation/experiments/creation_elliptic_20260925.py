"""Intrinsic squareclass-controlled binary certificates on the fruit curve."""

from fractions import Fraction as F
from math import gcd, isqrt, lcm

A, B = 109, 224
G, T = (F(-4), F(28)), (F(56), F(728))


def check(test, label):
    if not test:
        raise RuntimeError(label)


def on_curve(p):
    return p is None or p[1]**2 == p[0]**3+A*p[0]**2+B*p[0]


def neg(p):
    return None if p is None else (p[0], -p[1])


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x,y = p
    u,v = q
    if x == u and y == -v:
        return None
    slope = (3*x*x+2*A*x+B)/(2*y) if p == q else (v-y)/(u-x)
    z = slope*slope-A-x-u
    return z, slope*(x-z)-y


def mul(n, p):
    if n < 0:
        return mul(-n, neg(p))
    result = None
    while n:
        if n % 2:
            result = add(result, p)
        p = add(p, p)
        n //= 2
    return result


def sqrt_rational(x):
    x = F(x)
    if x < 0:
        return None
    a,b = isqrt(x.numerator), isqrt(x.denominator)
    return F(a,b) if a*a == x.numerator and b*b == x.denominator else None


def alpha(p):
    if p is None:
        return 1
    if p[0] == 0:
        return 14
    matches = [c for c in (1,-1,14,-14) if sqrt_rational(p[0]/c) is not None]
    check(len(matches) == 1, "four-class subgroup domain")
    return matches[0]


def halves(p):
    """All rational halves, by exact square roots and final group checks."""
    check(on_curve(p), "halving curve domain")
    if p is None:
        return (None, (F(0),F(0)))
    x,y = p
    w = sqrt_rational(x)
    if w is None or w == 0:
        return ()
    result = set()
    for z in (2*x+2*y/w, 2*x-2*y/w):
        d = sqrt_rational(z*z-4*B)
        if d is None:
            continue
        for u in ((z+d)/2, (z-d)/2):
            v = sqrt_rational(u**3+A*u*u+B*u)
            if v is None:
                continue
            for signed_v in (v, -v):
                q = (u,signed_v)
                if add(q,q) == p:
                    result.add(q)
    return tuple(sorted(result))


def decode_step(p):
    """No free coefficient is supplied to this function."""
    check(on_curve(p), "decoder curve domain")
    c = alpha(p)
    r, t = int(c < 0), int(abs(c) == 14)
    remainder = add(add(p, neg(mul(r,G))), neg(mul(t,T)))
    choices = halves(remainder)
    check(len(choices) == 2, "exact two rational halves")
    selected = [q for q in choices if abs(alpha(q)) == 1]
    check(len(selected) == 1, "unique even-torsion section")
    q = selected[0]
    check(add(add(mul(2,q),mul(r,G)),mul(t,T)) == p, "certificate edge")
    return q,r,t


def fruit(p):
    if p is None:
        return (1,-1,0)
    x,y = p
    values = (56-x+y,56-x-y,-56-12*x)
    d = lcm(*(v.denominator for v in values))
    ints = [int(v*d) for v in values]
    common = gcd(gcd(abs(ints[0]),abs(ints[1])),abs(ints[2]))
    ints = tuple(v//common for v in ints)
    return tuple(-v for v in ints) if sum(ints) < 0 else ints


def rank(n):
    return n if n >= 0 else -n-1


def main():
    check(mul(6,T) is None and mul(3,T) == (0,0), "torsion6")
    bank = {add(mul(n,G),mul(k,T)) for n in (-1,0) for k in (0,2,4)}
    check(len(bank) == 6, "six terminal states")
    points = {(n,k):add(mul(n,G),mul(k,T)) for n in range(-24,25) for k in range(6)}
    edges = 0
    longest = 0
    for (n,k), p in points.items():
        check(on_curve(p), "curve")
        check(alpha(p) == (-1 if n % 2 else 1) * (14 if k % 2 else 1), "intrinsic parity class")
        q,r,t = decode_step(p)
        next_n = n//2
        next_k = next(j for j in (0,2,4) if (2*j+t-k) % 6 == 0)
        check(q == add(mul(next_n,G),mul(next_k,T)), "coefficient-independent decoder")
        check(rank(next_n) == rank(n)//2, "integer rank halves")
        check(p in bank or 2*rank(next_n)+next_k % 2 < 2*rank(n)+k % 2,
              "strict full integer rank outside bank")
        edges += 1
        current = p
        trace = []
        while current not in bank:
            check(len(trace) <= 8, "bounded certificate length")
            current,r,t = decode_step(current)
            trace.append((r,t))
        for r,t in reversed(trace):
            current = add(add(mul(2,current),mul(r,G)),mul(t,T))
        check(current == p, "full certificate reconstruction")
        longest = max(longest,len(trace))
    print('INTRINSIC DECODER all294 points n=-24..24,k0..5:',edges,'edges PASS; maxtrace',longest)

    # Independent squareclass homomorphism controls, including exceptional torsion cases.
    small = {(n,k):points[n,k] for n in range(-4,5) for k in range(6)}
    count = 0
    for p in small.values():
        for q in small.values():
            actual = alpha(add(p,q))
            check(sqrt_rational(F(alpha(p)*alpha(q),actual)) is not None, "chord squareclass homomorphism")
            count += 1
    print('KUMMER HOMOMORPHISM all54^2=',count,'pairs including O/U PASS')

    for n in range(1,13):
        for k in range(6):
            p = points[n,k]
            doubled = add(p,p)
            x,y = p
            check(doubled[0] == (x*x-B)**2/(4*y*y) and doubled[0] >= 0,
                  "all real doubles outside fruit-positive negative-x window")
    print('DUPLICATION POSITIVITY OBSTRUCTION all72 non-torsion test points PASS')

    a = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
    b = 36875131794129999827197811565225474825492979968971970996283137471637224634055579
    c = 4373612677928697257861252602371390152816537558161613618621437993378423467772036
    check(fruit(points[9,0]) == (a,b,c), "recovered exact9G triple")
    check(sum((F(a,b+c),F(b,a+c),F(c,a+b))) == 4, "fruit equality")
    current = points[9,0]
    labels = []
    while current not in bank:
        current,r,t = decode_step(current)
        labels.append((r,t))
    check(labels == [(1,0),(0,0),(0,0),(1,0)] and current is None, "9G binary certificate")
    check(halves(points[9,0]) == (), "positive9G has no rational raw half")
    print('9G certificate digits',labels,'terminalO; chain9G->4G->2G->G->O')
    print('9G descendants fruit-positive flags:',[all(v>0 for v in fruit(points[n,0])) for n in (9,4,2,1)])
    for n in range(1,25):
        for sigma in (-1,1):
            numerator = add(mul(3,points[n,0]),mul(sigma,G)) if n % 2 else points[n,0]
            candidates = [q for q in halves(numerator) if abs(alpha(q)) == 1]
            target = (3*n+sigma)//2 if n % 2 else n//2
            check(candidates == [mul(target,G)], "actual Collatz guarded elliptic lift")
    check(mul(14,G) != points[4,0], "creation step differs from Collatz at9")
    print('COLLATZ LIFT all48 signed-sheet steps n1..24 PASS; plus9G->14G differs from creation9G->4G')
    for base in (2,3,6,8):
        # Independent image/coset calculation in Z/(base*6) x Z/6.
        quotient = {(n % base,k % gcd(base,6)) for n in range(base*6) for k in range(6)}
        check(len(quotient) == base*gcd(base,6), "finite division quotient")
        print('QUOTIENT L/',base,'L has',len(quotient),'classes')
    check(decode_step(neg(G))[0] == neg(G), "terminal sign boundary fixed point")
    print('HOSTILE -G is decoderfixed; terminalbank/selectedtarget cannot be omitted')
    print('ALL CHECKS PASS; certificates verify elliptic subgroup membership, not Collatz convergence')


if __name__ == '__main__':
    main()
