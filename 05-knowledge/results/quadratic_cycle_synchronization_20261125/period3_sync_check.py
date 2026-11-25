"""Exact replays for the integer classification, synchronized fibering and period 3."""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd, isqrt
from itertools import permutations


def fixed_graph(a, d, bound=2):
    """All affine rational preperiodic points when the absolute escape bound is bound."""
    assert gcd(a, d) == 1
    nodes = {x for x in range(-bound*d, bound*d+1) if gcd(x, d) == 1}
    edges = {x: (x*x-a)//d for x in nodes if (x*x-a) % d == 0}
    while True:
        nxt = {x: y for x, y in edges.items() if y in edges}
        if nxt == edges:
            return nxt
        edges = nxt


def cycle_data(a, b):
    assert a > 0 and b > 0 and gcd(a, b) == 1
    t = Q(a, b)
    u, v, w = -1-1/t, t, -1/(1+t)
    x, y, z = (u+w-v)/2, (u+v-w)/2, (v+w-u)/2
    param = x*x-z
    F = a**3+a*a*b-2*a*b*b-b**3
    d = 2*a*b*(a+b)
    assert x < y < z and (z*z-param, y*y-param) == (y, x)
    assert param == Q(7,4)+Q(F*F,d*d)
    assert param.denominator == d*d and all(p.denominator == d for p in (x,y,z))
    return param, d, (int(x*d), int(z*d), int(y*d)), F


@lru_cache(None)
def prime_count(k):
    n = 0
    p = 2
    while p*p <= k:
        if k % p == 0:
            n += 1
            while k % p == 0:
                k //= p
        p += 1
    return n + (k > 1)


def check_integer():
    def classifier(z):
        if z < 0:
            return set()
        n = (isqrt(4*z+1)-1)//2
        if z == n*(n+1):
            return {-n,n,-(n+1),n+1} | ({0} if n == 1 else set())
        if z == n*(n+1)+1:
            return {-n,n,-(n+1),n+1}
        return set()
    for z in range(-10, 10001):
        bound = 2+isqrt(abs(z))
        assert fixed_graph(z, 1, bound).keys() == classifier(z)
    for z in range(8):
        print("integer", z, fixed_graph(z,1,2+isqrt(z)))


def check_param():
    cases = [
        (1,1,Q(29,16),4,(-7,5,-1),-1,8),
        (5,4,Q(226801,129600),360,(-629,469,-19),1,6),
        (3,2,Q(6469,3600),60,(-107,83,7),13,6)]
    for a,b,param,d,cyc,F,count in cases:
        assert cycle_data(a,b) == (param,d,cyc,F)
        graph = fixed_graph(param.numerator,d)
        assert tuple(graph[x] for x in cyc) == cyc[1:] + cyc[:1]
        assert len(graph) == count
        print("period 3", param, "at denominator", d, ":", cyc, "preperiodic points:", sorted(graph))
    unit_examples = []
    for b in range(1,201):
        for a in range(1,2*b+1):
            if gcd(a,b) != 1:
                continue
            F = a**3+a*a*b-2*a*b*b-b**3
            if abs(F) <= 7:
                unit_examples.append((a,b,F))
            d = 2*a*b*(a+b)
            A = 7*(d//2)**2+F*F
            assert (F % 7 == 0) == ((a-2*b) % 7 == 0)
            if F % 7 == 0:
                assert F % 49 != 0
            if A <= 2*d*d and prime_count(d) <= 2:
                assert (a,b) == (1,1)
    assert unit_examples == [(1,1,-1),(2,1,7),(5,4,1)]
    print("bounded unit and 7-norm replay:", unit_examples)


def check_synchronization():
    a = (-5,-1,-3,1)
    d = 4
    tags = tuple(x*x-d*y for x,y in zip(a,a[1:]+a[:1]))
    assert tags == (29,13,5,21)
    for A in tags:
        assert -Q(1,4) <= Q(A,d*d) <= 2
        m = fixed_graph(A,d)
        for xs in permutations(m,4):
            assert tuple(m[x] for x in xs) != xs[1:]+xs[:1]
    print("untagged spurious 4-cycle",a,"requires four incompatible tags",tags)


def local_reduction_table(p):
    """Possible lift periods at good p>=5, by the reduced cycle and multiplier."""
    assert p >= 5
    possible = set()
    for z in range(p):
        for x in range(p):
            y, C, mult = x, set(), 1
            while y not in C:
                C.add(y)
                mult = 2*y*mult % p
                y = (y*y-z) % p
            if x != y:
                continue
            m = len(C)
            possible.add(m)
            if mult:
                r = next(k for k in range(1,p) if pow(mult,k,p) == 1)
                possible.add(m*r)
    return possible


def check_local():
    primes = [5,7,11,13,17,19,23,29,31,37,41,43,47]
    rows = {p:local_reduction_table(p) for p in primes}
    blocked = {n:[p for p in primes if n not in rows[p]] for n in range(6,11)}
    assert blocked == {
        6:[5], 7:[5,7,11,13,17,23], 8:[7,11,19,23],
        9:[5,7,11,13,17,23,29], 10:[5,7,13,17,19,29,37,47]}
    # For n>=7, primes 2 and 3 are independently forced by the cited local results.
    table = {n:[p for p in primes if n not in rows[p]] for n in range(7,35)}
    assert min(map(len,table.values())) == 4
    assert [n for n in table if len(table[n]) <= 4] == [8]
    extend = {n:[p for p in primes if n not in rows[p]] for n in range(7,67)}
    assert [n for n in extend if len(extend[n]) < 5] == [8]
    print("Primes that must be bad for a rational cycle:", blocked)
    print("For periods 7..66, at least 5 more bad primes each except 8, which has 4.")


if __name__ == "__main__":
    check_integer()
    check_param()
    check_synchronization()
    check_local()
