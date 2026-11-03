#!/usr/bin/env python3
"""Finite algebra controls for FI_INTERFACES.md; proofs are in the note."""
from itertools import product, permutations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def clean(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    return clean([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                  for i in range(max(len(a), len(b)))])


def mul(a, b):
    out = [0] * max(0, len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return clean(out)


def star(a, b):
    out = [0] * max(0, len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            for t in range(min(i, j)+1):
                out[i+j-t] += x*y*comb(i,t)*comb(j,t)*factorial(t)
    return clean(out)


def deriv(a):
    return clean([i*a[i] for i in range(1, len(a))])


def shift(a, k):
    for _ in range(k):
        a = add(a, deriv(a))
    return a


def fall(n, k):
    out = 1
    for j in range(k):
        out *= n-j
    return out


def evaluate(a, n):
    return sum(x*fall(n,i) for i, x in enumerate(a))


def hom(a, b):
    out = (1,)
    for k, exponent in enumerate(a):
        for _ in range(exponent):
            out = star(out, shift(b,k))
    return out


def nat_count(a, b):
    out = 1
    for k, exponent in enumerate(a):
        out *= evaluate(b,k)**exponent
    return out


def main():
    profiles = [clean(p) for p in product(range(2), repeat=3)]
    pair_checks = 0
    adjunction_checks = 0
    for a, b in product(profiles, repeat=2):
        require(shift(star(a,b),1) == star(shift(a,1),shift(b,1)),
                ('shift is a star homomorphism', a,b))
        require(deriv(star(a,b)) == add(add(star(deriv(a),b), star(a,deriv(b))),
                                      star(deriv(a),deriv(b))),
                ('twisted Leibniz', a,b))
        for r in range(7):
            expected = 1
            for k, exponent in enumerate(a):
                expected *= evaluate(b,r+k)**exponent
            require(evaluate(hom(a,b),r) == expected, ('hom formula',a,b,r))
            require(evaluate(star(a,b),r) == evaluate(a,r)*evaluate(b,r),
                    ('matching evaluation',a,b,r))
        require(evaluate(hom(a,b),0) == nat_count(a,b), ('map count',a,b))
        pair_checks += 1
        for c in profiles:
            require(hom(mul(a,b),c) == hom(a,hom(b,c)), ('currying',a,b,c))
            require(nat_count(mul(a,b),c) == nat_count(a,hom(b,c)),
                    ('adjunction cardinality',a,b,c))
            adjunction_checks += 1
    # Independent literal injection enumeration of the matching decomposition.
    injection_cases = 0
    for n in range(6):
        for k,l in product(range(4), repeat=2):
            actual = [0]*(min(k,l)+1)
            for left in permutations(range(n), k):
                for right in permutations(range(n), l):
                    actual[len(set(left)&set(right))] += 1
            for t, count in enumerate(actual):
                expected = comb(k,t)*comb(l,t)*factorial(t)*fall(n,k+l-t)
                require(count == expected, ('literal matching',n,k,l,t))
            injection_cases += 1
    # Independent pointed-complement decomposition: different old or new slots.
    pointer_cases = 0
    for n in range(6):
        for k in range(1,5):
            actual = sum(1 for inj in permutations(range(n),k)
                         for slot in range(n) if slot != inj[0])
            expected = (k-1)*fall(n,k)+fall(n,k+1)
            require(actual == expected, ('pointed complement',n,k))
            pointer_cases += 1
    require(hom((0,1),(0,3,3,1)) == (3,9,6,1), 'triangle hom example')
    print({'profile_pairs':pair_checks, 'currying_and_adjunction_triples':adjunction_checks,
           'literal_matching_cases':injection_cases, 'pointed_complement_cases':pointer_cases,
           'status':'PASS'})


if __name__ == '__main__':
    main()
