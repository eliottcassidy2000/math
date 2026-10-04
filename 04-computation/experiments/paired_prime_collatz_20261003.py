#!/usr/bin/env python3
"""Exact paired prime sieve and Collatz identities; no sampled or inherited filters."""
from collections import Counter
from itertools import product
from math import factorial, isqrt

LIMIT = 100_000


def prime_table(n):
    p = bytearray(b"\1") * (n + 1)
    p[0:2] = b"\0\0"
    for d in range(2, isqrt(n) + 1):
        if p[d]:
            start = d*d
            p[start::d] = b"\0" * ((n-start)//d+1)
    return p


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def star(a, b):
    return a+b+2*a*b


def shortcut(n, sign):
    return n//2 if n % 2 == 0 else (3*n+sign)//2


def ordinary(n, sign):
    return n//2 if n % 2 == 0 else 3*n+sign


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def main():
    primes = prime_table(6*LIMIT+1)
    assert all(bool(primes[n]) == trial_prime(n) for n in range(2, 40_002))
    print("Eratosthenes vs independent trial division: integers 2..40001")
    # Complete factor universe: every composite has a factor <= sqrt(maximum).
    CA, CB = set(), set()
    witnesses = 0
    for u in range(3, isqrt(4*LIMIT+1)+1, 2):
        for v in range(u, (4*LIMIT+1)//u+1, 2):
            s, t = (u+1)//4, (v+1)//4
            if u % 4 == v % 4 == 1:
                N, family = 4*s*t+s+t, CA
            elif u % 4 == v % 4 == 3:
                N, family = 4*s*t-s-t, CA
            elif u % 4 == 1:
                N, family = 4*s*t-s+t, CB
            else:
                N, family = 4*s*t+s-t, CB
            assert 4*N+(1 if family is CA else -1) == u*v
            if N <= LIMIT:
                family.add(N)
            witnesses += 1
    expected_A = {N for N in range(1, LIMIT+1) if not primes[4*N+1]}
    expected_B = {N for N in range(1, LIMIT+1) if not primes[4*N-1]}
    assert CA == expected_A and CB == expected_B
    print(f"Four family product laws exactly classify {2*LIMIT} labels; {witnesses} factor pairs")
    PA, PB = sorted(set(range(1,LIMIT+1))-CA), sorted(set(range(1,LIMIT+1))-CB)
    print("A prime indices first 20:", PA[:20])
    print("B prime indices first 20:", PB[:20])
    def gap_record(xs):
        return max((b-a, a, b) for a,b in zip(xs,xs[1:]))
    assert gap_record(sorted(CA))[0] == gap_record(sorted(CB))[0] == 3
    print("Largest consecutive composite-index gap in each row: 3 (sharp)")
    for label, xs in (("A",PA),("B",PB),("A union B", sorted(set(PA)|set(PB)))):
        print(label, "prime-index count / max gap (gap,left,right):", len(xs), gap_record(xs))
    states=Counter((bool(primes[4*N-1]),bool(primes[4*N+1])) for N in range(1,LIMIT+1))
    print("(B prime,A prime) state counts:", sorted(states.items()))
    both_composite = CA & CB
    assert min(both_composite) == 14
    print("First both-composite atom: N=14, (55,57)")
    run = best = 0
    start = best_start = 0
    for N in range(1,LIMIT+1):
        if N in both_composite:
            if run == 0:
                start = N
            run += 1
            if run > best:
                best, best_start = run, start
        else:
            run = 0
    print("Longest both-composite run in 1..100000:", (best_start,best_start+best-1,best))
    # Infinite gap construction: each listed proper divisor is explicit.
    for L in range(1, 21):
        K = factorial(4*L+2)
        for j in range(1,L+1):
            for sign in (-1,1):
                d = 4*j+sign
                value = K+d
                assert value % d == 0 and 1 < d < value
    print("Factorial simultaneous-gap witnesses: lengths 1..20, both rows")
    # Algebraic controls, including both signs, parity branches and n=1.
    for N in range(1, 10_001):
        assert [shortcut(2*N-1,1),shortcut(2*N,1)] == [3*N-1,N]
        assert [shortcut(2*N,-1),shortcut(2*N+1,-1)] == [N,3*N+1]
        for sign in (-1,1):
            phi = lambda n: 2*n+sign
            m = phi(N)
            transported = (m+sign)//2 if N % 2 == 0 else 3*m
            assert phi(ordinary(N,sign)) == transported
        assert odd_part(3*(4*N-1)+1) == 6*N-1
        assert odd_part(3*(4*N+1)-1) == 6*N+1
        assert odd_part(3*(4*N-1)-1) == odd_part(3*N-1)
        assert odd_part(3*(4*N+1)+1) == odd_part(3*N+1)
        # Multiplication by 3 in the z chart is the 3n+1 affine branch.
        assert star(1,N) == 3*N+1
        assert 2*star(1,N)+1 == 3*(2*N+1)
    print("Central-block shortcut pairs and both affine conjugacies: n=1..10000")
    assert [odd_part(3*n-1) for n in (1,5,7)] == [1,7,5]
    assert [odd_part(3*n+1) for n in (3,7,11,15)] == [5,11,17,23]
    # Neither sign descends to y alone, even with all outputs >=3.
    plus = [odd_part(3*m+1) for m in (7,9)]
    minus = [odd_part(3*m-1) for m in (7,9)]
    assert plus == [11,7] and minus == [5,13]
    assert [(m+1)//4 for m in plus] == [3,2]
    assert [(m+1)//4 for m in minus] == [1,3]
    print("Same-atom full-odd dynamics: + maps (7,9)->(11,7); - maps (7,9)->(5,13)")
    assert primes[7] and primes[23] and primes[47]
    assert odd_part(3*7+1) == 11 and odd_part(3*23+1) == 35
    assert odd_part(3*7-1) == 5 and odd_part(3*47-1) == 35
    print("Family+prime flag fails to determine next prime flag: + at 7/23, - at 7/47")
    print("3n-1 cycle control: odd accelerated cycle 5->7->5; 1 fixed")
    for a,b,c,d in product(range(6),repeat=4):
        defect = star(a+b,c+d)-star(a,c)-star(b,d)
        assert defect == 2*(a*d+b*c)
    assert star(2,2) == 12 and 2*star(1,1) == 8
    print("Eckmann-Hilton interchange defect 2(ad+bc): all 1296 inputs in 0..5")
    for I in range(1,51):
        for J in range(I+1,52):
            O,E = 2*(I+J)-1,4*J
            assert (O+1)//2 == I+J and E//2 == 2*J
            assert (E//2)-(O+1)//2 == J-I
    print("q atom-address transport: 1275 positive pairs with I<J<=51")
    # Prime-pair controls for the two different twin-prime placements.
    assert primes[11] and primes[13] and primes[5] and primes[7]
    vertical=sum(primes[4*N-1] and primes[4*N+1] for N in range(1,LIMIT+1))
    diagonal=sum(primes[4*N+1] and primes[4*N+3] for N in range(1,LIMIT))
    print("Twin-pair counts (vertical, diagonal within atom window):", (vertical,diagonal))
    quartets = [N for N in range(1,LIMIT+1)
                if all(primes[c*N+s] for c in (4,6) for s in (-1,1))]
    independent_quartets = [N for N in range(1,LIMIT+1)
                            if all(trial_prime(c*N+s) for c in (4,6) for s in (-1,1))]
    assert quartets == independent_quartets
    print("Four-label prime atoms (independent trial-division replay): count / first 20:", len(quartets), quartets[:20])
    for p in range(2,102):
        if primes[p]:
            roots = {N for N in range(p) if any((c*N+s)%p == 0
                      for c in (4,6) for s in (-1,1))}
            expected = 0 if p == 2 else 2 if p in (3,5) else 4
            assert len(roots) == expected and 0 not in roots
    print("Four-form local root counts checked at every prime <=101: 0,2,2,4 thereafter")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
