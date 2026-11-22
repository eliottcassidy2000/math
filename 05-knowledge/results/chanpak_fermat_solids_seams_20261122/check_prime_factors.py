"""Lucas certificates, exact order bounds, and a literal realizable prime factor."""
from math import factorial, gcd, isqrt, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_prime_frontier_selector_20261118'))
from join_loads import analyze
from small_arm import g

SMALL = 2432582681
LARGE = 1092853292237112554142488617
MID = 1468888833652033002879689
# Each row lists a full-order residue and the complete prime factors of n-1.
LUCAS = {
    SMALL: (3, {2: 3, 5: 1, 11: 1, 2017: 1, 2741: 1}),
    MID: (3, {2: 3, 59: 1, 137: 1, 3557: 1, 3761: 1, 155663: 1, 10908217: 1}),
    LARGE: (5, {2: 3, 3: 1, 31: 1, MID: 1}),
}


def prime_certificate(n):
    if n not in LUCAS:
        assert 2 <= n <= 10908217
        assert all(n % d for d in range(2, isqrt(n) + 1))
        return
    a, factors = LUCAS[n]
    assert prod(q**k for q, k in factors.items()) == n - 1
    assert pow(a, n-1, n) == 1
    for q in factors:
        prime_certificate(q)
        assert gcd(pow(a, (n-1)//q, n)-1, n) == 1


def blocks_bound(widths):
    maxima = {1: 1, 5: 15, 10: 15745}
    return factorial(sum(widths)) * prod(maxima[w] for w in widths) // prod(factorial(w) for w in widths)


def reachable(rows):
    vertices = 1
    while True:
        old = vertices
        for i, targets in enumerate(rows):
            if old >> i & 1:
                vertices |= targets
        if old == vertices:
            return vertices.bit_count()


def main():
    assert SMALL * LARGE == (1 << 60) * ((1 << 61)-1) + 1
    prime_certificate(SMALL)
    prime_certificate(LARGE)
    for p, cap, bounds in ((SMALL, 41, (10,5)), (LARGE, 117, (10,10,10,1))):
        assert 5**(cap-1) <= p**3 < 5**cap
        assert blocks_bound(bounds) < p
        print('prime', p, 'strong order range', sum(bounds)+1, cap)
    b, e = 49, 41466580
    witness = analyze(b,e)
    assert witness['H'] == 47*8525963 + b*e == SMALL
    assert (g(b-1),g(e-1)) == (47,8525963)
    rows = witness['rows']
    reverse = [sum(1 << j for j,r in enumerate(rows) if r >> i & 1) for i in range(len(rows))]
    assert len(rows) == reachable(rows) == reachable(reverse) == 34
    print('small prime strong F1 witness', b, e, 'on', len(rows), 'vertices')


if __name__ == '__main__':
    main()
