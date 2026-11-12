"""Exact support-cut Burnside intervals, with no floating-point arithmetic."""
from fractions import Fraction
from math import factorial, gcd, comb
from collections import Counter


def odd_nonunit_partitions(n, minimum=3):
    if n == 0:
        yield ()
    for j in range(minimum, n + 1, 2):
        for rest in odd_nonunit_partitions(n-j, j):
            yield (j,) + rest


def sector(n, mu):
    """Ordinary-class Burnside contribution including implicit fixed points."""
    m, r = sum(mu), len(mu)
    assert m <= n
    e = comb(n-m, 2) + (n-m)*r + sum((j-1)//2 for j in mu)
    e += sum(gcd(a,b) for i,a in enumerate(mu) for b in mu[i+1:])
    z = factorial(n-m)
    for j,c in Counter(mu).items():
        z *= j**c * factorial(c)
    return Fraction(1 << e, z)


def odd_derangement_weights(n):
    """c[m] = number of fixed-point-free odd-cycle permutations / m!."""
    c = [Fraction(1)]
    for m in range(1,n+1):
        c.append(sum((c[m-j] for j in range(3,m+1,2)),Fraction())/m)
    return c


def lower_and_tail(n, support_cut):
    """Return exact rational L,B such that L <= a(n) <= L+B."""
    assert n >= 0 and support_cut >= 0
    lower = sum((sector(n,mu) for m in range(min(n,support_cut)+1)
                 for mu in odd_nonunit_partitions(m)),Fraction())
    c = odd_derangement_weights(n)
    N = comb(n,2)
    tail = Fraction()
    for m in range(max(3,support_cut+1),n+1):
        d = (m*(2*n-m-1)+2)//3
        tail += c[m]*Fraction(1 << (N-d), factorial(n-m))
    return lower,tail


def integer_interval(n,support_cut):
    L,B = lower_and_tail(n,support_cut)
    U = L+B
    return -(-L.numerator//L.denominator), U.numerator//U.denominator


def sector_lower_log2(n):
    """A numerical display ONLY; proof uses the symbolic expression."""
    from math import log2
    return (n-(log2(n)+1)**2)/2

if __name__ == '__main__':
    for n in range(19):
        exact=lower_and_tail(n,n)[0]
        assert exact.denominator==1
        for s in (0,3,6,9):
            L,B=lower_and_tail(n,s)
            assert L<=exact<=L+B
    print('PASS: exact intervals at n=0..18, cuts 0,3,6,9')
