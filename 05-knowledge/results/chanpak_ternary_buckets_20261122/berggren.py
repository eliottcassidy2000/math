"""Berggren decoding, its two odd labels, and finite bases for complete count rays."""
from collections import Counter
from math import gcd, isqrt
from pathlib import Path
import sys

PACKETS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKETS / 'chanpak_decimal_kernel_20261122'))
from e15_pair import literal_count, star_rows, subpaths


def child(m, n, digit):
    return ((m+2*n, n), (2*m+n, m), (2*m-n, m))[digit]


def parameters(h):
    assert h > 0 and h & 1
    k = (h-1)//2
    digits = []
    while k:
        k, digit = divmod(k, 3)
        digits.append(digit)
    m, n = 1, 0
    for digit in reversed(digits):
        m, n = child(m, n, digit)
    return m, n


def labels(m, n):
    digits = []
    while (m, n) != (2, 1):
        if m > 3*n:
            digits.append(0)
            m, n = m-2*n, n
        elif m > 2*n:
            digits.append(1)
            m, n = n, m-2*n
        else:
            digits.append(2)
            m, n = n, 2*n-m
    low = sum(digit*3**i for i,digit in enumerate(digits))
    scale = 3**len(digits)
    return (2*(scale+low)+1, 2*(2*scale+low)+1)


SMALL = {11:(8,5,9,18,7), 13:(8,5,17,6,11), 17:(2,4,9,35,15,23),
         63:(2,16,3,23,197,31,47,111)}

# (core row masks, core apex mask, left address, right address).
LITERAL = {
    1549681951: ((2,8,3,5), 5, 91820640, 0),
    2324522929: ((0,), 0, 1166118272, 0),
    6973568801: ((0,), 0, 3955238204, 0),
    20920706405: ((0,), 0, 10625241056, 0),
    188286357653: ((2,4,1), 3, 6449010978, 2),
    376572715303: ((2,8,3,5), 6, 26073287319, 0),
    564859072961: ((0,), 0, 501719569403, 0),
}


def factor(n):
    factors = []
    p = 3
    while p*p <= n:
        while n%p == 0:
            factors.append(p)
            n //= p
        p += 2
    if n != 1:
        factors.append(n)
    return factors


def base_products(h):
    if h <= 10**9:
        return [h]
    primes = factor(h)
    sevens = primes.count(7)
    other = [p for p in primes if p != 7]
    if sevens == 1:
        # Combine the seven with a realized factor; 21 itself is forbidden.
        if other.count(3) >= 2:
            other.remove(3)
            other.remove(3)
            other.append(63)
        elif any(p != 3 and 7*p <= 10**9 for p in other):
            p = next(p for p in other if p != 3 and 7*p <= 10**9)
            other.remove(p)
            other.append(7*p)
        else:
            return []
    if sevens >= 2:
        if sevens%2:
            other.append(343)
            sevens -= 3
        other.extend([49]*(sevens//2))
    return other if all(p <= 10**9 for p in other) else []


def main():
    for h, rows in SMALL.items():
        assert subpaths(rows)[-1] == h
    for h,(rows, tag, x, y) in LITERAL.items():
        host = star_rows(rows, tag, x, y)
        sizes = (1,)*x.bit_length() + (len(rows),) + (1,)*y.bit_length()
        assert literal_count(host, sizes) == h
    for h in range(3,5000,2):
        pair = parameters(h)
        assert h in labels(*pair)
        m,n = pair
        a,b = m*m-n*n, 2*m*n
        for digit, delta in enumerate((a+b, b-a, -a-b)):
            M,N = child(m,n,digit)
            assert M*M-N*N-2*M*N == delta
    by_gap = {1: [3,5]}
    cap = 1200
    for m in range(2,isqrt(cap)+1):
        for n in range(1,m):
            if gcd(m,n)>1 or (m-n)%2==0:
                continue
            gap = m*m+2*m*n-n*n
            if gap <= cap:
                bucket = by_gap.setdefault(gap,[])
                for digit in (0,2):
                    bucket.extend(labels(*child(m,n,digit)))
    assert sorted(by_gap[7]) == [7,11,13,17]
    tally=Counter()
    for gap, seeds in by_gap.items():
        assert len(set(seeds)) == len(seeds) == 2**(1+len(set(factor(gap))))
        for h in seeds:
            if h in (7,21):
                tally['forbidden seed'] += 1
            elif h in LITERAL:
                tally['literal'] += 1
            else:
                fs=base_products(h)
                assert fs and all(p not in (7,21) for p in fs), (gap,h)
                value=1
                for p in fs:
                    value *= p
                assert value==h
                tally['finite window' if h<=10**9 else 'ordinal product'] += 1
    print('inclusive gap cap',cap,'distinct gaps',len(by_gap),'scalar seeds',sum(map(len,by_gap.values())),dict(tally))


if __name__ == '__main__':
    main()
