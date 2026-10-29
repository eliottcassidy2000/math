#!/usr/bin/env python3
"""Standalone exact finite audit; no assert, no producer imports."""
import itertools
import json
from fractions import Fraction
from math import gcd, isqrt, lcm


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def prime(p):
    return p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))


def factors(n):
    out = []
    q = 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def order(a, p, k):
    require(prime(p) and gcd(a, p) == 1, 'order domain')
    m = p ** k
    d = (p - 1) * p ** (k - 1)
    for q in factors(d):
        while d % q == 0 and pow(a, d // q, m) == 1:
            d //= q
    require(pow(a, d, m) == 1, 'order return')
    require(all(pow(a, d // q, m) != 1 for q in factors(d)), 'minimality')
    return d


def valuation(n, p):
    require(n != 0, 'finite valuation')
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def fq(a, p):
    z = pow(a % (p * p), p - 1, p * p)
    require(z % p == 1, 'Fermat little theorem')
    return ((z - 1) // p) % p


def native(n, word, sign):
    for a in word:
        z = 3 * n + sign
        require(z > 0 and valuation(z, 2) == a, 'native exact guard')
        n = z // (2 ** a)
    return n


def main():
    ps = [p for p in range(3, 44) if prime(p)] + [1093, 3511]
    rows = []
    lift_checks = 0
    for p in ps:
        for a in (2, 3, 5):
            if a % p == 0:
                continue
            d = order(a, p, 1)
            s = valuation(pow(a, d) - 1, p)
            for k in range(1, 5):
                require(order(a, p, k) == d * p ** max(0, k - s), 'lift formula')
                lift_checks += 1
            if p in (1093, 3511) and a in (2, 3):
                rows.append(dict(p=p, base=a, order=d, defect=s,
                                 quotient=fq(a, p), square_order=order(a, p, 2)))
    expected = [(1093, 2, 364, 2, 0, 364), (1093, 3, 7, 1, 312, 7651),
                (3511, 2, 1755, 2, 0, 1755), (3511, 3, 270, 1, 7, 947970)]
    require([tuple(r[k] for k in ('p','base','order','defect','quotient','square_order'))
             for r in rows] == expected, 'published exact rows')
    joint_subgroups = []
    for p, expected_h1 in ((1093,364),(3511,3510)):
        d2,d3 = order(2,p,1),order(3,p,1)
        H1 = {pow(2,i,p)*pow(3,j,p)%p for i in range(d2) for j in range(d3)}
        require(len(H1)==expected_h1==lcm(d2,d3), 'joint subgroup base')
        for k in range(1,4):
            # Generator of the full kernel of reduction p^(k+1) -> p^k.
            g = pow(3,order(3,p,k),p**(k+1))
            require(g!=1 and g%p**k==1 and pow(g,p,p**(k+1))==1,
                    'joint subgroup full kernel')
        joint_subgroups.append(dict(p=p,level_one=len(H1),level_two=p*len(H1)))
    clock_checks = 0
    for p in ps:
        if p == 3:
            continue
        m = p * p
        for J in range(13):
            for D in range(13):
                alpha = pow(3, J, m) * pow(pow(2, D, m), -1, m) % m
                require(fq(alpha, p) == (J * fq(3, p) - D * fq(2, p)) % p,
                        'two clock quotient')
                if alpha % p == 1:
                    require(fq(alpha, p) == -((alpha - 1) // p) % p, 'identity lift')
                clock_checks += 1
    # Actual complete exponent periods, not presumed independent CRT coordinates.
    period_checks = 0
    for p in (3, 5, 7, 11, 1093, 3511):
        for k in (1, 2):
            L = order(2, p, k)
            xs = [pow(2, t, p ** k) for t in range(L)]
            require(len(set(xs)) == L, 'period injectivity')
            for h in (0, 1, L - 1, L, 2 * L):
                hits = sum(xs[(t + h) % L] == xs[t] for t in range(L))
                require(hits == (L if h % L == 0 else 0), 'lag equality correlation')
            period_checks += 1
    for m1, m2 in ((9,25),(49,121)):
        L1 = next(t for t in range(1,m1+1) if pow(2,t,m1)==1)
        L2 = next(t for t in range(1,m2+1) if pow(2,t,m2)==1)
        orbit = {(pow(2,t,m1),pow(2,t,m2)) for t in range(lcm(L1,L2))}
        require(len(orbit)==lcm(L1,L2), 'CRT lcm period')
    word_checks = source_checks = 0
    for sign in (-1, 1):
        for J in range(1, 5):
            for word in itertools.product(range(1, 4), repeat=J):
                C, D = 0, 0
                for a in word:
                    C = 3 * C + 2 ** D
                    D += a
                A, B = 3 ** J, 2 ** D
                M = 2 * B
                r = ((B - sign * C) * pow(A, -1, M)) % M
                require(r % 2 == 1 and 0 < r < M, 'canonical odd source')
                require((A*r+sign*C) % B == 0, 'integral endpoint')
                h = (A*r+sign*C)//B
                for t in range(12):
                    n = r + M*t
                    y = native(n,word,sign)
                    require(y == h + 2*A*t, 'chart endpoint')
                    source_checks += 1
                for m in (5,7,25,35,49):
                    ns = [(r+M*t)%m for t in range(m)]
                    ys = [(h+2*A*t)%m for t in range(m)]
                    require(len(set(ns)) == m and len(set(ys)) == m, 'uniform marginals')
                    require(len(set(zip(ns,ys))) == m, 'graph support')
                    for N in (1,m-1,m,m+1,2*m+3):
                        counts = [0]*m
                        for t in range(3,3+N):
                            counts[(h+2*A*t)%m] += 1
                        require(all(c in (N//m,(N+m-1)//m) for c in counts), 'interval bound')
                for p in (1093,3511):
                    require(gcd(M,p*p)==gcd(2*A,p*p)==1, 'Wieferich CRT permutation certificate')
                word_checks += 1
    affine_checks = 0
    for p in (3,5,7):
        m = p*p
        for u,v,n in itertools.product(range(p),range(p),range(m)):
            alpha,beta = 1+p*u,p*v
            z=n
            first=None
            for t in range(1,p+1):
                z=(alpha*z+beta)%m
                require(z==(n+t*p*(u*n+v))%m,'affine iterate')
                if z==n and first is None:
                    first=t
            require(first==(1 if (u*n+v)%p==0 else p),'affine period')
            affine_checks += 1
    return dict(status='PASS',order_lift_checks=lift_checks,numerical_rows=rows,
                joint_subgroups=joint_subgroups,
                clock_checks=clock_checks,complete_period_checks=period_checks,
                signed_exact_words=word_checks,native_sources=source_checks,
                affine_point_checks=affine_checks,
                scope='finite checks only; all-parameter proofs in ORDERS_CARRIES_AND_CRT.md')


if __name__ == '__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
