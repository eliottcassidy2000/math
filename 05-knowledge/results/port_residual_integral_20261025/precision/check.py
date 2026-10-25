"""Fixed 2-adic precision from bounded odd-cycle partitions.

The production path uses odd denominator inverses, not rational partitions.
Full Fraction sums appear only in the independent bounded verification.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import factorial, gcd, isqrt, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def v2(n):
    need(n > 0, 'v2 requires positive integer')
    return (n & -n).bit_length()-1


def factorial_v2(m):
    return m-m.bit_count()


def baseline(n):
    return 0 if n <= 2 else n//2


def beta(m, odd):
    if odd:
        return (m-1)**2//2-factorial_v2(m)
    if m == 2:
        return 0  # separately proved two-cycle bound for even n>=4
    return m*(m-2)//2-factorial_v2(m)


def cycle_counts(n, precision):
    cap = min(n,2+isqrt(2*precision+2))
    return [m for m in range(1 if n%2 else 2,cap+1,2)
            if beta(m,bool(n%2)) < precision]


def fixed_length_partitions(total, count, least=1):
    if count == 0:
        if total == 0:
            yield ()
        return
    if total < count*least or total%2 != count%2:
        return
    if count == 1:
        if total >= least and total%2:
            yield (total,)
        return
    for first in range(least,total//count+1,2):
        for tail in fixed_length_partitions(total-first,count-1,first):
            yield (first,)+tail


def term(parts):
    mult = Counter(parts)
    n,m = sum(parts),len(parts)
    G = sum(gcd(a,b) for a,b in combinations(parts,2))
    exponent = (n-m)//2+G
    z = prod(d**k*factorial(k) for d,k in mult.items())
    twos = sum(factorial_v2(k) for k in mult.values())
    return exponent,z,twos


def normalized_residue(n, precision):
    if type(n) is not int or n < 0 or type(precision) is not int or precision < 1:
        raise ValueError('n>=0 and precision>=1 must be integers')
    modulus = 1 << precision
    if n <= 2:
        return 1 % modulus, {'cycle_counts':[], 'visited':0, 'active':0}
    counts = cycle_counts(n,precision)
    answer = visited = active = 0
    B = baseline(n)
    for m in counts:
        for parts in fixed_length_partitions(n,m):
            visited += 1
            exponent,z,twos = term(parts)
            relative = exponent-twos-B
            need(relative >= 0, ('baseline',n,parts))
            if relative >= precision:
                continue
            odd = z >> twos
            need(odd%2 == 1, 'denominator must be a unit')
            answer = (answer + (1 << relative)*pow(odd,-1,modulus)) % modulus
            active += 1
    return answer, {'cycle_counts':counts,'visited':visited,'active':active}


def all_odd_partitions(n, least=1):
    if n == 0:
        yield ()
    for first in range(least,n+1,2):
        for tail in all_odd_partitions(n-first,first):
            yield (first,)+tail


def totient(n):
    return sum(gcd(a,n) == 1 for a in range(1,n+1))


def leading_class(n):
    if n == 4:
        return True
    if n%2 or n%4 == 0:
        return False
    m = n//2
    p = next((d for d in range(3,isqrt(m)+1,2) if m%d == 0),m)
    while m%p == 0:
        m //= p
    return m == 1 and p%4 == 3


def main():
    total_partitions = coefficient_checks = 0
    samples = []
    known = [1,1,1,2,4,12,56,456,6880,191536,9733056,903753248]
    for n in range(51):
        exact = Fraction(0)
        partitions = list(all_odd_partitions(n))
        total_partitions += len(partitions)
        for parts in partitions:
            e,z,t = term(parts)
            exact += Fraction(1 << e,z)
            if n >= 3:
                relative = e-t-baseline(n)
                need(relative >= 0, ('termwise baseline',n,parts))
                need(relative >= beta(len(parts),bool(n%2)), ('cycle bound',n,parts))
        need(exact.denominator == 1, ('integer Burnside count',n))
        exact = exact.numerator
        if n < len(known):
            need(exact == known[n], ('inherited prefix',n))
        need(exact % (1 << baseline(n)) == 0, ('whole baseline',n))
        normalized = exact >> baseline(n)
        if n >= 3 and n%2:
            need(normalized%2 == 1, ('inherited odd valuation',n))
            need(normalized%4 == (n+2)%4, ('odd second bit',n))
        elif n >= 4:
            need(normalized%2 == (totient(n)//2)%2, ('even leading bit',n))
            need((v2(exact) == n//2) == leading_class(n), ('even equality classification',n))
            predicted = (-1)**(n//2+1)*(totient(n)//2)+2*int(n in (4,6))
            need(normalized%4 == predicted%4, ('even second bit',n))
        for K in range(1,7):
            residue,stats = normalized_residue(n,K)
            need(residue == normalized % (1 << K), ('precision comparison',n,K))
            allowed = set(stats['cycle_counts'])
            if n >= 3:
                for parts in partitions:
                    if len(parts) not in allowed:
                        e,z,t = term(parts)
                        need(e-t >= baseline(n)+K, ('omitted contribution',n,K,parts))
            coefficient_checks += 1
            if n in (20,30,40,49,50) and K in (1,2,6):
                samples.append(dict(n=n,K=K,all_partitions=len(partitions),
                                    normalized_residue=residue,**stats))
    for K in range(1,101):
        cap = 2+isqrt(2*K+2)
        for m in range(1,401):
            if beta(m,bool(m%2)) < K:
                need(m <= cap, ('uniform square-root budget',K,m))
    e,z,t = term((3,3))
    need(e-t == 4 and baseline(6) == 3, 'equal-cycle hostile')
    e,z,t = term((1,1))
    need(e-t == 0, 'order-two exception')
    for bad in ((-1,1),(3,0),(True,1),(3,1.5)):
        try:
            normalized_residue(*bad)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid input accepted')
    print(json.dumps({'status':'PASS','orders':[0,50],'precision_bits':[1,6],
        'modular_coefficient_checks':coefficient_checks,'full_odd_partitions':total_partitions,
        'samples':samples,'cycle_budgets_K1_to_K6':{
            'odd':[cycle_counts(101,K) for K in range(1,7)],
            'even':[cycle_counts(100,K) for K in range(1,7)]},
        'scope':'low-precision Burnside filter, no enumeration of tournament cores'},
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
