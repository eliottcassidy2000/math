"""Exact residue-power transforms and the first mod-four prime-race reversal.

No floating roots, asymptotic estimators, imported prime tables, or GRH gates.
Run with python -B [-O] and compare the JSON with the saved result.
"""
from bisect import bisect_right
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
import json

CAP = 1_000_000


def integer(n, low=0, high=None):
    if type(n) is not int or n < low or (high is not None and n > high):
        raise ValueError('exact integer outside the declared domain')
    return n


def iroot(x, k):
    integer(x)
    integer(k, 1)
    if x < 2 or k == 1:
        return x
    lo, hi = 0, 1 << ((x.bit_length()+k-1)//k)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**k <= x:
            lo = mid
        else:
            hi = mid
    return hi if hi**k <= x else lo


def mobius(n):
    integer(n, 1)
    sign, d = 1, 2
    while d*d <= n:
        if n % d == 0:
            n //= d
            sign = -sign
            if n % d == 0:
                return 0
            while n % d == 0:
                n //= d
        d += 1
    return -sign if n > 1 else sign


@lru_cache(maxsize=256)
def _primes(limit):
    marks = bytearray(b'\1')*(limit+1)
    marks[:min(2,limit+1)] = b'\0'*min(2,limit+1)
    for p in range(2,isqrt(limit)+1):
        if marks[p]:
            marks[p*p:limit+1:p] = b'\0'*((limit-p*p)//p+1)
    return tuple(i for i in range(2,limit+1) if marks[i])


def primes(limit):
    integer(limit, 0, CAP)
    return _primes(limit)


def prime_counts(x, modulus):
    integer(x, 0, CAP)
    integer(modulus, 2)
    out = [0]*modulus
    for p in primes(x):
        out[p % modulus] += 1
    return tuple(out)


def push_power(values, k):
    if type(values) is not tuple or len(values) < 2:
        raise ValueError('a full residue vector is required')
    if any(type(v) not in (int,Fraction) for v in values):
        raise ValueError('exact rational vector required')
    integer(k, 1)
    modulus = len(values)
    out = [Fraction(0)]*modulus
    for b, v in enumerate(values):
        out[pow(b,k,modulus)] += v
    return tuple(out)


def weighted_counts(x, modulus):
    """J_q(x,a): one weight 1/k for each prime power p^k <= x in class a."""
    integer(x, 0, CAP)
    integer(modulus, 2)
    out = [Fraction(0)]*modulus
    for p in primes(x):
        value, k = p, 1
        while value <= x:
            out[value % modulus] += Fraction(1,k)
            value *= p
            k += 1
    return tuple(out)


def forward_transform(x, modulus):
    integer(x, 0, CAP)
    integer(modulus, 2)
    out = [Fraction(0)]*modulus
    for k in range(1,x.bit_length()):
        row = push_power(prime_counts(iroot(x,k),modulus),k)
        out = [a+b/k for a,b in zip(out,row)]
    return tuple(out)


def inverse_transform(x, modulus):
    integer(x, 0, CAP)
    integer(modulus, 2)
    out = [Fraction(0)]*modulus
    for k in range(1,x.bit_length()):
        mu = mobius(k)
        if mu:
            row = push_power(weighted_counts(iroot(x,k),modulus),k)
            out = [a+Fraction(mu,k)*b for a,b in zip(out,row)]
    return tuple(out)


def mod4_split(x):
    integer(x, 0, CAP)
    counts = prime_counts(x,4)
    delta = counts[1]-counts[3]
    j = weighted_counts(x,4)
    squares = Fraction(max(0,len(primes(isqrt(x)))-int(x>=4)),2)
    remainder, absolute_bill = Fraction(0), Fraction(0)
    terms = []
    for k in range(3,x.bit_length()):
        root = iroot(x,k)
        row = prime_counts(root,4)
        numerator = row[1]-row[3] if k%2 else row[1]+row[3]
        value = Fraction(numerator,k)
        remainder += value
        absolute_bill += Fraction(row[1]+row[3],k)
        terms.append((k,root,value))
    return delta,j[1]-j[3],squares,remainder,absolute_bill,tuple(terms)


@dataclass(frozen=True)
class RaceReceipt:
    bound: int
    counts_one: int
    counts_three: int
    weighted_difference: Fraction
    square_correction: Fraction
    higher_power_correction: Fraction


def receipt(x):
    integer(x, 0, CAP)
    counts = prime_counts(x,4)
    _,j,s,r,_,_ = mod4_split(x)
    return RaceReceipt(x,counts[1],counts[3],j,s,r)


def audit_receipt(value, expected_bound):
    integer(expected_bound, 0, CAP)
    if type(value) is not RaceReceipt:
        raise ValueError('typed receipt required')
    for n in (value.bound,value.counts_one,value.counts_three):
        integer(n)
    for n in (value.weighted_difference,value.square_correction,value.higher_power_correction):
        if type(n) is not Fraction:
            raise ValueError('rational fields must remain exact')
    if value.bound != expected_bound or value != receipt(expected_bound):
        raise ValueError('changed source or arithmetic')
    return value


def first_lead(limit):
    integer(limit, 0, CAP)
    delta = 0
    for p in primes(limit):
        if p%4 == 1:
            delta += 1
        elif p%4 == 3:
            delta -= 1
        if delta > 0:
            return p
    return None  # only the supplied finite interval has been tested


def semiprime_split(x, modulus, negative_class):
    """Inherited real-character diagonal identity, for moduli 4 and 6."""
    integer(x, 0, CAP)
    if type(modulus) is not int or modulus not in (4,6):
        raise ValueError('only the stated two-class groups are supported')
    if type(negative_class) is not int or negative_class != modulus-1:
        raise ValueError('the retained negative class is required')
    ps = tuple(p for p in primes(x) if gcd(p,modulus)==1)
    signed_unordered = ordered = 0
    for i,p in enumerate(ps):
        if p*p > x:
            break
        for q in ps[i:]:
            if p*q > x:
                break
            signed_unordered += 1 if (p*q)%modulus==1 else -1
    for p in ps:
        for q in ps:
            if p*q > x:
                break
            ordered += (1 if p%modulus==1 else -1)*(1 if q%modulus==1 else -1)
    diagonal = sum(p*p<=x for p in ps)
    return signed_unordered,ordered,diagonal


CHECKS = 0


def check(truth, label):
    global CHECKS
    CHECKS += 1
    if not truth:
        raise RuntimeError(label)


def rejects(call):
    try:
        call()
    except (ValueError,TypeError):
        check(True,'exact rejection')
        return
    check(False,'accepted malformed input')


def main():
    # Independent primality path for every integer relevant to the first lead.
    ps = set(primes(30000))
    for n in range(30001):
        trial = n>=2 and all(n%d for d in range(2,isqrt(n)+1))
        check((n in ps)==trial,'sieve versus exhaustive trial division')
    for x in range(513):
        for k in range(1,11):
            r = iroot(x,k)
            check(r**k<=x<(r+1)**k,'integer-root endpoints')
        for modulus in range(2,13):
            check(forward_transform(x,modulus)==weighted_counts(x,modulus),'forward residue transform')
            check(inverse_transform(x,modulus)==prime_counts(x,modulus),'Mobius inverse with residue relabelling')
        delta,j,s,r,bill,_ = mod4_split(x)
        check(j==delta+s+r,'complete signed mod-four decomposition')
        check(abs(r)<=bill,'explicit higher-power interval')
    for modulus in range(2,13):
        values = tuple(Fraction(2*b-3,b+1) for b in range(modulus))
        for a in range(1,8):
            for b in range(1,8):
                check(push_power(push_power(values,a),b)==push_power(values,a*b),'power semigroup')
    for x in range(1,201):
        for modulus in (4,6):
            semis,ordered,diag = semiprime_split(x,modulus,modulus-1)
            check(2*semis==ordered+diag,'unordered diagonal identity')
    check(first_lead(26860) is None,'no earlier strict lead')
    check(first_lead(30000)==26861,'first strict lead')
    for x,one,three in ((26860,1472,1472),(26861,1473,1472),
                        (26862,1473,1472),(26863,1473,1473),(26879,1473,1474)):
        row=receipt(x)
        audit_receipt(row,x)
        check((row.counts_one,row.counts_three)==(one,three),'crossing counts')
    check(bisect_right(primes(30000),26861)==2946,'full prime index differs from class count')
    # The scalar inverse fails if the square map is replaced by the identity.
    wrong = weighted_counts(9,4)[1]-weighted_counts(3,4)[1]/2
    check(wrong==Fraction(3,2) and prime_counts(9,4)[1]==1,'discarded square preimage hostile')
    check(mod4_split(27)[3]==Fraction(-1,3),'squares alone are not the exact full correction')
    check(push_power((0,1,0,0),2)==push_power((0,0,0,1),2),'square map forgets the base class')
    check(prime_counts(9,4)[1]-prime_counts(9,4)[3] != weighted_counts(9,4)[1]-weighted_counts(9,4)[3],
          'prime count is not weighted prime-power count')
    # The mod-six race removes prime two from the mod-three negative class.
    for x in range(2,513):
        p3,p6 = prime_counts(x,3),prime_counts(x,6)
        check(p6[1]-p6[5]==p3[1]-p3[2]+1,'ramified-prime sidecar')
    sample=receipt(26861)
    for bad in (True,9.0,-1,CAP+1):
        rejects(lambda bad=bad: receipt(bad))
    rejects(lambda: audit_receipt(sample,26862))
    rejects(lambda: audit_receipt(replace(sample,counts_one=1473.0),26861))
    rejects(lambda: audit_receipt(replace(sample,bound=True),26861))
    rejects(lambda: push_power((0,True,0,0),2))
    rejects(lambda: iroot(9,False))
    rejects(lambda: first_lead(True))
    snapshots=[]
    for x in (9,27,100,1000,26861,26863,100000,1000000):
        delta,j,s,r,bill,terms=mod4_split(x)
        snapshots.append({'x':x,'pi_one':prime_counts(x,4)[1],'pi_three':prime_counts(x,4)[3],
                          'prime_difference':delta,'J_difference':str(j),'squares':str(s),
                          'higher_powers':str(r),'higher_absolute_bill':str(bill),
                          'higher_terms':[[k,root,str(value)] for k,root,value in terms]})
    result={'status':'PROVED elementary finite transforms; FINITE-EXACT first lead; analytic bias results separately CITED',
            'checks':CHECKS,'first_strict_mod4_lead':26861,'class_counts':[1473,1472],
            'full_prime_index':2946,'transform_universe':{'x':[0,512],'moduli':[2,12],'all_residue_classes':True},
            'trial_division_universe':[0,30000],'snapshots':snapshots,
            'hostiles':{'wrong_scalar_inverse_at9':str(wrong),'higher_power_at27':'-1/3',
                        'square_class_collision':[1,3],'first_lead_is_not_permanent':26863},
            'scope':'No finite census proves a limiting density or an all-source Collatz rule.'}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
