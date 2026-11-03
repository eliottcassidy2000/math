#!/usr/bin/env python3
"""Exact prime clocks and CRT obstructions; no stochastic assertions."""
import bisect
import json
import math
from fractions import Fraction


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def sieve(limit):
    flags = bytearray(b"\1") * (limit + 1)
    flags[:2] = b"\0\0"
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b"\0" * ((limit-p*p)//p+1)
    return flags, [p for p in range(limit + 1) if flags[p]]


def factors(n):
    out = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1
    if n > 1:
        out.append((n, 1))
    return out


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def first_hit(m, flags):
    k = 1
    while k*m+1 < len(flags):
        for p in (k*m-1, k*m+1):
            if flags[p]:
                return p
        k += 1
    raise RuntimeError(f"sieve too small for {m}")


def main():
    flags, primes = sieve(4_000_000)
    table = []
    rows = []
    record = Fraction(0)
    records = []
    expected = {262:(54,263),278:(57,277),292:(60,293),
                120:(49,239),144:(81,431),208:(173,1039),290:(227,1451)}
    for m in range(4, 10001, 2):
        fs = factors(m)
        phi = m
        for p, _ in fs:
            phi = phi//p*(p-1)
        hit = first_hit(m, flags)
        rank = bisect.bisect_right(primes, hit)-len(fs)
        dead = bisect.bisect_right(primes, m-2)-len(fs)
        ratio = Fraction(2*rank, phi)
        row = dict(m=m, phi=phi, first_prime=hit, eligible_rank=rank,
                   zero_hazard_prefix=dead, residual_prime_rank=rank-dead,
                   iid_ratio=str(ratio))
        rows.append(row)
        if ratio > record:
            records.append(row)
            record = ratio
        if m in expected or m in (210,9970):
            direct = [p for p in range(2,hit+1)
                      if trial_prime(p) and math.gcd(p,m)==1]
            check(len(direct)==rank, f"independent rank: {m}")
            check(next(p for p in direct if p%m in (1,m-1))==hit,
                  f"independent hit: {m}")
            if m in expected:
                check((rank,hit)==expected[m], f"pasted table: {m}")
            table.append(row)
    small_divisors = [17,3,3,7,11,13,19,3,3]
    candidates = [k*290+s for k in range(1,6) for s in (-1,1)][:-1]
    obstruction = []
    M = math.prod([2,3,7,11,13,17,19])
    for i,(n,q) in enumerate(zip(candidates,small_divisors)):
        check(n>q and n%q==0 and M%q==0, "CRT certificate")
        obstruction.append(dict(candidate=n, divisor=q, factorization=factors(n)))
    # Independent finite probes of the entire infinite CRT certificate's pattern.
    for t in (0,1,2,17,10**6):
        m = 290+M*t
        for i,q in enumerate(small_divisors):
            k,s = i//2+1, (-1 if i%2==0 else 1)
            check((k*m+s)%q==0 and k*m+s>q, "progression hostile")
    units = [a for a in range(290) if math.gcd(a,290)==1]
    reps = {min(a,290-a) for a in units}
    def order(a):
        x = a
        for k in range(1,113):
            if x in (1,289):
                return k
            x = x*a%290
        raise RuntimeError("quotient order")
    orders = {}
    for a in reps:
        k = order(a)
        orders[k] = orders.get(k,0)+1
    check(len(reps)==56 and max(orders)==28 and orders.get(2)==3,
          "C2 x C28 quotient invariants")
    stern = [(a,227-2*a*a,factors(227-2*a*a))
             for a in range(1,math.isqrt(226//2)+1)]
    check(all(not trial_prime(n) for _,n,_ in stern), "227 Stern control")
    check(all(any((227-2*a*a)%q==0 for q in (3,5,11))
              for a in range(1,11)), "Stern small-prime cover")
    check(trial_prime(557) and trial_prime(269) and 557-2*12**2==269,
          "Stern continuation hostile")
    check(all(not trial_prime(557-2*a*a) for a in range(1,12)),
          "Stern first new success")
    # Exact means of the MRB coefficient masks, retaining parity.
    mask_means = {
        'unit_alternating': Fraction(sum((-1)**n for n in range(1,291)
                                       if math.gcd(n,290)==1),290),
        'plus_minus_alternating': Fraction(sum((-1)**n for n in range(1,291)
                                             if n%290 in (1,289)),290),
        'oriented_plus_minus': Fraction(sum((n%290==1)-(n%290==289)
                                           for n in range(1,291)),290)}
    check(mask_means==dict(unit_alternating=Fraction(-56,145),
                          plus_minus_alternating=Fraction(-1,145),
                          oriented_plus_minus=Fraction(0)), 'MRB mask means')
    payload = dict(universe="all even 4 <= m <= 10000; primes <= 4000000",
                   table=table, record_ratios=records,
                   max_first_prime=max(row['first_prime'] for row in rows),
                   crt_modulus=M, crt_obstructions=obstruction,
                   quotient_U290_by_sign_order_histogram=orders,
                   square_roots_of_plus_one=[a for a in units if a*a%290==1],
                   square_roots_of_minus_one=[a for a in units if a*a%290==289],
                   stern_227=stern,
                   mrb_mask_means={k:str(v) for k,v in mask_means.items()},checks="PASS")
    print(json.dumps(payload,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
