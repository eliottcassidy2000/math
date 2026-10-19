#!/usr/bin/env python3
"""Independent bounded controls, not a proof of prime infinitude.
Uses factor-product sigma and an Eratosthenes sieve, without importing the verifier.
"""
from fractions import Fraction
import json

LIMIT = 20000

def check(ok, detail):
    if not ok:
        raise RuntimeError(detail)

def sigma(n):
    check(isinstance(n, int) and n >= 1, ('domain', n))
    remaining, result, divisor = n, 1, 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            power, geometric = 1, 1
            while remaining % divisor == 0:
                remaining //= divisor
                power *= divisor
                geometric += power
            result *= geometric
        divisor = 3 if divisor == 2 else divisor + 2
    if remaining > 1:
        result *= remaining + 1
    return result

def aliquot(n):
    return sigma(n) - n

flags = bytearray(b'\x01') * (LIMIT + 5)
flags[0:2] = b'\x00\x00'
for p in range(2, len(flags)):
    if flags[p] and p * p < len(flags):
        flags[p*p::p] = b'\x00' * (((len(flags)-1-p*p)//p)+1)
primes = [p for p in range(2, LIMIT+1) if flags[p]]

def prime(n):
    return n >= 2 and n < len(flags) and bool(flags[n])

# Cross-check this factor-product engine against direct proper-divisor sums.
for n in range(1, 1001):
    check(aliquot(n) == sum(d for d in range(1, n) if n % d == 0), ('sigma', n))

# Two-prime/common-shift uniqueness, with exact rational shifts, all positive cores <=1000.
perfect_small = []
for core in range(1, 1001):
    p, q = next((p,q) for p,q in [(1009,1013)] if core % p and core % q)
    hp = Fraction(aliquot(core*p), core) - p
    hq = Fraction(aliquot(core*q), core) - q
    is_perfect = sigma(core) == 2*core
    check((hp == hq) == is_perfect, ('two-prime rigidity', core, hp, hq))
    if is_perfect:
        check(hp == 2, ('unique shift', core, hp))
        perfect_small.append(core)
check(perfect_small == [6,28,496], ('small perfect list', perfect_small))

# Each perfect-core pair is tested with both endpoint guards retained.
perfect_counts = {}
for core in (6,28,496,8128):
    count = 0
    for p in primes:
        if p < 5 or core % p == 0:
            continue
        value = aliquot(core*p)
        check(value == core*(p+2), ('translation', core,p,value))
        if prime(p+2) and core % (p+2):
            count += 1
            check(value//core == p+2, ('target',core,p))
    perfect_counts[str(core)] = count

# Count actual edges from independently evaluated aliquot values, and direct orbit exits.
actual, twins = [], []
exits = {1:0, 2:0}
for p in primes:
    if p < 5:
        continue
    value = aliquot(6*p)
    q, residue = divmod(value, 6)
    if residue == 0 and q >= 5 and prime(q):
        actual.append((p,q))
    if prime(p+2):
        twins.append((p,p+2))
    current = 6*p
    for step in (1,2):
        current = aliquot(current)
        marked, residue = divmod(current,6)
        if residue or marked < 5 or not prime(marked):
            exits[step] += 1
            break
    else:
        raise RuntimeError(('lane exceeds two steps',p))
check(actual == twins, ('bijection', actual, twins))
for cutoff in (0,29,30,41,42,100,1000,6*LIMIT):
    check(sum(6*p <= cutoff for p,q in actual) ==
          sum(p <= cutoff//6 for p,q in twins), ('finite cutoff',cutoff))

# Native-domain and finite-exception controls.
check(aliquot(18) == 21 != 30, 'source 3 divides core 6')
check(aliquot(54) == 66 == 6*(9+2) and not prime(9), 'accidental composite agreement')
check(aliquot(150) == 222 != 162, 'composite disagreement')
check(aliquot(28*5) == 28*7 and 28 % 7 == 0, 'target guard finite exception')
check(aliquot(30) == 42 and aliquot(42) == 54, 'two-step exit bound attained')
check(aliquot(66) == 78 and aliquot(78) == 90, 'second two-step control')

print(json.dumps({'status':'PASS','prime_source_limit':LIMIT,
                  'positive_cores_tested_for_rigidity':1000,
                  'perfect_cores_up_to_1000':perfect_small,
                  'guarded_twin_edges_by_perfect_core':perfect_counts,
                  'fixed_6_twin_edges':len(actual),
                  'fixed_6_exit_steps':exits,
                  'native_domain_controls':6}, indent=2, sort_keys=True))
