"""Bounded-prime, all-initial-state control for the squarefree selector.

Higher primes never affect the tested coordinates. Monotonicity makes the
empty and full prime banks exact lower/upper controls for every initial bank.
No claim of global finite-time stabilization is made.
"""
from hashlib import sha256
from pathlib import Path
import argparse
import json

checks = 0
def require(test, message):
    global checks
    checks += 1
    if not test:
        raise ValueError(message)

def control(limit):
    spf = list(range(limit + 1))
    for p in range(2, int(limit ** .5) + 1):
        if spf[p] == p:
            for n in range(p * p, limit + 1, p):
                if spf[n] == n:
                    spf[n] = p
    primes = [p for p in range(2, limit + 1) if spf[p] == p]
    dependencies = {}
    for p in primes:
        n = p - 1
        factors = []
        squarefree = True
        while n > 1:
            q = spf[n]
            n //= q
            factors.append(q)
            if n % q == 0:
                squarefree = False
                break
        if squarefree:
            dependencies[p] = frozenset(factors)
    core = {2, 3, 7, 43}.intersection(primes)
    update = lambda bank: {p for p, factors in dependencies.items() if factors <= bank}
    require(update(core) == core, 'fixed core in the tested window')
    lower, upper = set(), set(primes)
    rows = []
    for clock in range(limit.bit_length() + 2):
        require(lower <= core <= upper, 'monotone sandwich covers every initial bank')
        bad = upper - lower
        if clock >= 1:
            require(not bad or min(bad) >= (1 << (clock + 1)) - 1,
                    'uniform least discrepancy doubles after the first clock')
        rows.append(dict(clock=clock, lower_count=len(lower), upper_count=len(upper),
                         least_possible_discrepancy=min(bad) if bad else None))
        if not bad:
            break
        lower, upper = update(lower), update(upper)
    require(lower == upper == core, 'every initial bank stabilizes in this window')
    # Independent complete divisors for the actual finite fixed integer.
    divisors = [d for d in range(1, 1807) if 1806 % d == 0]
    def isprime(n):
        return n >= 2 and all(n % d for d in range(2, int(n ** .5) + 1))
    selected = [d + 1 for d in divisors if isprime(d + 1)]
    require(selected == [2, 3, 7, 43], 'complete divisor selector, not only n+1')
    return dict(limit=limit, prime_count=len(primes), rows=rows,
                fixed_integer=1806, selected_primes=selected,
                status='FINITE-EXACT bounded coordinates, all initial banks',
                explicit_checks=checks,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=1000000)
    args = parser.parse_args()
    if args.limit < 43:
        raise ValueError('use limit >= 43')
    print(json.dumps(control(args.limit), indent=2, sort_keys=True))
