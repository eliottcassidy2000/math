#!/usr/bin/env python3
"""Exact controls for the clock/factor guard; no mass mutation search."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def tournament(rows):
    n = len(rows)
    return all(not (rows[i] >> i & 1) for i in range(n)) and all(
        ((rows[i] >> j & 1)+(rows[j] >> i & 1)) == 1
        for i in range(n) for j in range(i))

def hp(rows):
    n = len(rows)
    incoming = tuple(sum(1 << u for u in range(n) if rows[u] >> v & 1)
                     for v in range(n))
    @lru_cache(None)
    def end(mask, v):
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        possible = rest & incoming[v]
        return sum(end(rest, u) for u in range(n) if possible >> u & 1)
    return 1 if not n else sum(end((1 << n)-1, v) for v in range(n))

def literal_hp(rows):
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(rows)-1))
               for p in permutations(range(len(rows))))

def ordinal(left, right):
    a, b = len(left), len(right)
    upper = ((1 << b)-1) << a
    return tuple(row | upper for row in left) + tuple(row << a for row in right)

def strong(rows):
    n = len(rows)
    for start in range(n):
        reached = 1 << start
        while True:
            extended = reached
            for v in range(n):
                if reached >> v & 1:
                    extended |= rows[v]
            if extended == reached:
                break
            reached = extended
        if reached != (1 << n)-1:
            return False
    return True

def modular_step(m, epsilon):
    if m % 8 == (4+epsilon) % 8:
        return 3*m+2*epsilon
    require(m % 8 == epsilon % 8, 'moving-modulus domain')
    return (m+epsilon)//2

def main():
    cases = Counter()
    for epsilon in (-1, 1):
        for x in range(1, 512, 2):
            m = 4*x+epsilon
            y = 3*x+epsilon
            moving = modular_step(m, epsilon)
            require(moving == 4*y+epsilon, 'odd conjugacy')
            if m % 3 == 0:
                require(moving % 3 != 0, 'odd exit from triangle scalar guard')
            j = 0
            while y % 2 == 0:
                y //= 2
                moving = modular_step(moving, epsilon)
                j += 1
                numerator = 3*m+(2**j+1)*epsilon
                require(numerator % (2**j) == 0, 'clock formula integrality')
                require(moving == numerator//(2**j) == 4*y+epsilon, 'clock formula')
                require((moving % 3 == 0) == (j % 2 == 1), 'alternating triangle residue')
                if j == 1:
                    require(moving == 3*(m+epsilon)//2, 'one-halved formula')
                cases['legal_halving_endpoints'] += 1
            cases['odd_sources'] += 1
        for x in range(2, 513, 2):
            m = 4*x+epsilon
            if m % 3 == 0:
                require(modular_step(m, epsilon) % 3 != 0, 'even exit from guard')
                cases['even_guard_exits'] += 1
    plus = [13]
    minus = [43]
    for _ in range(2):
        plus.append(modular_step(plus[-1], 1))
        minus.append(modular_step(minus[-1], -1))
    require(plus == [13, 41, 21], 'plus forbidden endpoint')
    require(minus == [43, 127, 63], 'minus factor hostile')

    # Inherited literal witness, not newly found or a fresh all-class census.
    rows63 = (162, 252, 249, 17, 1, 152, 185, 24)
    require(tournament(rows63) and strong(rows63), 'strong inherited H63 witness')
    require(hp(rows63) == literal_hp(rows63) == 63, 'two independent H63 counters')
    triangle = (2, 4, 1)
    rows189 = ordinal(triangle, rows63)
    require(tournament(rows189) and hp(rows189) == 189, 'one exposed triangle')
    budget = []
    for r in range(2, 17):
        n = 7*3**r
        k = r-2
        require(n//3**k == 63, 'constructive budget residual')
        require(n//3**(k+1) == 21 and n//3**(k+2) == 7, 'forbidden residuals')
        budget.append({'r':r, 'H':n, 'max_ordinal_C3_factors':k})
    for h in range(1, 502, 2):
        for delta in (2, 4):
            require((3*h+delta) % 3 == delta % 3 != 0, 'nonzero ternary digit guard')
            cases['digit_guard_checks'] += 1

    sources = [
        '01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md',
        '01-canon/theorems/THM-344-n8-h63-complete-omega.md',
        '04-computation/h63_counterexample_audit_s8.py',
        '05-knowledge/results/triples_ports_sieve_20261112/ports/FOUR_RETURNS.md',
        '05-knowledge/results/central_modular_chanpak_20261110/central/PROOF.md',
        '05-knowledge/results/odd_tables_signed_carries_20261111/tournament/CONTENT.md',
    ]
    print(json.dumps({
        'status':'PASS; FINITE-EXACT controls with inherited all-order dependencies',
        'universe':{'odd_x_min':1,'odd_x_max':511,'signs':[-1,1],
                    'halving_depth':'all legal depths after each sampled odd source',
                    'digit_H_min':1,'digit_H_max':501,'factor_r':[2,16]},
        'checks':dict(cases),
        'clock_hostiles':{'plus':plus,'minus':minus},
        'literal_witnesses':[
            {'rows':rows63,'H':63,'n':8,'strong':True,'counters':['subset DP','all 8! permutations']},
            {'rows':rows189,'H':189,'n':11,'counters':['subset DP'],'ordinal_triangle':True}],
        'factor_budget_controls':budget,
        'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'limits':['No all-length mutation-invariance claim',
                  'No inference that divisibility implies an ordinal factor',
                  'No Collatz convergence or Chan-Pak completeness conclusion'],
    }, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
