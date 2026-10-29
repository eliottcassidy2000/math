"""Exact modular response quotients; no exact-singularity observation.

All (z,a,b) modulo each declared modulus, all formal translations. A 5x5
grid determines the complete integer polynomial response function by finite
differences. At small moduli compare against every residue query as well.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def q(a,b):
    return a*a+a*b+b*b


def ternary_exponent(m):
    k,n = 0,m
    while n > 1 and n % 3 == 0:
        k,n = k+1,n//3
    return k if n == 1 else None


def census(m):
    grid = list(product(range(min(m,5)),repeat=2))
    full_grid = list(product(range(m),repeat=2)) if m <= 9 else None
    k = ternary_exponent(m)
    by_signature, by_key, full_by_signature = {}, {}, {}
    for a,b in product(range(m),repeat=2):
        values = [q(a+u,b+v)%m for u,v in grid]
        full_values = ([q(a+u,b+v)%m for u,v in full_grid]
                       if full_grid is not None else None)
        for z in range(m):
            signature = bytes((z-v)**2 % m for v in values)
            by_signature[signature] = by_signature.get(signature,0)+1
            if k is not None and k >= 2:
                key = (z,a%(m//3),b%(m//3),None if z%3 == 1 else (a-b)%m)
                need(by_key.setdefault(key,signature) == signature,
                     'predicted ternary equivalence implies identical responses')
            if full_values is not None:
                full = bytes((z-v)**2%m for v in full_values)
                need(full_by_signature.setdefault(signature,full) == full,
                     'grid equality gives full residue-query equality')
    if k is not None and k >= 2:
        need(len(by_key) == len(by_signature), 'ternary classes have no further collisions')
    histogram = dict(sorted(Counter(by_signature.values()).items()))
    need(sum(size*classes for size,classes in histogram.items()) == m**3,
         'complete source universe conserved')
    if k:
        need(len(by_signature) == 7*27**(k-1), 'all-order ternary recurrence')
        need(histogram == {3:6*27**(k-1),9:27**(k-1)}, 'ternary fiber histogram')
    if m in (5,7):
        need(histogram == {1:m**3}, 'prime-to-six carrier is faithful')
    return {'modulus':m,'source_states':m**3,'response_classes':len(by_signature),
            'fiber_histogram':histogram,'full_residue_queries_compared':full_grid is not None,
            'ternary_coordinate_classification_checked':bool(k and k>=2)}


def main():
    moduli = (1,2,3,4,5,7,8,9,12,15,16,20,27,76,81)
    rows = [census(m) for m in moduli]
    counts = {r['modulus']:r['response_classes'] for r in rows}
    for m,n in [(3,4),(3,5),(4,5)]:
        need(counts[m*n] == counts[m]*counts[n], 'compatible CRT product')
    need(counts[76] == 8*19**3 == 54872, 'named modulus76 exact census')
    need([counts[m] for m in (2,4,8,16)] == [8,8,40,320], 'dyadic finite census only')
    need((2-q(0,0))**2%3 == (2-q(1,0))**2%3 == 1,
         'z2 base-level exceptional collapse')
    # A claimed one-bit-for-all-primes decoder would miss these exact fibers.
    need(counts[3] != 3**3 and counts[9] != 9**3, 'ternary faithfulness hostile')
    out = {'gates':GATES,'rows':rows,'crt_proved_count_228':8*7*19**3,
           'scope':'complete modular determinant functions; no exact singularity flag',
           'dyadic_all_order_status':'PROVED in DYADIC_AND_FIELD_PORTS.md; this script retains the original finite census'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
