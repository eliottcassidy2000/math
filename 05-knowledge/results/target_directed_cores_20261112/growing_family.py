"""Exact membership in the explicitly entered H35 growing-core family.

The parameter selector is complete for this family. This is not a complete
Hamilton spectrum selector, or an enumeration of all entered modules.
"""
from math import gcd, isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
D = (8, 5, 9, 18, 71, 31, 47)


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def state(a, b):
    need(a >= 0 and b >= 0 and a + b > 0, 'nonempty module')
    h, f = 35**b, 676**b * 2**a
    H, F = 1386*f - 13*h, 2592*f - 129*h
    u, B = 3978*f - 157*h, 14346*f - 703*h
    return dict(a=a, b=b, target=H, module_pair=[h, f],
                module_recipe=dict(ordinal_copies_of=list(D), copies=b,
                                   transitive_suffix=a),
                source_order=7*b+a+11, source_at_x=[H, 15*h, 15*h, F],
                RNR_order=7*b+a+14, balanced_state=[B, u, u, B],
                controller='C3 then inherited mod11 bridge / C5')


def select(N):
    """All (a,b) in this family, with their full physical source recipes."""
    answers = []
    b, hb, fb = 0, 1, 1
    while 1373*fb <= N:
        value, rem = divmod(N+13*hb, 1386*fb)
        if rem == 0 and value > 0 and value & (value-1) == 0:
            a = value.bit_length()-1
            if a+b > 0:
                answers.append(state(a, b))
        b += 1
        hb *= 35
        fb *= 676
    return answers


def main():
    p = 5039
    need(all(p % d for d in range(2, isqrt(p)+1)), 'prime5039')
    modular = {str(n): pow(n, 2519, p) for n in (2, 35, 676, 2365)}
    need(modular == {'2': 1, '35': 1, '676': 1, '2365': 5038}, 'modular gates')
    need(3978*2365 % p == 157, 'exceptional ratio')
    # A declared arithmetic family head, with no exhaustive graph claim.
    head = [state(a, b) for a in range(25) for b in range(25) if a+b]
    for row in head:
        h, f = row['module_pair']
        B, u, _, _ = row['balanced_state']
        need(gcd(h, f) == 1 and gcd(B, u) == 1, 'primitive family')
        need((B-2*u) % 3 == h % 3 != 0, 'directC3guard')
        recovered = select(row['target'])
        need((row['a'], row['b']) in {(r['a'], r['b']) for r in recovered}, 'exact inverse')
    # This fixed N bound admits only b=0,1,2 and at most24 values of a.
    cap = 1_000_000_000
    explicit = {}
    for b in range(3):
        for a in range(25):
            if a+b == 0:
                continue
            row = state(a, b)
            if row['target'] <= cap:
                explicit.setdefault(row['target'], []).append([a, b])
    need(1373*676**3 > cap and 1386*2**25-13 > cap, 'complete head cutoffs')
    for N, pairs in explicit.items():
        need(sorted([r['a'], r['b']] for r in select(N)) == sorted(pairs), 'complete parameter head')
    hostile = [7, 21, 881, 1691, 53084231, 159252695]
    need(all(not select(N) for N in hostile), 'scoped negative controls')
    # At b=0,a=1 the chart source must be the inherited middle462 value.
    need(state(1, 0)['target'] == 2759, 'positive base')
    receipt = dict(status='PASS', universe='624 parameter pairs; all family targets <=10^9',
                   modular_certificate=modular, checked_pairs=len(head),
                   complete_head_bound=cap, complete_head_targets=len(explicit),
                   complete_head=explicit, scoped_negative_targets=hostile,
                   examples=[state(a, b) for a, b in ((1, 0), (0, 1), (2, 1), (0, 2))],
                   scope='Family membership only; all-module and all-parameter proofs in GROWING_FAMILY.md')
    (HERE/'growing_family_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('PASS', len(head), 'parameter pairs;', len(explicit), 'targets in complete head')


if __name__ == '__main__':
    main()
