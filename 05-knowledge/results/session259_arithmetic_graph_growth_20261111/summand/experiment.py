#!/usr/bin/env python3
"""Small exact probes of new Session259 summand extensions; no old large census."""
from math import gcd
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def closure(a, b, limit):
    live = {a, b}
    for n in range(1, limit + 1):
        if any(x in live and n-x in live for x in range(1, (n+1)//2)):
            live.add(n)
    return {x for x in live if x <= limit}


def fibre(live, n):
    return sum(x in live and n-x in live for x in range(1, (n+1)//2))


def main():
    limit = 500
    pairs = 0
    coprime_cases = 0
    examples = []
    for a in range(1, 16):
        for b in range(a+1, 17):
            pairs += 1
            live = closure(a, b, limit)
            expected = {a, b} | {i*a+j*b for i in range(1, limit//a+1)
                                            for j in range(1, limit//b+1)
                                            if i*a+j*b <= limit}
            require(live == expected, ('closure', a, b))
            d = gcd(a, b)
            A, B = a//d, b//d
            relative_holes = [n for n in range(1, limit//d+1) if d*n not in live]
            if A >= 2:
                require(max(relative_holes) == A*B, ('final scaled hole', a, b))
                require(len(relative_holes) == (A*B+A+B-5)//2, ('genus', a, b))
            else:
                require(relative_holes == list(range(2, B)), ('seed one', a, b))
            if d == 1 and a >= 2:
                coprime_cases += 1
                holes = set(relative_holes)
                h, q, F = len(holes), sum(holes), max(holes)
                total = 0
                for n in range(1, limit+1):
                    r = fibre(live, n)
                    total += r
                    loss = sum(x < n for x in holes) - sum(
                        n-x in holes for x in holes if x <= n-x)
                    require(r == (n-1)//2-loss, ('deletion', a, b, n))
                    if n >= 2*F:
                        predicted = (n-1)**2//4-h*n+q+h*(h+1)//2
                        require(total == predicted, ('cumulative', a, b, n))
                if (a, b) in [(2, 3), (3, 5), (4, 7)]:
                    examples.append({'seeds':[a,b], 'holes': sorted(holes),
                                     'genus':h, 'last_hole':F, 'hole_sum':q})
    fronts = []
    for b in [3, 5, 7, 9, 11, 15, 21]:
        live = {2,b}
        T = (b+5)//2
        row = []
        for t in range(9):
            M = max(live)
            if t >= 3:
                MT = (3*b+4)*2**(T-3)+2
                predicted = ((3*b+4)*2**(t-3)+2 if t <= T
                             else (MT-1)*2**(t-T)+1)
                require(M == predicted, ('frontier', b, t))
                if t <= T:
                    opposite = max(x for x in live if (M-x)%2)
                    require(M-opposite == b-2*(t-3), ('parity deficit', b, t))
                    require({M,M-2,M-4} <= live, ('same-parity top', b,t))
                if t >= T:
                    require({M,M-1,M-2} <= live, ('consecutive top', b,t))
            row.append(M)
            if t < 8:
                ordered = sorted(live)
                new = live | {x+y for i,x in enumerate(ordered) for y in ordered[i+1:]}
                D = {M-x for x in live}
                gap = min(D-{0})
                predicted_D = {x+y-gap for x in D for y in D if x<y} | {M-gap+x for x in D}
                require(max(new) == 2*M-gap, ('deficit maximum',b,t))
                require(predicted_D == {max(new)-x for x in new}, ('deficit law',b,t))
                live = new
        fronts.append({'b':b, 'T_proved':T, 'maxima_t0_to_t8':row})
    print(json.dumps({'status':'FINITE-EXACT PASS', 'numeric_limit':limit,
                      'all_seed_pairs_checked':pairs,
                      'coprime_deletion_cases':coprime_cases,
                      'examples':examples, 'frontiers':fronts,
                      'symbolic_b101': {'T':53,'M_T':307*2**50+2,
                          'status':'formula evaluation, not enumerated state'}}, indent=2))

if __name__ == '__main__':
    main()
