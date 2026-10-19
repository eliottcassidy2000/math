#!/usr/bin/env python3
"""Exact controls for T2_GLOBAL_CLASSIFICATION.md; not a substitute for its proof."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

E = {2, 4, 6, 8, 10, 12, 16, 24, 28, 32, 40, 44, 48}
CYCLES = [
    [16, 78, 46, 30, 22, 18],
    [24, 52, 54, 34],
    [28, 36, 42],
    [32, 158, 86, 50],
    [40, 66],
    [48, 108, 96, 220, 180, 150, 82],
]
CHAINS = [
    [20, 30, 22, 18],
    [56, 80, 138, 76, 72, 94, 54],
    [88, 108, 96, 220, 180, 150, 82],
    [184, 192, 444, 348, 276, 222, 118],
    [192, 444, 348, 276, 222, 118],
    [224, 344, 332, 264, 262, 138],
    [384, 892, 684, 528, 558, 286],
    [448, 696, 640, 1146, 580, 450, 232],
    [752, 768, 1788, 1356, 1032, 934, 474],
    [768, 1788, 1356, 1032, 934, 474],
    [1536, 3580, 2700, 2040, 1816, 1620, 1230],
    [1792, 2808, 2488, 2208, 2266, 1140],
]
SMALL_STEPS = {row[0]: len(row) - 1 for row in CHAINS}


def core(n: int) -> int:
    assert n > 0
    return n & -n


def g(n: int) -> int:
    assert n > 0 and n % 2 == 0
    c = core(n)
    p = n // c
    return (c - 1) * p + 4 * c - 1


def orbit(n: int, steps: int = 6) -> list[int]:
    result = [n]
    for _ in range(steps):
        result.append(g(result[-1]))
    return result


def check_chain(row: list[int]) -> None:
    assert all(g(a) == b for a, b in zip(row, row[1:])), row


def symbolic_controls() -> int:
    c, d, e, h, p, t, k = sp.symbols('c d e h p t k')
    identities = []

    def eq(lhs, rhs):
        assert sp.factor(lhs - rhs) == 0, (lhs, rhs)
        identities.append(str(rhs))

    def step(n, C, parameter=2):
        return n - n / C + 2 * parameter * C - 1

    n = c * p
    n1 = step(n, c, t)
    n2 = step(n1, d, t)
    eq(d * (n2 - n), (c + d) * (2*t*(d - 1) - p) + 2*(t - 1)*d + p + 1)

    n0 = e*d*(d*k - 1)
    n1 = step(n0, e*d)
    eq(n1, d*(e*(d*k + 3) - k))
    n2 = step(n1, d)
    eq(n2, e*(d - 1)*(d*k + 3) + d*(4 - k) + k - 1)

    # k=3, d<c.
    n0_3 = n0.subs(k, 3)
    n2_3 = n2.subs(k, 3)
    eq(n2_3, 3*e*(d*d - 1) + d + 2)
    eq(step(n2_3, 2) - n0_3, -e*(3*d*d - 2*d + 3)/2 + d/2 + 8)
    eq(step(n2_3.subs(e, 2), 4) - n0_3.subs(e, 2), -3*d*d/2 + 11*d/4 + 12)

    # k=3 diagonal, with proved successive cores.
    z = c*(3*c - 1)
    closed = [3*c*c, 7*c*c-4, 21*c*c/4+12, 63*c*c/16+24,
              441*c*c/128+52, 1323*c*c/512+54]
    for C, expected in zip([c, c*c, 4, 4, 8, 4], closed):
        z = step(z, C)
        eq(z, expected)
    eq(z - c*(3*c - 1), -213*c*c/512+c+54)

    # k=1, d=2c and d=c.
    z = 2*c*c-c
    for C, expected in zip([c, c, 2], [c*(2*c+1), 2*c*c+3*c-2, c*c+3*c/2+6]):
        z = step(z, C)
        eq(z, expected)
    eq(z-(2*c*c-c), -c*c+5*c/2+6)
    z = c*(c-1)
    for C, expected in zip([c, 2*c, 2], [c*(c+2), c*c+19*c/2-2, c*c/2+19*c/4+6]):
        z = step(z, C)
        eq(z, expected)
    eq(z-c*(c-1), -c*c/2+23*c/4+6)

    # k=1, e<d.
    n0_1 = n0.subs(k, 1)
    n2_1 = n2.subs(k, 1)
    eq(n2_1, e*(d*d+2*d-3)+3*d)
    n3_1 = step(n2_1, e)
    eq(n3_1, e*(d+1)**2-d*d+d+2-3*d/e)
    eq(n3_1-n0_1, e*(3*d+1)-d*d+d+2-3*d/e)
    eq(n3_1.subs(e, d/2), d**3/2+3*d/2-4)
    eq(step(n3_1.subs(e, d/2), 4)-n0_1.subs(e, d/2),
       -d**3/8+d*d/2+9*d/8+12)

    # k=1, e=d.
    eq(n2_1.subs(e, d), d*d*(d+2))
    z = step(n2_1.subs(e, d), 2*d*d)
    eq(z, d**3+10*d*d-d/2-2)
    eq(step(z, 2)-n0_1.subs(e, d), -d**3/2+6*d*d-d/4+6)

    # k=1, e>d, e=hd.
    z = step(n2_1.subs(e, h*d), d)
    eq(z, h*(d-1)**2*(d+3)+7*d-4)
    baseline = h*d*d*(d-1)
    eq((step(z, 2)-baseline).subs(h, 2), -d**3+3*d*d-3*d/2+8)
    eq(step(z, 4)-baseline, -h*(d-1)*(d-3)**2/4+21*d/4+12)
    eq((step(z, 8)-baseline).subs(h, 4), -d**3/2+15*d*d/2-91*d/8+38)
    z4 = step(z.subs(d, 4), 8)
    eq(z4, 441*h/8+52)
    eq(step(z4, 4)-baseline.subs(d, 4), -213*h/32+54)
    eq(step(z.subs(d, 2), 2)-baseline.subs(d, 2), 12-3*h/2)
    return len(identities)


def prediction(n: int) -> tuple[str, int | None, list[tuple[int, int]]]:
    """Proved case partition. A guard pair (i,C) means core(G^i(n))=C."""
    c = core(n)
    q = n//c + 1
    if q > 4*c:
        return 'immediate', 1, []
    if q == 4*c:
        return 'fixed', None, []
    d = core(q)
    k = q//d
    if k >= 5:
        return 'general_oddpart', 2, [(1, d)]
    if n in E:
        return 'exception', None, []
    if n in SMALL_STEPS:
        return 'finite_chain', SMALL_STEPS[n], []
    if k == 3:
        if d == c:
            assert c >= 32
            return 'k3_diagonal', 6, [(1,c*c),(2,4),(3,4),(4,8),(5,4)]
        assert d < c
        e = c//d
        if d >= 4 and e >= 4:
            return 'k3_e_large', 3, [(1,d),(2,2)]
        if e == 2 and d >= 8:
            return 'k3_e2', 3, [(1,d),(2,4)]
        assert d == 2 and e >= 8
        return 'k3_d2', 2, [(1,d)]
    assert k == 1
    if d == 2*c:
        assert c >= 8
        return 'k1_d2c', 3, [(1,c),(2,2)]
    if d == c:
        assert c >= 16
        return 'k1_dc', 3, [(1,2*c),(2,2)]
    assert d < c
    e = c//d
    if e < d:
        if e <= d//4:
            return 'k1_e_small', 3, [(1,d),(2,e)]
        assert e == d//2 and d >= 16
        return 'k1_e_half', 4, [(1,d),(2,e),(3,4)]
    if e == d:
        assert d >= 16
        return 'k1_e_diagonal', 4, [(1,d),(2,2*d*d),(3,2)]
    h = e//d
    assert e == h*d and h >= 2
    if h == 2 and d >= 4:
        return 'k1_h2', 4, [(1,d),(2,d),(3,2)]
    if h >= 8 and d >= 8:
        return 'k1_h_large', 4, [(1,d),(2,d),(3,4)]
    if h == 4 and d >= 16:
        return 'k1_h4', 4, [(1,d),(2,d),(3,8)]
    if d == 4:
        assert h >= 64
        return 'k1_d4', 5, [(1,d),(2,d),(3,8),(4,4)]
    assert d == 2 and h >= 16
    return 'k1_d2', 4, [(1,d),(2,d),(3,2)]


def check_n(n: int, counts: Counter) -> None:
    label, j, guards = prediction(n)
    counts[label] += 1
    if j is None:
        assert (g(n) == n) if label == 'fixed' else (n in E)
        return
    path = orbit(n, j)
    for i, C in guards:
        assert core(path[i]) == C, (n,label,i,C,core(path[i]))
    assert path[-1] < n, (n,label,j,path)
    assert 1 <= j <= 6


def guard_boundary_controls() -> list[dict]:
    # Each invalid lowered threshold would assert a false exact core.
    controls = [
        ('k1_e_half requires d>=16', 224, 3, 4),
        ('k1_e_diagonal small d=4', 48, 3, 2),
        ('k1_h4 requires d>=16', 1792, 3, 8),
        ('k3_diagonal n5 requires c>=32', 752, 5, 4),
        ('k3_diagonal n3 requires c>=8', 44, 3, 4),
        ('k1_dc n2 requires c>=8', 12, 2, 2),
    ]
    result = []
    for description,n,i,false_core in controls:
        observed = core(orbit(n,i)[-1])
        assert observed != false_core
        result.append({'control':description,'start':n,'step':i,
                       'invalid_claimed_core':false_core,'actual_core':observed})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=1000000)
    parser.add_argument('--max-exponent', type=int, default=128)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    identities = symbolic_controls()
    for row in CHAINS:
        check_chain(row)
        assert row[-1] < row[0] and len(row)-1 <= 6
    cycle_members = set()
    for cycle in CYCLES:
        check_chain(cycle+[cycle[0]])
        assert not (cycle_members & set(cycle))
        cycle_members.update(cycle)
    for start in E:
        z = start
        seen = set()
        while z not in cycle_members:
            assert z not in seen
            seen.add(z)
            z = g(z)
    grid_counts = Counter()
    for a in range(1,args.max_exponent+1):
        c = 1 << a
        for b in range(1,args.max_exponent+1):
            d = 1 << b
            for k in [1,3,5,7,9,11]:
                check_n(c*(d*k-1), grid_counts)
    interval_counts = Counter()
    for n in range(2,args.max_n+1,2):
        check_n(n, interval_counts)
    result = {
        'status':'PASS',
        'scope':'Finite controls supporting the written all-exponent proof, not a census proof.',
        'symbolic_identities_checked':identities,
        'explicit_descent_chains_checked':len(CHAINS),
        'nontrivial_cycles_checked':CYCLES,
        'exceptional_starts_checked':sorted(E),
        'exponent_grid':{'a_max':args.max_exponent,'b_max':args.max_exponent,
                         'odd_coefficients':[1,3,5,7,9,11],
                         'inputs_checked':sum(grid_counts.values()),
                         'branch_counts':dict(sorted(grid_counts.items()))},
        'initial_interval':{'max_n':args.max_n,'even_inputs_checked':sum(interval_counts.values()),
                            'branch_counts':dict(sorted(interval_counts.items()))},
        'guard_boundary_controls':guard_boundary_controls(),
    }
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
