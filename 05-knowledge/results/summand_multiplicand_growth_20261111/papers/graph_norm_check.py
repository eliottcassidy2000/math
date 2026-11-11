"""Literal finite controls for the cyclic-set norm and pair-growth chart."""
from collections import Counter
from itertools import product
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def perm(lengths):
    out = []
    for d in lengths:
        base = len(out)
        out.extend(base + (i + 1) % d for i in range(d))
    return tuple(out)


def cycles(P):
    seen = set()
    out = Counter()
    for i in range(len(P)):
        if i in seen:
            continue
        j, n = i, 0
        while j not in seen:
            seen.add(j)
            n += 1
            j = P[j]
        out[n] += 1
    return out


def add(P, Q):
    return tuple(P) + tuple(len(P) + j for j in Q)


def mul(P, Q):
    return tuple(P[i] * len(Q) + Q[j] for i in range(len(P)) for j in range(len(Q)))


def power(P, r):
    out = (0,)
    for _ in range(r):
        out = mul(out, P)
    return out


def frobenius(P, r):
    out = []
    for i in range(len(P)):
        j = i
        for _ in range(r):
            j = P[j]
        out.append(j)
    return tuple(out)


def verschiebung(P, r):
    return tuple(i*r+j+1 if j+1 < r else P[i]*r
                 for i in range(len(P)) for j in range(r))


def norm(P, r):
    words = list(product(range(len(P)), repeat=r))
    index = {word: i for i, word in enumerate(words)}
    return tuple(index[(P[word[-1]],) + word[:-1]] for word in words)


def ghost(P, m):
    return sum(j == i for i, j in enumerate(frobenius(P, m)))


def copies(P, n):
    out = ()
    for _ in range(n):
        out = add(out, P)
    return out


def main():
    family = [perm(ds) for ds in ((), (1,), (2,), (3,), (4,), (1, 1), (1, 2), (2, 2))]
    cycle_product_controls = 0
    for a in range(1, 7):
        for b in range(1, 7):
            assert cycles(mul(perm((a,)), perm((b,)))) == {lcm(a, b): gcd(a, b)}
            cycle_product_controls += 1

    operator_controls = 0
    ghost_controls = 0
    for P in family:
        for r in (1, 2, 3, 4):
            N = norm(P, r)
            assert cycles(frobenius(verschiebung(P, r), r)) == cycles(copies(P, r))
            assert cycles(frobenius(N, r)) == cycles(power(P, r))
            operator_controls += 2
            for m in range(1, 19):
                d = gcd(r, m)
                assert ghost(N, m) == ghost(P, m//d)**d
                ghost_controls += 1
            for Q in family[:4]:
                assert cycles(verschiebung(mul(P, frobenius(Q, r)), r)) == cycles(mul(verschiebung(P, r), Q))
                assert cycles(norm(mul(P, Q), r)) == cycles(mul(N, norm(Q, r)))
                operator_controls += 2

    carry_controls = 0
    for P in family[1:]:
        for Q in family[1:4]:
            for p in (2, 3):
                rhs = add(norm(P, p), norm(Q, p))
                for a in range(1, p):
                    assert comb(p, a) % p == 0
                    piece = verschiebung(mul(power(P, a), power(Q, p-a)), p)
                    rhs = add(rhs, copies(piece, comb(p, a)//p))
                assert cycles(norm(add(P, Q), p)) == cycles(rhs)
                carry_controls += 1

    birth_controls = 0
    for n in range(13):
        P = perm((1,)*n)
        assert cycles(norm(P, 2)) == Counter({1: n, 2: n*(n-1)//2})
        rhs = add(add(norm(P, 2), (0,)), verschiebung(P, 2))
        assert cycles(norm(add(P, (0,)), 2)) == cycles(rhs)
        birth_controls += 1

    # The two formulas distinguish actual graph identities from false upgrades.
    vf = cycles(verschiebung(frobenius((0,), 2), 2))
    twice = cycles(copies((0,), 2))
    assert vf != twice
    assert cycles(norm(perm((2,)), 2)) == {4: 1}
    assert ghost((), 1) == ghost(perm((2,)), 1) == 0
    assert ghost(frobenius(perm((2,)), 2), 1) == 2
    assert cycles(norm(perm((1, 1)), 6)) == {1: 2, 2: 1, 3: 2, 6: 9}

    chart = {}
    chart_controls = 0
    for n in range(1, 41):
        for a in range(1, n+1):
            s, t = a+n, a*n
            delta = s*s-4*t
            q = isqrt(delta)
            assert q*q == delta and (s-q)//2 == a and (s+q)//2 == n
            assert t == n*s-n*n
            assert (delta == 0) == (a == n)
            assert (s, t) not in chart
            chart[s, t] = (a, n)
            chart_controls += 1
    stabilization_controls = 0
    for M in range(2, 31):
        sums = {(a,b) for a in range(1,2*M+1) for b in range(a,2*M+1) if a+b <= M}
        products = {(a,b) for a in range(1,2*M+1) for b in range(a,2*M+1) if a*b <= M}
        assert max(b for a,b in sums) == M-1
        assert max(b for a,b in products) == M
        stabilization_controls += 2

    receipt = {
        'status': 'FINITE-EXACT controls for all-order proofs in PROOF.md',
        'finite_universe': {'cycle_factors_a_b': [1,6], 'permutation_cycle_types': [dict(cycles(P)) for P in family],
                            'r_values': [1,2,3,4], 'ghost_m_values': [1,18], 'carry_primes': [2,3],
                            'integer_birth_n': [0,12], 'joint_chart_max_label': 40, 'output_bounds': [2,30]},
        'controls': {'cycle_products': cycle_product_controls, 'operator_identities': operator_controls,
                     'ghost_norm': ghost_controls, 'positive_norm_carries': carry_controls,
                     'one_integer_births': birth_controls, 'joint_chart': chart_controls,
                     'stabilization': stabilization_controls},
        'hostiles': {'VF_on_one_point': dict(vf), 'two_fixed_points': dict(twice),
                     'norm2_of_C2': {4: 1}, 'hidden_C2_initial_g1': 0, 'hidden_C2_after_F2_g1': 2,
                     'nonnegative_ghosts_2_0_require_c2': -1},
        'composite_norm_example': {'N6_two_fixed_points': {1: 2, 2: 1, 3: 2, 6: 9}},
        'joint_birth_example_n5': [{'pair': [a,5], 'sum': a+5, 'product': a*5, 'fixed': a==5} for a in range(1,6)],
    }
    (HERE/'graph_norm_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'pass', **receipt['controls']}))


if __name__ == '__main__':
    main()
