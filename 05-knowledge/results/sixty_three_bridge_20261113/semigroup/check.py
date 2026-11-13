"""Small exact controls for activation ideals and deleted-atom countermodels.

No tournament census: the physical source counts 3,5,35,49,63,343 are
inherited from THM-4094 and its finite-prefix input. Assertions use need(),
so optimized Python retains every check.
"""
from itertools import product
from math import isqrt, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def in_m(n):
    return n > 0 and n % 2 == 1 and n not in (7, 21)


def atom_formula(n):
    return (n % 2 == 1 and prime(n) and n != 7
            or n % 7 == 0 and n // 7 >= 5 and n // 7 != 7 and prime(n // 7)
            or n in (49, 63, 343))


def atom_divisors(n):
    return n > 1 and in_m(n) and not any(
        n % d == 0 and in_m(d) and in_m(n // d)
        for d in range(3, isqrt(n) + 1, 2))


def deleted(n, cutoff=13):
    return (prime(n) and n % 2 == 1 and n > cutoff and n != 7
            or n % 7 == 0 and n // 7 > cutoff and n // 7 != 7 and prime(n // 7))


def in_s(n, cutoff=13):
    return in_m(n) and not deleted(n, cutoff)


def in_conductor(n, cutoff=13):
    return (n > 0 and n % 2 == 1 and n not in (1, 3, 7, 21)
            and not deleted(n, cutoff)
            and not (prime(n) and deleted(7 * n, cutoff)))


def tau_formula(n, cutoff=13):
    if n == 7:
        return 2
    return int(deleted(n, cutoff))


def source_recipe(a, b, c):
    if b:
        return [35] + [3] * a + [5] * (b - 1) + [49] * c
    if c:
        return [343] + [3] * a + [49] * (c - 1)
    if a >= 2:
        return [63] + [3] * (a - 2)
    return None


def main():
    counts = dict(atom_checks=0, closure_pairs=0, conductor_products=0,
                  square_cube_inputs=0, localization_numerators=0,
                  activation_inputs=0, subadditivity_pairs=0,
                  activation_addresses=0)
    for n in range(1, 1002, 2):
        need(atom_formula(n) == atom_divisors(n), 'inherited arithmetic atom list')
        counts['atom_checks'] += 1

    # An explicit infinite deleted-atom family: all ordinary p and seven-p
    # with p>13. The bounded tests below use the exact untruncated predicate.
    inputs = range(1, 202, 2)
    for a, b in product(inputs, repeat=2):
        if in_s(a) and in_s(b):
            need(in_s(a * b), 'deleted-atom monoid closure')
            counts['closure_pairs'] += 1
        if in_conductor(a):
            need(in_s(a * b), 'conductor positive direction')
            counts['conductor_products'] += 1
    for a in inputs:
        if not in_conductor(a):
            if a in (1, 3, 7, 21):
                witness = {1: 7, 3: 7, 7: 1, 21: 1}[a]
            elif deleted(a):
                witness = 1
            else:
                witness = 7
            need(not in_s(a * witness), 'explicit conductor exclusion witness')

    for n in range(1, 1002, 2):
        need(in_s(n * n) and in_s(n ** 3), 'all odd squares and cubes survive')
        need(in_s(9 * n) and in_s(49 * n), 'universal conductor multipliers')
        need(in_s(3 * n) == (n != 7), 'complete three-divisible slice')
        counts['square_cube_inputs'] += 1
    for a in range(3, 52, 2):
        if not in_s(a):
            continue
        for n in range(1, 202, 2):
            need(in_s(n * a * a), 'uniform localization numerator')
            counts['localization_numerators'] += 1

    free = [n for n in range(1, 202, 2) if n % 3]
    for n in free:
        brute_tau = next(k for k in range(3) if in_s(n * 3 ** k))
        need(brute_tau == tau_formula(n), 'exact activation threshold')
        counts['activation_inputs'] += 1
    for a, b in product(free, repeat=2):
        need(tau_formula(a * b) <= tau_formula(a) + tau_formula(b),
             'activation subadditivity')
        counts['subadditivity_pairs'] += 1

    examples = []
    for a, b, c in product(range(6), repeat=3):
        value = 7 * 3 ** a * 5 ** b * 49 ** c
        recipe = source_recipe(a, b, c)
        active = a >= 2 or b >= 1 or c >= 1
        need(active == in_m(value), 'exact actual smooth activation ideal')
        need(active == (recipe is not None), 'source recipe exists exactly when active')
        if active:
            need(prod(recipe) == value, 'ordinal source recipe exact target')
            need(all(q in (3, 5, 35, 49, 63, 343) for q in recipe), 'known source bank')
        else:
            need(value in (7, 21), 'finite activation box contains only known holes')
        if (a, b, c) in [(0, 0, 0), (1, 0, 0), (2, 0, 0),
                          (0, 1, 0), (0, 0, 1), (2, 2, 1)]:
            examples.append(dict(address=[a, b, c], value=value,
                                 active=active, ordinal_factors=recipe))
        counts['activation_addresses'] += 1

    # Finite prefix preservation is tested for several cutoff choices;
    # the theorem permits any finite prefix and infinitely many later holes.
    for cutoff in (7, 13, 31, 101):
        need(all(in_s(n, cutoff) == in_m(n) for n in range(1, cutoff + 1, 2)),
             'arbitrary prescribed prefix preserved')

    result = dict(
        status='PASS', counts=counts,
        scope='Arithmetic controls only; no tournament census and no claim that S_E is an actual tournament spectrum.',
        infinite_countermodel='E={p,7p:p prime>13}; S_E=positive odds minus {7,21} minus E',
        named_missing=[n for n in range(1, 130, 2) if not in_s(n)],
        activation_ideal=dict(base=7, multipliers=[3, 5, 49],
                              minimal_addresses=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                              source_values=[63, 35, 343], examples=examples),
        controls=dict(known_holes=[7, 21],
                      first_ternary_activation=[7, 21, 63],
                      removed_ordinary_atom=dict(value=17, first_activation=51),
                      removed_seven_prime_atom=dict(value=119, first_activation=357),
                      conductor_multipliers=[9, 49]),
        inherited_source='THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md; Session193 prime_induction/PROOF.md',
    )
    (HERE / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status='PASS', counts=counts, examples=examples), indent=2))


if __name__ == '__main__':
    main()
