"""Independent arithmetic checks for the row-multiset information obstruction."""
import argparse
import decimal
import fractions
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def need(test, message):
    if not test:
        raise ValueError(message)


def composition_count(mass, bins):
    # Dynamic programming, independent of the closed stars-and-bars formula.
    counts = [1] + [0] * mass
    for _ in range(bins):
        for i in range(1, mass + 1):
            counts[i] += counts[i - 1]
    return counts[mass]


def rational_entropy_comparison(value, target):
    # f(p/q) = log_2((p+q)^(p+q)/(p^p q^q))/q.
    p, q = value.numerator, value.denominator
    numerator, denominator = target.numerator, target.denominator
    left = (p + q) ** ((p + q) * denominator)
    right = p ** (p * denominator) * q ** (q * denominator) * 2 ** (q * numerator)
    return (left > right) - (left < right)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    stars_checks = 0
    for mass in range(31):
        for bins in range(1, 25):
            need(composition_count(mass, bins) == math.comb(mass + bins - 1, bins - 1), 'fixed-mass composition count')
            stars_checks += 1
    need(17 ** 17 < 2 ** 70, 'exact row-entropy bound')
    need(64 ** 128 == 2 ** (15 * 64 - 192), 'uniform lower-bound base')
    need(65 ** 128 < 64 ** 128 * 2 ** 15, 'uniform lower-bound induction')
    checked = []
    for n in range(64, 513):
        L = n // 16
        N = n * (n - 1) // 2
        row_upper = math.comb(n + L, L) ** n
        mass_upper = math.comb(N + n * L - 1, n * L - 1)
        strong_labelled_lower = 1 << (N - n)
        need(row_upper * math.factorial(n) < strong_labelled_lower, 'row format fails class separation')
        need(mass_upper * math.factorial(n) < strong_labelled_lower, 'mass format fails class separation')
        checked.append(n)
    F = fractions.Fraction
    brackets = [
        ('row', F(1048, 10000), F(1049, 10000), F(1, 2)),
        ('twice_mass', F(2938, 10000), F(2940, 10000), F(1)),
    ]
    for name, lo, hi, target in brackets:
        need(rational_entropy_comparison(lo, target) < 0, name + ' lower bracket')
        need(rational_entropy_comparison(hi, target) > 0, name + ' upper bracket')
    with decimal.localcontext() as ctx:
        ctx.prec = 65
        D = decimal.Decimal
        log2 = D(2).ln()
        lo, hi = D(0), D(1)
        for _ in range(210):
            x = (lo + hi) / 2
            entropy = ((1 + x) * (1 + x).ln() - x * x.ln()) / log2
            if entropy < 1:
                lo = x
            else:
                hi = x
        constant = format((lo + hi) / 4, '.30f')
    result = {
        'status': 'PASS',
        'scope': 'Necessary class-separation information bound for a single unaugmented outgoing-row multiset signature, even restricted to strong tournaments. Not a bound on full authenticated inventory reconstruction.',
        'fixed_mass_upper_bound': 'binom(n(n-1)/2+nL-1,nL-1)',
        'strong_class_lower_bound': '2^(n(n-3)/2)/n!',
        'necessary_liminf_L_over_n': {'equation': 'f(2 alpha)=1', 'decimal': constant, 'rigorous_rational_bracket': ['0.1469', '0.1470']},
        'stars_and_bars_DP_checks': stars_checks,
        'exact_finite_pigeonhole_checks': {'n_min': min(checked), 'n_max': max(checked), 'orders': len(checked), 'L': 'floor(n/16)', 'both_row_and_mass_upper_bounds': True},
        'all_order_threshold': 'For every n>=64 and integer 1<=L<=n/16, the signature cannot separate all strong classes.',
        'all_order_threshold_proof_checked': '17^17<2^70, exact n=64 equality, and 65^128<2^15*64^128 induction step; not inferred merely from finite sweep.',
    }
    dest = HERE / 'independent_alphabet_counting_receipt.json'
    if args.write:
        dest.write_text(json.dumps(result, indent=2) + '\n')
    else:
        need(json.loads(dest.read_text()) == result, 'frozen receipt mismatch')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
