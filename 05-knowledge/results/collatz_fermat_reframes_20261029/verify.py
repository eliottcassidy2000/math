"""Exact finite audits for README.md; no assert statements or conjecture solver.

Run: python -B path/to/verify.py --output path/to/results.json
Dependencies: Python 3.10+, SymPy (polynomial arithmetic and small factorization).
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
import importlib.util
from math import comb, gcd, prod
import argparse
import json
from pathlib import Path

import sympy as sp

GATES = Counter()


def check(ok, name):
    GATES[name] += 1
    if not ok:
        raise ArithmeticError(name)


def valuation2(n):
    if n == 0:
        raise ValueError("zero has no finite valuation")
    n = abs(n)
    return (n & -n).bit_length() - 1


def step(a, x, b=1):
    y = a * Q(x) + b
    v = valuation2(y.numerator) - valuation2(y.denominator)
    check(v >= 1, "odd rational domain")
    return y / (1 << v), v


def numerator(a, word):
    c, h = 0, 0
    for v in word:
        c = a * c + (1 << h)
        h += v
    return c


def fixed(a, word):
    return Q(numerator(a, word), (1 << sum(word)) - a ** len(word))


def rotations(word):
    return [word[i:] + word[:i] for i in range(len(word))]


def primitive(word):
    k = len(word)
    return all(word != word[:d] * (k // d)
               for d in range(1, k) if k % d == 0)


def compositions(h, k):
    for cuts in combinations(range(1, h), k - 1):
        ends = (0,) + cuts + (h,)
        yield tuple(ends[i + 1] - ends[i] for i in range(k))


def marked_count(a, b, k, h):
    """Residue DP independent of rational-cycle enumeration."""
    modulus = (1 << h) - a ** k
    if modulus <= 0:
        return 0
    state = {(0, 0): 1}
    for j in range(k):
        nxt = defaultdict(int)
        for (height, residue), count in state.items():
            new_residue = (a * residue + (1 << height)) % modulus
            for v in range(1, h - height - (k - j - 1) + 1):
                nxt[height + v, new_residue] += count
        state = nxt
    return sum(count for (height, r), count in state.items()
               if height == h and b * r % modulus == 0)


def cycle_count(a, b, k, h):
    m = sum(int(sp.mobius(d)) * marked_count(a, b, k // d, h // d)
            for d in sp.divisors(gcd(k, h)))
    check(m % k == 0 and m >= 0, "necklace Mobius integrality")
    return m // k


def audit_cycles():
    examples = []
    for word, expected5, expected3 in [
        ((1, 1, 5), (13, 33, 83), (Q(19, 101), Q(79, 101), Q(169, 101))),
        ((1, 3, 3), (17, 43, 27), (Q(31, 101), Q(97, 101), Q(49, 101))),
        ((1, 4), (1, 3), (Q(5, 23), Q(19, 23))),
    ]:
        x3, x5 = [fixed(3, w) for w in rotations(word)], [fixed(5, w) for w in rotations(word)]
        check(tuple(x3) == expected3 and tuple(x5) == expected5, "explicit transports")
        d = x3[0].denominator
        integer_forced = [int(d * x) for x in x3]
        for i, x in enumerate(integer_forced):
            y, v = step(3, x, d)
            check(y == integer_forced[(i + 1) % len(word)] and v == word[i], "scaled forcing")
        examples.append({"word": word, "five": list(expected5), "three": list(map(str, x3)), "forcing": d})

    spectrum = defaultdict(Counter)
    native = {3: [], 5: []}
    total_words = primitive_necklaces = cone_tests = 0
    for k in range(1, 6):
        for h in range(k, 21):
            for word in compositions(h, k):
                total_words += 1
                if word != min(rotations(word)) or not primitive(word):
                    continue
                primitive_necklaces += 1
                rows = {}
                for a in (3, 5, 7):
                    points = [fixed(a, w) for w in rotations(word)]
                    x = points[0]
                    for i, v in enumerate(word):
                        y, actual_v = step(a, points[i])
                        check(y == points[(i + 1) % k] and v == actual_v, "rational phase recurrence")
                        check(points[i].denominator == x.denominator, "rotation denominator")
                    check(gcd(x.denominator, 2 * a) == 1, "denominator prime to multiplier")
                    check(prod(a + 1 / t for t in points) == 1 << h, "cycle product")
                    rows[a] = points
                    if a in (3, 5) and x > 0:
                        spectrum[a, k, h][x.denominator] += 1
                        if x.denominator == 1:
                            native[a].append({"word": word, "cycle": [int(t) for t in points]})
                for a, b in ((3, 5), (3, 7), (5, 7)):
                    if rows[b][0] > 0:
                        cone_tests += 1
                        check(min(rows[a]) < Q(1, b - a), "multiplier cone separation")
                    if rows[a][0] > 0 and rows[a][0].denominator == 1:
                        check(rows[b][0] < 0, "integer cone sign reversal")

    check(total_words == sum(comb(20, k) for k in range(1, 6)), "composition universe")
    check(len(native[3]) == 1 and len(native[5]) == 3, "native controls in finite box")
    # Independent exact recurrence, modest universe to keep it quick.
    dp_cases = 0
    for a, b, k in product((3, 5), (1, 3, 5, 7, 23, 101, 303), range(1, 5)):
        for h in range(k, 13):
            expected = sum(c for d, c in spectrum[a, k, h].items() if b % d == 0)
            check(cycle_count(a, b, k, h) == expected, "residue DP vs conductor divisor sum")
            dp_cases += 1
    # Invert the forcing divisor relation at all divisors of 303, for each clock.
    for d in sp.divisors(303):
        for k in range(1, 5):
            for h in range(k, 13):
                inverted = sum(int(sp.mobius(d // e)) * cycle_count(3, e, k, h)
                               for e in sp.divisors(d))
                check(inverted == spectrum[3, k, h][d], "forcing Mobius inversion")
    # Same k,H,denominator but different marked and unmarked cycles: ports matter.
    check(fixed(3, (1, 1, 5)).denominator == fixed(3, (1, 3, 3)).denominator == 101,
          "denominator collision")
    check(set(rotations((1, 1, 5))).isdisjoint(rotations((1, 3, 3))), "word selector needed")
    attached = [fixed(3, w + (1,)) for w in ((1, 1, 5), (1, 3, 3))]
    check(attached == [Q(37, 35), Q(221, 175)], "one shared attachment separates denominators")
    return {"universe": "positive valuation words k<=5, H<=20; primitive rotation orbits",
            "all_words": total_words, "primitive_necklaces": primitive_necklaces,
            "cone_tests": cone_tests, "dp_cases": dp_cases, "native_cycles": native,
            "explicit_transports": examples,
            "denominator_101_cycles_k3_H7": spectrum[3, 3, 7][101],
            "common_one_step_attachment": list(map(str, attached))}


def audit_fuel():
    word = (1, 1, 5)
    x, h = fixed(3, word), sum(word)
    n, d = x.numerator, x.denominator
    rows = []
    for repetitions in range(1, 21):
        modulus = 1 << (h * repetitions + 1)
        residue = n * pow(d, -1, modulus) % modulus
        carry = (d * residue - n) // modulus
        y = Q(residue)
        for expected in word * repetitions:
            y, v = step(3, y)
            check(v == expected, "finite prefix exact clock")
        check(residue > 0, "all finite prefixes inhabited")
        if rows:
            check(carry % d == rows[-1]["carry"] * pow(1 << h, -1, d) % d, "carry phase clock")
        rows.append({"repetitions": repetitions, "least_source": residue, "carry": carry})
    return {"period_mod_101": int(sp.n_order(pow(2, 7, 101), 101)), "rows": rows}


def audit_fermat():
    rows = []
    for n in range(11):
        m = 1 << n
        f = (1 << m) + 1
        check(pow(2, m, f) == f - 1, "base two half clock")
        # Base 2 belongs to the Hall-projection kernel, including F0=3.
        check(pow(2, f - 1, f) == 1, "base two blindness")
        if n >= 1:
            check(gcd(3, f) == 1, "base three unit")
            r = pow(3, f - 1, f)
            check((r == 1) == (n <= 4), "base three finite Fermat certificates")
            # Independent repeated-squaring implementation of the exponent 2^m.
            sq = 3
            for _ in range(m):
                sq = sq * sq % f
            check(r == sq, "independent modular powering")
            rows.append({"n": n, "F_bits": f.bit_length(), "three_projection": str(r)})
    factors = (641, 6700417)
    f5 = (1 << 32) + 1
    check(prod(factors) == f5 and all(sp.isprime(p) for p in factors), "F5 factor certificate")
    odd_hall = prod((p - 1) >> valuation2(p - 1) for p in factors)
    kernel = prod(1 << valuation2(p - 1) for p in factors)
    check(odd_hall == 261735 and kernel == 16384, "F5 Hall orders")
    x = 3303966555
    check(x != 1 and pow(x, 3, f5) == 1, "inherited explicit order three")
    # Generic prime-power hostile cases: no squarefreeness assumption.
    for modulus in (9, 25, 45, 65, 85, 91, 257):
        units = [a for a in range(1, modulus) if gcd(a, modulus) == 1]
        decomposition = sp.factorint(modulus)
        exponent = max(valuation2(int(p) - 1) for p in decomposition)
        kernel_size = prod(1 << valuation2(int(p) - 1) for p in decomposition)
        image = Counter(pow(a, 1 << exponent, modulus) for a in units)
        expected_odd = prod(int(p) ** (int(e) - 1) * ((int(p) - 1) >> valuation2(int(p) - 1))
                            for p, e in decomposition.items())
        check(len(image) == expected_odd and set(image.values()) == {kernel_size}, "Hall fibers including prime powers")
        check(all(pow(a, expected_odd, modulus) == 1 for a in image), "odd image orders")
    return {"finite_certificates": rows, "F5_odd_Hall": odd_hall, "F5_kernel": kernel,
            "F5_fraction_of_units_detected": str(Q(odd_hall - 1, odd_hall)),
            "uniform_base_three_claim": "OPEN: not inferred from this finite box"}


def audit_charts():
    X = sp.Symbol("X")
    checked = successes = 0
    for s in range(3):
        degree = 1 << s
        fpoly = sp.Poly(X ** degree + 1, X)
        for p, q, d in product((X, X + 1, X - 1, X ** 2 + X + 1, 3, 3303966555),
                               (sp.Integer(1), X + 1, X ** 2 + 1), (3, 5, 7)):
            P, Qp = sp.Poly(p, X), sp.Poly(q, X)
            remainder = (P ** d - Qp ** d).rem(fpoly)
            size = sum(abs(int(c)) for c in remainder.all_coeffs())
            if remainder.is_zero and not Qp.rem(fpoly).is_zero:
                check((P - Qp).rem(fpoly).is_zero, "cyclotomic odd root rigidity samples")
            for n in range(s, 9):
                t = 1 << (1 << (n - s))
                f = t ** degree + 1
                pv, qv = int(P.eval(t)), int(Qp.eval(t))
                if gcd(qv, f) != 1:
                    continue
                checked += 1
                residue = pv * pow(qv, -1, f) % f
                success = residue != 1 and pow(residue, d, f) == 1
                successes += int(success)
                if success:
                    check(not remainder.is_zero and t <= 1 + size, "fixed chart effective cutoff")
                if t > 1 + size:
                    check(not success, "fixed chart exclusion above cutoff")
        # Exact even-order counterexample to any attempted removal of oddness.
        for n in range(s, 9):
            t = 1 << (1 << (n - s))
            modulus = t ** degree + 1
            check(pow(t, degree, modulus) == modulus - 1 and pow(t, 2 * degree, modulus) == 1,
                  "even clock hostile control")
    return {"admissible_fixed_chart_samples": checked, "nontrivial_odd_roots_in_samples": successes,
            "theorem_is_symbolic": True, "sample_orders": [3, 5, 7], "sample_lags": [0, 1, 2]}


def audit_inherited_hostile():
    path = Path(__file__).resolve().parent.parent / "four_vertex_glue_fermat_20261029" / "verify.py"
    spec = importlib.util.spec_from_file_location("inherited_four_vertex", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.hostile_control()
    check(result["automorphisms"] == [81, 27] and result["two_run_mod3_ranks"] == [24, 22],
          "inherited tournament hostile control")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = {"cycles": audit_cycles(), "repeated_word_fuel": audit_fuel(),
            "fermat": audit_fermat(), "charts": audit_charts(), "hostile_81_27": audit_inherited_hostile()}
    data["gates"] = dict(sorted(GATES.items()))
    data["total_gates"] = sum(GATES.values())
    output = json.dumps(data, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(json.dumps({"total_gates": data["total_gates"], "cycles": data["cycles"],
                      "charts": data["charts"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
