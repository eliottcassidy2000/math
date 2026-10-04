"""Independent Fibonacci-digit verification of fixed arithmetic guards.

Run: python3 04-computation/experiments/zeckendorf_guard_automaton_20261003.py
All arithmetic is exact. No external packages, tables, or inherited filters.
The owner confirmed red/black/blue diagonals with extra copies of 1. Their
typed representation reader is separate from the auxiliary XOR charge.
The supplied finite stripe does not yet fix a unique infinite row rule.
"""

from importlib.util import module_from_spec, spec_from_file_location
from collections import Counter
from itertools import product
from math import isqrt
from pathlib import Path


SUPPLIED_STRIPE = "RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB"
MODULI = (3, 9, 27, 81, 4, 16, 64)


def zeckendorf_digits(n):
    """Canonical digits, high weight first; weights low-to-high are 1,2,3,..."""
    if not isinstance(n, int) or n < 0:
        raise ValueError("source must be a nonnegative integer")
    if n == 0:
        return (0,)
    weights = [1, 2]
    while weights[-1] <= n:
        weights.append(weights[-1] + weights[-2])
    digits = []
    for weight in reversed(weights[:-1]):
        digit = int(weight <= n)
        digits.append(digit)
        n -= digit * weight
    assert n == 0
    return tuple(digits)


def beatty_register(n):
    m = n + 1
    return (3 * m - isqrt(5 * m * m) - 1) // 2


def charge(n):
    return n % 2, beatty_register(n) % 2


def xor(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def digit_automaton(digits, modulus):
    """State is value X, Fibonacci shift Y, previous digit, and XOR charge."""
    if modulus < 1:
        raise ValueError("modulus must be positive")
    X = Y = previous = ca = cb = 0
    for digit in digits:
        if digit not in (0, 1) or (previous and digit):
            raise ValueError("not a binary Fibonacci normal form")
        X, Y = (Y + digit) % modulus, (X + Y + 2 * digit) % modulus
        ca, cb = cb ^ digit, ca ^ cb
        previous = digit
    return X, Y, previous, (ca, cb)


def direct_digit_value(digits, shifted=False):
    a, b = (2, 3) if shifted else (1, 2)
    result = 0
    for digit in reversed(digits):
        result += digit * a
        a, b = b, a + b
    return result


def wythoff_candidate(n):
    if n < 1:
        raise ValueError("stripe positions start at 1")
    trailing = next(i for i, digit in enumerate(reversed(zeckendorf_digits(n))) if digit)
    return "R" if trailing == 0 else "K" if trailing % 2 else "B"


def valuation_prefix(n, count):
    """Direct arithmetic path, independent of the route encoding functions."""
    valuations = []
    for _ in range(count):
        n = 3 * n + 1
        k = 0
        while n % 2 == 0:
            n //= 2
            k += 1
        valuations.append(k)
    return tuple(valuations)


def ordinary_support(n):
    return frozenset(i + 2 for i, d in enumerate(reversed(zeckendorf_digits(n))) if d)


def extra_unit_fibre(n):
    if n <= 0:
        raise ValueError("fibre formula here is for positive integers")
    before = ordinary_support(n - 1)
    result = [ordinary_support(n), before | {0}]
    if 2 not in before:
        result.append(before | {1})
    return tuple(result)


def extra_unit_reader(support, modulus):
    """Read 1_K,1_B,1_R,F3,F4,... without identifying labels with charge."""
    support = frozenset(support)
    if any(i < 0 or i + 1 in support for i in support):
        raise ValueError("extended atoms must occupy nonadjacent indices")
    highest = max(support | {2})
    digits = tuple(int(i in support) for i in range(highest, 1, -1))
    X, Y, red, color = digit_automaton(digits, modulus)
    marker = "K" if 0 in support else "B" if 1 in support else "ordinary"
    defect = 0
    if marker != "ordinary":
        if marker == "B":
            assert red == 0
        X, Y = (X + 1) % modulus, (Y + 2 - red) % modulus
        color = color[0] ^ 1, color[1] ^ red
        defect = 1 - red if marker == "K" else 0
    return X, Y, color, marker, defect


def multiply_root_pair(z, w):
    """Products in Z[omega], omega^2+omega+1=0."""
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c - b * d


def color_fourier(vector):
    powers = ((1, 0), (0, 1), (-1, -1))
    return tuple(tuple(sum(vector[c] * powers[(-r * c) % 3][k] for c in range(3))
                       for k in range(2)) for r in range(3))


def load_sibling_module(filename, name):
    path = Path(__file__).with_name(filename)
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    controls = crt_controls = 0
    for n in range(10001):
        digits = zeckendorf_digits(n)
        shift = direct_digit_value(digits, shifted=True)
        assert direct_digit_value(digits) == n
        assert shift == 2 * n - beatty_register(n)
        for modulus in MODULI:
            X, Y, _, color = digit_automaton(digits, modulus)
            assert (X, Y) == (n % modulus, shift % modulus)
            assert color == charge(n)
            if modulus % 2:
                xx, yy, _, cc = digit_automaton(digits, 2 * modulus)
                assert (xx % modulus, yy % modulus) == (X, Y)
                assert (xx % 2, yy % 2) == color == cc
                crt_controls += 1
            else:
                assert (X % 2, Y % 2) == color
            controls += 1
    try:
        digit_automaton((1, 1), 3)
    except ValueError:
        pass
    else:
        raise AssertionError("adjacent occupied digits must be rejected")
    print("Universe: every n=0..10000; moduli 3,9,27,81,4,16,64")
    print(f"Exact digit/value/shift/charge checks: {controls}; odd-modulus CRT checks: {crt_controls}")
    print("Hostile normal form 11 rejected; zero and leading-zero states allowed")

    states = [digit_automaton(zeckendorf_digits(n), 3) for n in (2, 8)]
    assert states == [(2, 0, 0, (0, 1)), (2, 1, 0, (0, 1))]
    shifted_residues = [digit_automaton(zeckendorf_digits(n) + (0,), 3)[0] for n in (2, 8)]
    assert shifted_residues == [0, 1]
    print("Shift sidecar hostile: 2 and 8 share value mod3, charge, and last digit; shifts 3 and 13 split")

    assert len(SUPPLIED_STRIPE) == 35
    candidate = "".join(wythoff_candidate(n) for n in range(1, 36))
    mismatches = [i + 1 for i, (a, b) in enumerate(zip(SUPPLIED_STRIPE, candidate)) if a != b]
    assert mismatches == [35] and SUPPLIED_STRIPE[-1] == "B" and candidate[-1] == "R"
    assert charge(2) == charge(4) == (0, 1)
    assert SUPPLIED_STRIPE[1] == "K" and SUPPLIED_STRIPE[3] == "R"
    assert charge(6) == (0, 0)
    print("Supplied stripe:", SUPPLIED_STRIPE)
    print("Wythoff candidate mismatch: position35, supplied B / candidate R")
    print("XOR charge differs from supplied stripe: 2 and4 share (0,1), but are K/R; neutral at6")

    # Complete independent subset enumeration, not generated by the fibre formula.
    fib = [0, 1]
    while len(fib) < 14:
        fib.append(fib[-1] + fib[-2])
    weights = [1, 1] + fib[2:]
    enumerated = {n: set() for n in range(1, 201)}
    for mask in range(1, 1 << len(weights)):
        if mask & (mask << 1):
            continue
        support = frozenset(i for i in range(len(weights)) if (mask >> i) & 1)
        value = sum(weights[i] for i in support)
        if 1 <= value <= 200:
            enumerated[value].add(support)
    # The next omitted Fibonacci weight is 377>200, so this universe is complete.
    extended_checks = representations = 0
    for n in range(1, 201):
        fibre = extra_unit_fibre(n)
        assert set(fibre) == enumerated[n]
        r = int(2 in ordinary_support(n))
        before_r = int(2 in ordinary_support(n - 1))
        polynomial = Counter((int(0 in s), int(1 in s), int(2 in s)) for s in fibre)
        expected = Counter({(0, 0, r): 1, (1, 0, before_r): 1})
        if not before_r:
            expected[(0, 1, 0)] += 1
        assert polynomial == expected
        charge_preserving = 0
        for support in fibre:
            raw_a = raw_b = 0
            for i in support:
                atom_charge = (-1, 1) if i == 0 else (1, 0) if i == 1 else (
                    weights[i] - 2 * beatty_register(weights[i]), beatty_register(weights[i]))
                raw_a += atom_charge[0]
                raw_b += atom_charge[1]
            q = n - 2 * beatty_register(n), beatty_register(n)
            charge_preserving += (raw_a, raw_b) == q
            for modulus in (3, 9, 27):
                X, Y, color, marker, defect = extra_unit_reader(support, modulus)
                assert (X, Y) == (n % modulus, (2 * n - beatty_register(n)) % modulus)
                assert color == charge(n)
                assert (raw_a - q[0], raw_b - q[1]) == (-2 * defect, defect)
                recovered = ordinary_support(n) if marker == "ordinary" else (
                    ordinary_support(n - 1) | ({0} if marker == "K" else {1}))
                assert recovered == support
                extended_checks += 1
            representations += 1
        assert charge_preserving == 2
    words = ["K", "B", "R", "RK", "RKB"]
    while len(words) < 14:
        words.append(words[-1] + words[-2])
    expand = lambda support: "".join(words[i] for i in sorted(support, reverse=True))
    assert expand(frozenset((9, 1))) == SUPPLIED_STRIPE
    row36 = {expand(support) for support in extra_unit_fibre(36)}
    assert len(row36) == 1 and all(word[34] == "R" for word in row36)
    assert frozenset((0, 2)) in extra_unit_fibre(2)
    try:
        extra_unit_reader((0, 1), 3)
    except ValueError:
        pass
    else:
        raise AssertionError("cyclic relabelling of the units is not an index-legality symmetry")
    print(f"Confirmed extra-unit convention: all fibres n=1..200, {representations} representations, {extended_checks} native reader checks")
    print("Marker polynomial and exactly two charge-preserving sheets pass; (value,marker) recovers every support")
    print("Supplied35-blue row = atom34 + blue1; both legal row36 supports force red at position35 in this grammar")
    print("Unit-label rotation is not an index symmetry: black+red is legal, black+blue is not")
    powers = ((1, 0), (0, 1), (-1, -1))
    for label in SUPPLIED_STRIPE:
        vector = tuple(int(label == c) for c in "KBR")
        transformed = color_fourier(vector)
        for c in range(3):
            terms = [multiply_root_pair(transformed[r], powers[(r * c) % 3]) for r in range(3)]
            assert tuple(sum(term[k] for term in terms) for k in range(2)) == (3 * vector[c], 0)
    blue, red = color_fourier((0, 1, 0)), color_fourier((0, 0, 1))
    assert blue != red
    assert [2 * a - b for a, b in blue] == [2 * a - b for a, b in red]
    assert [a * a - a * b + b * b for a, b in blue] == [a * a - a * b + b * b for a, b in red]
    print("Actual-color C3 Fourier basis: all35 positionwise inversions exact; cosine-only and magnitudes conflate blue/red")

    defects = set()
    diagonal_checks = 0
    for a in range(1, 201):
        for b in range(1, 201):
            delta = beatty_register(a) + beatty_register(b) - beatty_register(a + b)
            defects.add(delta)
            assert xor(xor(charge(a), charge(b)), (0, delta % 2)) == charge(a + b)
            diagonal_checks += 1
    assert defects == {-1, 0, 1}
    assert xor(charge(1), charge(4)) == (1, 1)
    assert xor(charge(2), charge(3)) == (1, 0)
    assert charge(5) == (1, 0)
    print(f"Diagonal carry: {diagonal_checks} positive pairs; defects -1,0,1; strict sum5 splits before repair")

    routes = load_sibling_module("collatz_route_tournaments_20261003.py", "route_certificates")
    valid = rejected = macro_guards = negative_guards = valid_changed_heads = 0
    for length in range(1, 6):
        for word in product(range(5), repeat=length):
            try:
                source = routes.decode_odd_word(word)
            except ValueError:
                rejected += 1
                continue
            valid += 1
            for p, j in enumerate(word):
                modulus = 3 ** p
                phase = digit_automaton(zeckendorf_digits(j), modulus)[0]
                changed_j = j + modulus
                changed_phase = digit_automaton(zeckendorf_digits(changed_j), modulus)[0]
                assert phase == changed_phase == j % modulus
                changed = word[:p] + (changed_j,) + word[p + 1:]
                new_source = routes.decode_odd_word(changed)
                assert valuation_prefix(source, p) == valuation_prefix(new_source, p)
                macro_guards += 1
                if p:
                    wrong_j = j + 1
                    assert digit_automaton(zeckendorf_digits(wrong_j), modulus)[0] != phase
                    wrong_word = word[:p] + (wrong_j,) + word[p + 1:]
                    try:
                        wrong_source = routes.decode_odd_word(wrong_word)
                    except ValueError:
                        pass
                    else:
                        assert valuation_prefix(source, p) != valuation_prefix(wrong_source, p)
                        valid_changed_heads += 1
                    negative_guards += 1
    assert (valid, rejected, macro_guards, negative_guards) == (694, 3211, 3174, 2480)
    print("Route universe: j=0..4, lengths1..5; 694 valid, 3211 rejected")
    print(f"Fibonacci-digit macro phase gates: {macro_guards}; direct prefix valuations agree")
    print(f"Wrong-phase controls: {negative_guards}; {valid_changed_heads} remain valid but change the head")

    guards = load_sibling_module("collatz_fourier_guards_20261003.py", "compiled_guards")
    compiled_cases = root_exceptions = 0
    for p in range(1, 5):
        modulus = 3 ** p
        for prefix in product(range(3), repeat=p):
            phases = guards.compile_guard(prefix, target=1)
            assert len(phases) == 2 ** p
            for index in range(2 * modulus):
                phase = digit_automaton(zeckendorf_digits(index), modulus)[0]
                arithmetic_accepts = phase in phases
                predicted = arithmetic_accepts and index > 0
                try:
                    source = routes.decode_odd_word(prefix + (index,))
                except ValueError:
                    actual = False
                else:
                    actual = True
                    assert valuation_prefix(source, p) == phases[phase]
                assert actual == predicted
                root_exceptions += int(arithmetic_accepts and index == 0)
                compiled_cases += 1
    assert (compiled_cases, root_exceptions) == (14760, 45)
    assert guards.compile_guard((1,))[0] == (4,)
    assert 0 in guards.compile_guard((1,))
    assert routes.decode_odd_word((1, 3)) == 453
    print(f"Compiler/DFA/strict-route interface: {compiled_cases} cases; {root_exceptions} isolated index-zero padding exclusions")
    print("First-hit guard removes t=0 only, not its whole residue class: head(1), t=3 is valid source453")
    print("Scope: fixed-modulus finite automata on unbounded digit words; no finite color clock or universal coverage claimed")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
