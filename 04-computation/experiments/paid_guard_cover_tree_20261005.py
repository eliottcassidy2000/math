"""Exact, scoped paid residue-tree refinement inside literal n = 27 (mod 64).

Standard library only.  Run normally or with -O; checks are not assertions.
The infinite claims are proved in the companion note.  No orbit-to-root
oracle is used, and an open tree leaf is not a nonconvergence claim.
"""
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math


ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(message)


def positive_odd(n):
    if type(n) is not int or n <= 0 or n % 2 != 1:
        raise ValueError("source must be an exact positive odd integer")


def valuation(n):
    if type(n) is not int or n <= 0:
        raise ValueError("valuation requires an exact positive integer")
    return (n & -n).bit_length() - 1


def step(n):
    positive_odd(n)
    z = 3*n + 1
    a = valuation(z)
    return z // 2**a, a


def coefficients(word):
    p, q, b = 1, 1, 0
    for a in word:
        if type(a) is not int or a < 1:
            raise ValueError("word letters must be exact positive integers")
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def exact_cell(word):
    p, q, b = coefficients(word)
    return ((q-b)*pow(p, -1, 2*q)) % (2*q), 2*q


def coarse_cell(prefix, minimum):
    """Exact prefix followed by a final valuation at least minimum."""
    p, q, b = coefficients(tuple(prefix) + (minimum,))
    return (-b*pow(p, -1, q)) % q, q


def rank(n):
    positive_odd(n)
    if n == 1:
        return 0, 0
    a = valuation(n-1)
    return 3**a*((n-1)//2**a)**2, a


def paid_fifth(n):
    """Return the actual fifth endpoint/exponent, or None outside this guard.

    This certifies a paid dependency, not a complete certificate to root 1.
    """
    positive_odd(n)
    if n % 256 != 219:
        return None
    t = (n-219)//256
    z = 209 + 243*t
    extra = valuation(z)
    return z // 2**extra, 3+extra


@dataclass(frozen=True)
class Leaf:
    label: str


@dataclass(frozen=True)
class Split:
    zero: object
    one: object


# Addresses read the bits of t=(n-27)/64 from least significant to most.
TREE = Split(Leaf("fifth=1"),
             Split(Leaf("fifth=2"),
                   Split(Leaf("fifth=3: new"), Leaf("fifth>=4: old"))))


def leaf_cells(tree, residue=0, bits=0):
    if isinstance(tree, Leaf):
        return [(residue, 2**bits, tree.label)]
    if not isinstance(tree, Split):
        raise ValueError("invalid typed tree")
    return (leaf_cells(tree.zero, residue, bits+1) +
            leaf_cells(tree.one, residue+2**bits, bits+1))


def evaluate(tree, parameter):
    if type(parameter) is not int or parameter < 0:
        raise ValueError("parameter must be an exact nonnegative integer")
    if isinstance(tree, Leaf):
        return tree.label
    if not isinstance(tree, Split):
        raise ValueError("invalid typed tree")
    return evaluate(tree.one if parameter % 2 else tree.zero, parameter//2)


def locate(tree, parameter):
    """Return the actual low-bit address and the unused high-bit parameter."""
    if type(parameter) is not int or parameter < 0:
        raise ValueError("parameter must be an exact nonnegative integer")
    address = []
    while isinstance(tree, Split):
        bit = parameter % 2
        address.append(bit)
        tree = tree.one if bit else tree.zero
        parameter //= 2
    if not isinstance(tree, Leaf):
        raise ValueError("invalid typed tree")
    return tree.label, tuple(address), parameter


def restore(address, quotient):
    if type(quotient) is not int or quotient < 0 or any(type(b) is not int or b not in (0, 1) for b in address):
        raise ValueError("invalid bit address or quotient")
    return sum(b*2**j for j, b in enumerate(address)) + 2**len(address)*quotient


def disjoint(c, d):
    return (c[0]-d[0]) % math.gcd(c[1], d[1]) != 0


def replay_controls():
    check(coefficients((1, 2, 1, 1)) == (81, 32, 85), "forced prefix carry")
    check(exact_cell((1, 2, 1, 1)) == (27, 64), "literal source domain")
    check(coefficients((1, 2, 1, 1, 3)) == (243, 256, 287), "five-step carry")
    check(exact_cell((1, 2, 1, 1, 3)) == (219, 512), "new exact half")
    check(coarse_cell((1, 2, 1, 1), 3) == (219, 256), "merged paid tail")
    check(coarse_cell((1, 2, 1, 1), 4) == (475, 512), "old switched tail")
    expected = [(0, 2, "fifth=1"), (1, 4, "fifth=2"),
                (3, 8, "fifth=3: new"), (7, 8, "fifth>=4: old")]
    check(leaf_cells(TREE) == expected, "typed low-bit decoder")
    check(sum(Fraction(1, m) for _, m, _ in expected) == 1, "complete local partition")
    for i, c in enumerate(expected):
        for d in expected[:i]:
            check(disjoint(c, d), "pairwise disjoint local leaves")
    for t in range(4096):
        n = 27+64*t
        x = n
        for a in (1, 2, 1, 1):
            x, actual = step(x)
            check(a == actual and x > n, "forced proper prefix grows")
        child, a = step(x)
        label = ("fifth=1" if a == 1 else "fifth=2" if a == 2 else
                 "fifth=3: new" if a == 3 else "fifth>=4: old")
        check(evaluate(TREE, t) == label, "literal/tree independent agreement")
        located_label, address, quotient = locate(TREE, t)
        check(located_label == label and restore(address, quotient) == t,
              "address plus unused bits is lossless")
        if a >= 3:
            check(paid_fifth(n) == (child, a), "compressed actual endpoint")
            check(0 < child < n and rank(child) < rank(n), "original-source payment")
        else:
            check(child > n, "two low-reset leaves still grow at step five")
    for t in range(1024):
        n = 219+512*t
        check(paid_fifth(n) == (209+486*t, 3), "whole exact affine progression")
        check(n-(209+486*t) == 10+26*t, "strict affine decrease")
    check(paid_fifth(27) is None and paid_fifth(91) is None, "unpaid-by-this-tree controls")
    check(paid_fifth(155) is None and (729*155+669)//1024 == 111,
          "tree-open does not mean uncovered by the inherited checkpoint")
    swapped = Split(TREE.one, TREE.zero)
    check(sorted(label for _, _, label in leaf_cells(swapped)) ==
          sorted(label for _, _, label in expected) and evaluate(swapped, 0) != evaluate(TREE, 0),
          "leaf label counts alone lose guard placement")
    check(step(3)[1] == 1 and exact_cell((1, 2, 1, 3)) == (187, 256),
          "old four-slot source cell retained")
    for bad in (True, 3.0, 0, -1, 2):
        try:
            paid_fifth(bad)
        except ValueError:
            pass
        else:
            raise ArithmeticError("source type guard")
    return [dict(parameter_residue=c, parameter_modulus=m, label=label,
                 source_residue=27+64*c, source_modulus=64*m)
            for c, m, label in expected]


def bank_controls():
    paths = [ROOT/'05-knowledge/results'/name for name in (
        'collatz_paid_portrait_controllers_20261004.json',
        'collatz_four_slot_compression_20261004.json',
        'collatz_binary_ternary_guard_fusion_20261004.json')]
    portrait, four, binary = [json.loads(path.read_text()) for path in paths]
    cells = portrait['rational_bank']['rows'] + portrait['switched_bank']['rows']
    cells += [row for row in portrait['safe_cylinders'] if row['c'] == 17]
    cells += [portrait['switched_bank']['checkpoint']]
    cells += four['permuted_bank']['rows']
    cells += [dict(residue=187, modulus=256)]
    for row in cells:
        check(disjoint((219, 512), (row['residue'], row['modulus'])),
              "new exact half disjoint from saved old dyadic rows")
    check([(row['s'], row['source_word'][:-1]) for row in binary['debt_rows'] if row['e'] == 1]
          == [(1, [1, 6]), (2, [6])], "binary16 e1 obstruction prefixes retained")
    old_switch = [r for r in portrait['switched_bank']['rows'] if (r['q'], r['r']) == (1, 2)]
    check(len(old_switch) == 1 and
          (old_switch[0]['L'], old_switch[0]['residue'], old_switch[0]['modulus']) == (4, 475, 512),
          "frozen old factor-two threshold")
    previous = tuple(Fraction(x) for x in four['coverage']['exponent_3_mod8_coverage_interval'])
    subtotal = tuple(x+Fraction(1, 16) for x in previous)
    return dict(saved_rows_checked=len(cells),
                provenance_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                exponent_increment="1/16", exponent_subtotal_interval=list(map(str, subtotal)),
                exponent_subtotal_decimal=float(subtotal[0]),
                scope="Intermediate subtotal over this frozen bank; not the current larger threshold-bank union.")


def power_controls():
    check([k for k in range(128) if pow(3, k, 512) == 219] == [83], "exact new exponent phase")
    check([k for k in range(64) if pow(3, k, 256) == 219] == [19], "coarse exponent phase")
    check([k for k in range(64) if pow(3, k, 256) == 187] == [27], "old exponent27 source distinction")
    check(pow(3, 128, 512) == 1 and pow(3, 64, 512) != 1, "exact dyadic power order")
    check(pow(3, 27, 64) == 59 and exact_cell((1, 2, 1, 1)) == (27, 64),
          "literal and exponent27mod64 must not be identified")
    for j in range(9):
        n = 3**(83+128*j)
        x = n
        for a in (1, 2, 1, 1, 3):
            x, actual = step(x)
            check(a == actual, "literal powers retain exact word")
        check(x == (243*n+287)//256 and 0 < x < n, "power source identity retained")
    return dict(new_exponents="83 mod128", merged_exponents="19 mod64",
                inherited_half="19 mod128", inherited_four_slot="27 mod64",
                literal_power_controls=9, maximum_exponent=83+128*8)


def main():
    tree = replay_controls()
    coverage = bank_controls()
    powers = power_controls()
    print("PROVED: exact guarded first descent on n219mod256; new frozen-bank half n219mod512.")
    print("LOCAL TREE: " + json.dumps(tree, sort_keys=True))
    print("NEW AFFINE ROW: n=219+512t ->209+486t; word12113; (P,Q,B)=(243,256,287).")
    print("COARSE ROW: n=219+256t ->oddpart(209+243t); fifth valuation=3+v2(209+243t).")
    print("POWERS: " + json.dumps(powers, sort_keys=True))
    print("FROZEN COMPARISON: " + json.dumps(coverage, sort_keys=True))
    print("FINITE UNIVERSE: 4096 local tree parameters; 1024 exact-half parameters; 9 source powers.")
    print("HOSTILES: literal27 differs from exponent27; source155 is tree-open but old-checkpoint paid.")
    print("SCOPE: paid dependencies only; other local leaves grow at step five; universal root coverage OPEN.")
    print("CHECKS: " + str(CHECKS))


if __name__ == '__main__':
    main()
