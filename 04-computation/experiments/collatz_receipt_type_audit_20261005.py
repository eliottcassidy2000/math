#!/usr/bin/env python3
"""Exact hostile controls for the concurrent four-receipt-question note."""
import ast
from fractions import Fraction
from pathlib import Path
import json
from collatz_fusion_helpers_20261005 import U, v2, value, boundary, D, K, add, prime_image


def main():
    checks = 0

    def check(test, message):
        nonlocal checks
        checks += 1
        if not test:
            raise ValueError(message)

    packet = {3: 1, 13: -1}
    check(U(3) == U(13) == 5, "common future")
    check(value(packet) == Fraction(3, 13), "signed packet is not a weak receipt for3")
    error = D(3, packet)
    check(error == {1: 1, 13: -1}, "labelled open endpoint")
    check(K(13) == {} and prime_image(error) == {13: -1}, "prime child still carries charge")
    check(v2(3*85+1) == 8 and U(U(85)) == U(5) == 1, "excluded two-step join")
    check(8+v2(3*U(85)+1)-v2(3*5+1) == 6, "old A bound fails")
    wild5 = {7: 2, 11: 1, 17: 1, 55: 1, 65: 1, 83: 1}
    check(value(wild5) == Fraction(1, 5), "wild reciprocal")
    insertion = add({35: 1, 7: -1}, boundary(wild5))
    check(not prime_image(insertion), "full insertion defect prime-balanced")
    check(insertion != K(35), "moving insertion label is not the full defect")
    source = Path(__file__).with_name('collatz_receipt_joins_addendum_20261005.py').read_text()
    node = next(n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'classes' for t in n.targets))
    classes = ast.literal_eval(node.value)
    check(len(classes) == 27 and sum(not row[2] for row in classes) == 8,
          "eight plain and nineteen multiplier classes")
    for n in range(1, 4096, 2):
        count = sum(n % modulus == residue for residue, modulus, _ in classes)
        check(count == (0 if n == 4095 else 1), "exact disjoint prefix-code coverage")
    check(799 % 64 == 31 and 1024 % 64 == 0, "J lies in multiplier11 cell")
    print(json.dumps({"status": "PASS; targeted type and scope controls only", "checks": checks,
                      "prime_child_error": error, "insertion5_at7": insertion,
                      "depth_two_witness": {"a": 5, "b": 17, "alpha": 8, "A": 6},
                      "semigroup_classes": {"total": 27, "plain": 8, "multiplier": 19}},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
