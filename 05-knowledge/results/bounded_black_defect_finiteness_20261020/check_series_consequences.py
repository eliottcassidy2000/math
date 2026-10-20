#!/usr/bin/env python3
"""Exact checks of the rational-series corollaries, not a center census.

D=z^4 is the known center control. A second hypothetical finite D includes
both even and odd unknown-center terms to test that no parity/classification
assumption is hidden in the recurrence or asymptotic error estimate.
"""
import json
from pathlib import Path

N = 1000
a = [0] * (N + 1)
a[0] = 1
for n in range(1, N + 1):
    a[n] = a[n - 1] + (a[n - 3] if n >= 3 else 0)
for n in range(N + 1):
    assert a[n] >= a[n // 2] * a[(n + 1) // 2]
    if n:
        assert a[n] >= a[n - 1]
    if n % 3 == 0:
        assert a[n] >= 2 ** (n // 3)

def coefficient(values, n):
    return values[n] if n >= 0 else 0

def half_coeff(n):
    return a[n // 2] if n >= 0 and n % 2 == 0 else 0

cases = []
for name, D in (("known_S4_only_control", {4: 1}),
                ("hypothetical_finite_alphabet_not_a_claim", {4: 1, 10: 2, 15: 1})):
    mixed = [sum(c * half_coeff(n - s) for s, c in D.items())
             for n in range(N + 1)]
    palindrome = [sum(half_coeff(n - s) for s in (0, 1, 3))
                  for n in range(N + 1)]
    assert all((a[n] - palindrome[n]) % 2 == 0 for n in range(N + 1))
    leaf = [(a[n] - palindrome[n]) // 2 + mixed[n] for n in range(N + 1)]
    assert all(x >= 0 for x in leaf)
    two_F = {s: 2 * c for s, c in D.items()}
    for s in (0, 1, 3):
        two_F[s] = two_F.get(s, 0) - 1
    two_C = sum(abs(c) for c in two_F.values())
    for n in range(N + 1):
        assert abs(2 * leaf[n] - a[n]) <= two_C * a[n // 2]
        if n > max(D):
            assert mixed[n] == coefficient(mixed, n - 2) + coefficient(mixed, n - 6)
        if n > max(max(D), 3) + 3:
            assert leaf[n] == (coefficient(leaf, n - 1) + coefficient(leaf, n - 2) - coefficient(leaf, n - 5)
                               + coefficient(leaf, n - 6) - coefficient(leaf, n - 7) - coefficient(leaf, n - 9))
    cases.append({"name": name, "D_coefficients": D, "checked_through": N,
                  "initial_leaf_counts_n3_through_n9": leaf[3:10],
                  "twice_absolute_error_constant": two_C})
result = {"status": "PASS", "scope": "Exact formal-series arithmetic, not center classification",
          "composition_checks_through": N, "cases": cases}
print(json.dumps(result, indent=2, sort_keys=True))
