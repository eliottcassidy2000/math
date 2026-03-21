#!/usr/bin/env python3
"""
test_shadow_cert.py — Rigorous tests for ShadowCert
opus-2026-03-21-S103

Tests:
  1. Basic functionality on valid tournaments
  2. Edge case: n=2 (degenerate)
  3. Edge case: n=1 (trivial)
  4. Incomplete data (not all pairs compared)
  5. Adversarial: fabricated scores that violate Landau condition
  6. Adversarial: scores that sum to wrong total
  7. Known tournaments at n=5 (verify against exact H)
  8. Large n (n=20, estimation only)
  9. Stress test: all 1024 tournaments at n=5
"""

import sys
sys.path.insert(0, '..')
from tournament_toolkit.shadow_cert import ShadowCert
from itertools import permutations, combinations
import numpy as np

passed = 0
failed = 0

def check(name, condition):
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✓ {name}")
    else:
        failed += 1
        print(f"  ✗ FAIL: {name}")

print("=" * 60)
print("  TEST SUITE: ShadowCert")
print("=" * 60)
print()

# Test 1: Basic valid tournament (n=5 transitive)
print("TEST 1: Basic valid tournament (n=5 transitive)")
cert = ShadowCert(5, [4, 3, 2, 1, 0])
check("n=5 transitive: sum correct", cert._integrity['sum_correct'])
check("n=5 transitive: Landau valid", cert._integrity['landau_valid'])
check("n=5 transitive: H is odd", cert._integrity['H_is_odd'])
check("n=5 transitive: verdict DECISIVE or DOMINATED", cert.verdict in ('DECISIVE', 'DOMINATED'))
check("n=5 transitive: verify passes", cert.verify())
print()

# Test 2: n=5 regular (cyclic)
print("TEST 2: n=5 regular tournament")
cert2 = ShadowCert(5, [2, 2, 2, 2, 2])
check("n=5 regular: sum correct", cert2._integrity['sum_correct'])
check("n=5 regular: regularity = 100%", abs(cert2.regularity - 1.0) < 0.01)
check("n=5 regular: S₂ = 0", cert2.S2 == 0)
check("n=5 regular: verify passes", cert2.verify())
print()

# Test 3: Edge case n=2
print("TEST 3: Edge case n=2")
cert3 = ShadowCert(2, [1, 0])
check("n=2: verify passes", cert3.verify())
check("n=2: H estimate = 1", cert3.H_estimate == 1)
print()

# Test 4: Edge case n=1
print("TEST 4: Edge case n=1")
cert4 = ShadowCert(1, [0])
check("n=1: H estimate = 1", cert4.H_estimate == 1)
print()

# Test 5: Incomplete data
print("TEST 5: Incomplete data (only 3 of 6 pairs compared)")
cert5 = ShadowCert.from_comparisons([("A", "B"), ("B", "C"), ("C", "D")])
check("Incomplete: n=4", cert5.n == 4)
check("Incomplete: completeness < 1", cert5.completeness < 1.0)
check("Incomplete: completeness = 0.5", abs(cert5.completeness - 0.5) < 0.01)
check("Incomplete: not marked as complete", not cert5._integrity.get('is_complete', False))
print()

# Test 6: Adversarial — scores don't sum to C(n,2)
print("TEST 6: Adversarial — wrong score sum")
cert6 = ShadowCert(5, [4, 4, 4, 4, 4])  # Sum = 20 ≠ C(5,2)=10
check("Wrong sum: detected", not cert6._integrity['sum_correct'])
check("Wrong sum: verify fails", not cert6.verify())
print()

# Test 7: Adversarial — violates Landau condition
print("TEST 7: Adversarial — Landau violation")
cert7 = ShadowCert(5, [0, 0, 0, 4, 6])  # sum=10 but Landau fails
# Actually 6 > n-1=4, so range check fails first
check("Landau/range: detected", not cert7._integrity['range_valid'] or not cert7._integrity['landau_valid'])
check("Landau/range: verify fails", not cert7.verify())
print()

# Test 8: Known tournament at n=5, verify exact H
print("TEST 8: Known n=5 tournament — exact H verification")
# Build transitive tournament adjacency
adj_trans = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        adj_trans[i][j] = 1
cert8 = ShadowCert.from_adjacency(adj_trans)
H_exact = cert8.compute_exact_H(adj_trans)
check("Transitive H = 1", H_exact == 1)
check("Scores = (4,3,2,1,0)", cert8.scores == [4, 3, 2, 1, 0])

# Build cyclic C_5
adj_cyc = [[0]*5 for _ in range(5)]
for i in range(5):
    adj_cyc[i][(i+1) % 5] = 1
    adj_cyc[i][(i+2) % 5] = 1
cert8b = ShadowCert.from_adjacency(adj_cyc)
H_exact_cyc = cert8b.compute_exact_H(adj_cyc)
check("Cyclic C_5 H = 15", H_exact_cyc == 15)
check("Cyclic scores = (2,2,2,2,2)", cert8b.scores == [2, 2, 2, 2, 2])
print()

# Test 9: Large n (estimation only)
print("TEST 9: Large n = 20 (estimation only)")
# Create a near-transitive score sequence
scores_20 = list(range(19, -1, -1))  # 19, 18, ..., 0
cert9 = ShadowCert(20, scores_20)
check("n=20: verify passes", cert9.verify())
check("n=20: S₂ > 0", cert9.S2 > 0)
check("n=20: regularity < 50%", cert9.regularity < 0.5)
check("n=20: H estimate is odd", cert9.H_estimate % 2 == 1)
print()

# Test 10: Stress test — all n=5 tournaments
print("TEST 10: Stress test — all 1024 n=5 tournaments")
n_valid = 0
n_invalid = 0
for mask in range(1024):
    adj = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (mask >> idx) & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    cert = ShadowCert.from_adjacency(adj)
    if cert.verify():
        n_valid += 1
    else:
        n_invalid += 1

check("All 1024 n=5 tournaments pass verification", n_valid == 1024 and n_invalid == 0)
print()

# Test 11: Certificate report format
print("TEST 11: Certificate report generation")
cert11 = ShadowCert.from_comparisons([
    ("MIT", "Stanford"), ("MIT", "Harvard"), ("Stanford", "Harvard"),
    ("MIT", "Caltech"), ("Stanford", "Caltech"), ("Harvard", "Caltech"),
])
report = cert11.report()
check("Report contains VERDICT", "VERDICT" in report)
check("Report contains S₂", "S₂" in report or "S2" in report or "irregularity" in report.lower())
check("Report contains INTEGRITY", "INTEGRITY" in report)
print()
print(cert11.report())
print()

# Summary
print("=" * 60)
print(f"  RESULTS: {passed} passed, {failed} failed out of {passed + failed}")
print("=" * 60)
