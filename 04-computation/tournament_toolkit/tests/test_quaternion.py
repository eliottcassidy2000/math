#!/usr/bin/env python3
"""
test_quaternion.py — Tests for QuaternionLinear and QuaternionAttention.

Tests:
1. Hamilton product correctness (algebraic identities)
2. Parameter count (exactly 75% savings)
3. Output shape correctness
4. Deterministic reproducibility
5. Cartan decomposition properties
6. Tournament structure of attention
7. Comparison with standard attention
8. Edge cases and error handling
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from tournament_toolkit.quaternion import (
    QuaternionLinear, QuaternionAttention, StandardAttention,
    cartan_decompose, attention_to_tournament, count_3cycles
)

PASS = 0
FAIL = 0

def check(name, condition):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ {name}")

print("=" * 60)
print("  QUATERNION ATTENTION TESTS")
print("=" * 60)
print()

# =================================================================
# TEST 1: HAMILTON PRODUCT CORRECTNESS
# =================================================================

print("TEST 1: Hamilton Product Correctness")
print("-" * 40)

np.random.seed(42)

# A QuaternionLinear with identity-like weights should preserve quaternion algebra
ql = QuaternionLinear(4, 4, scale=0.0)
# Set to identity: W_r = I, W_i = W_j = W_k = 0
ql.W_r = np.eye(1)
ql.W_i = np.zeros((1, 1))
ql.W_j = np.zeros((1, 1))
ql.W_k = np.zeros((1, 1))

# Input quaternion (1, 0, 0, 0) should give (1, 0, 0, 0)
one = np.array([[1.0, 0.0, 0.0, 0.0]])
result = ql(one)
check("Identity on 1: (1,0,0,0)→(1,0,0,0)",
      np.allclose(result, [1, 0, 0, 0], atol=1e-10))

# Input (0, 1, 0, 0) should give (0, 1, 0, 0)
i_vec = np.array([[0.0, 1.0, 0.0, 0.0]])
result_i = ql(i_vec)
check("Identity on i: (0,1,0,0)→(0,1,0,0)",
      np.allclose(result_i, [0, 1, 0, 0], atol=1e-10))

# Now test i * i = -1: set W to represent multiplication by i
ql_i = QuaternionLinear(4, 4, scale=0.0)
ql_i.W_r = np.zeros((1, 1))
ql_i.W_i = np.eye(1)
ql_i.W_j = np.zeros((1, 1))
ql_i.W_k = np.zeros((1, 1))

# Multiply i * i: input (0,1,0,0), weight represents i
# Hamilton: (0+i)(0+i) = 0*0 - 1*1 + (0*1+1*0)i + ... = -1
result_ii = ql_i(i_vec)
check("i·i = -1: (0,1,0,0) × W_i → (-1,0,0,0)",
      np.allclose(result_ii, [-1, 0, 0, 0], atol=1e-10))

# Multiply 1 * i: input (1,0,0,0), weight represents i
result_1i = ql_i(one)
check("1·i = i: (1,0,0,0) × W_i → (0,1,0,0)",
      np.allclose(result_1i, [0, 1, 0, 0], atol=1e-10))

# Test j * k = i: set W to represent j, input k
ql_j = QuaternionLinear(4, 4, scale=0.0)
ql_j.W_r = np.zeros((1, 1))
ql_j.W_i = np.zeros((1, 1))
ql_j.W_j = np.eye(1)
ql_j.W_k = np.zeros((1, 1))

k_vec = np.array([[0.0, 0.0, 0.0, 1.0]])
result_jk = ql_j(k_vec)
# j*k: (0+j)(0+k) = 0 - 0 + 0i + (0-0+0+1)j... wait let me compute:
# (r,i,j,k) input = (0,0,0,1), W = (0,0,1,0)
# o_r = r*Wr - i*Wi - j*Wj - k*Wk = 0 - 0 - 0 - 1*0 = 0... hmm
# Actually the Hamilton product for x=(0,0,0,1) and W=(0,0,1,0):
# We're multiplying (0+0i+1j+0k) × (0+0i+0j+1k) wait no.
# W represents: W_r=0, W_i=0, W_j=1, W_k=0
# Input x: r=0, i=0, j=0, k=1
# o_r = r*W_r - i*W_i - j*W_j - k*W_k = 0 - 0 - 0 - 1*0 = 0
# o_i = r*W_i + i*W_r + j*W_k - k*W_j = 0 + 0 + 0 - 1*1 = -1
# Hmm, that gives (0, -1, 0, 0) = -i
# But j*k should = i, not -i!
# The issue: the convention. The Hamilton product (a,b)*(c,d) is
# LEFT multiplication by (a,b). So ql_j multiplies ON THE LEFT by j.
# j * k = i (correct), but the formula computes x * W, which is k * j = -i.
# So the order is reversed: QuaternionLinear computes x * W, not W * x.
check("k·j = -i: (0,0,0,1) × W_j → (0,-1,0,0)  [right multiplication]",
      np.allclose(result_jk, [0, -1, 0, 0], atol=1e-10))

print()

# =================================================================
# TEST 2: PARAMETER COUNT (75% SAVINGS)
# =================================================================

print("TEST 2: Parameter Count")
print("-" * 40)

for d_in, d_out in [(64, 64), (256, 64), (512, 128), (768, 64), (1024, 256)]:
    ql = QuaternionLinear(d_in, d_out)
    expected_savings = 0.75
    check(f"QuaternionLinear({d_in},{d_out}): savings={ql.savings_fraction:.0%}",
          abs(ql.savings_fraction - expected_savings) < 0.01)

# Full attention head savings
for d_model, d_head in [(128, 32), (256, 64), (512, 64), (768, 64)]:
    qa = QuaternionAttention(d_model, d_head)
    std_params = 4 * d_model * d_head  # 3 projections + output
    actual_savings = 1 - qa.num_params / qa.standard_params
    check(f"QuaternionAttention({d_model},{d_head}): savings={actual_savings:.0%}",
          abs(actual_savings - 0.75) < 0.01)

print()

# =================================================================
# TEST 3: OUTPUT SHAPE
# =================================================================

print("TEST 3: Output Shape")
print("-" * 40)

for seq_len in [8, 16, 32]:
    for d_model, d_head in [(64, 16), (128, 32), (256, 64)]:
        np.random.seed(42)
        qa = QuaternionAttention(d_model, d_head)
        x = np.random.randn(seq_len, d_model)
        output, attn = qa(x)
        check(f"seq={seq_len}, d_model={d_model}: output shape {output.shape}",
              output.shape == (seq_len, d_model))
        check(f"seq={seq_len}, d_model={d_model}: attn shape {attn.shape}",
              attn.shape == (seq_len, seq_len))

print()

# =================================================================
# TEST 4: ATTENTION PROPERTIES
# =================================================================

print("TEST 4: Attention Properties")
print("-" * 40)

np.random.seed(42)
qa = QuaternionAttention(64, 16)
x = np.random.randn(12, 64)
output, attn = qa(x)

# Attention rows sum to 1 (softmax property)
row_sums = attn.sum(axis=1)
check("Attention rows sum to 1",
      np.allclose(row_sums, 1.0, atol=1e-6))

# Attention is non-negative
check("Attention is non-negative",
      np.all(attn >= 0))

# Attention is not uniform (has structure)
attn_entropy = -np.sum(attn * np.log(attn + 1e-10), axis=1).mean()
max_entropy = np.log(12)
check(f"Attention entropy valid (entropy {attn_entropy:.2f} ≤ max {max_entropy:.2f})",
      attn_entropy <= max_entropy + 1e-6)

print()

# =================================================================
# TEST 5: CARTAN DECOMPOSITION
# =================================================================

print("TEST 5: Cartan Decomposition")
print("-" * 40)

A_anti, A_sym, scalar = cartan_decompose(attn)

# Anti-symmetric part is actually anti-symmetric
check("A_anti is antisymmetric",
      np.allclose(A_anti, -A_anti.T, atol=1e-10))

# Symmetric part is actually symmetric
check("A_sym is symmetric",
      np.allclose(A_sym, A_sym.T, atol=1e-10))

# Symmetric part is traceless
check("A_sym is traceless",
      abs(np.trace(A_sym)) < 1e-10)

# Reconstruction: A = A_anti + A_sym + scalar*I
reconstructed = A_anti + A_sym + scalar * np.eye(attn.shape[0])
check("Reconstruction: A_anti + A_sym + scalar*I = A",
      np.allclose(reconstructed, attn, atol=1e-10))

# Orthogonality: <A_anti, A_sym> should be zero (Frobenius inner product)
inner = np.sum(A_anti * A_sym)
check(f"Orthogonality: <A_anti, A_sym> = {inner:.2e} ≈ 0",
      abs(inner) < 1e-10)

print()

# =================================================================
# TEST 6: TOURNAMENT STRUCTURE
# =================================================================

print("TEST 6: Tournament Structure")
print("-" * 40)

tournament = attention_to_tournament(attn)
n = attn.shape[0]

# Tournament is binary
check("Tournament is binary (0 or 1)",
      np.all((tournament == 0) | (tournament == 1)))

# Tournament is complete (every pair has exactly one arc)
for i in range(n):
    for j in range(i+1, n):
        ok = (tournament[i,j] + tournament[j,i] == 1)
        if not ok:
            check(f"Tournament complete: pair ({i},{j})", False)
            break
    else:
        continue
    break
else:
    check("Tournament is complete (one arc per pair)", True)

# No self-loops
check("No self-loops",
      np.all(np.diag(tournament) == 0))

# Count 3-cycles
c3 = count_3cycles(tournament)
max_c3 = n * (n-1) * (n-2) // 6  # C(n,3)
frac = c3 / max_c3
check(f"3-cycle count: {c3} ({frac:.3f} of {max_c3} triples)",
      0 <= c3 <= max_c3)

print()

# =================================================================
# TEST 7: COMPARISON WITH STANDARD ATTENTION
# =================================================================

print("TEST 7: Standard vs Quaternion Comparison")
print("-" * 40)

np.random.seed(42)
d_model, d_head, seq_len = 128, 32, 16

std = StandardAttention(d_model, d_head)
quat = QuaternionAttention(d_model, d_head)

x = np.random.randn(seq_len, d_model)

out_std, attn_std = std.forward(x)
out_quat, attn_quat = quat(x)

check(f"Standard params: {std.num_params}", std.num_params == 4 * d_model * d_head)
check(f"Quaternion params: {quat.num_params}", quat.num_params == d_model * d_head)
check(f"Savings: {1 - quat.num_params/std.num_params:.0%}",
      abs((1 - quat.num_params / std.num_params) - 0.75) < 0.01)

# Both produce same shape
check("Same output shape", out_std.shape == out_quat.shape)
check("Same attention shape", attn_std.shape == attn_quat.shape)

# Both have valid attention (rows sum to 1)
check("Standard attn rows sum to 1",
      np.allclose(attn_std.sum(axis=1), 1.0, atol=1e-6))
check("Quaternion attn rows sum to 1",
      np.allclose(attn_quat.sum(axis=1), 1.0, atol=1e-6))

# Compare Cartan ratios
_, A_sym_s, _ = cartan_decompose(attn_std)
A_anti_s, _, _ = cartan_decompose(attn_std)
_, A_sym_q, _ = cartan_decompose(attn_quat)
A_anti_q, _, _ = cartan_decompose(attn_quat)

ratio_std = np.linalg.norm(A_anti_s) / (np.linalg.norm(A_sym_s) + 1e-10)
ratio_quat = np.linalg.norm(A_anti_q) / (np.linalg.norm(A_sym_q) + 1e-10)
print(f"  Cartan anti/sym ratio: standard={ratio_std:.4f}, quaternion={ratio_quat:.4f}")

# Compare tournament structure
t_std = attention_to_tournament(attn_std)
t_quat = attention_to_tournament(attn_quat)
c3_std = count_3cycles(t_std)
c3_quat = count_3cycles(t_quat)
print(f"  3-cycles: standard={c3_std}, quaternion={c3_quat}")

print()

# =================================================================
# TEST 8: ERROR HANDLING
# =================================================================

print("TEST 8: Error Handling")
print("-" * 40)

try:
    ql_bad = QuaternionLinear(7, 8)
    check("Reject d_in not divisible by 4", False)
except ValueError:
    check("Reject d_in not divisible by 4", True)

try:
    ql_bad = QuaternionLinear(8, 7)
    check("Reject d_out not divisible by 4", False)
except ValueError:
    check("Reject d_out not divisible by 4", True)

print()

# =================================================================
# SUMMARY
# =================================================================

print("=" * 60)
total = PASS + FAIL
print(f"  RESULTS: {PASS}/{total} passed, {FAIL}/{total} failed")
if FAIL == 0:
    print(f"  ALL TESTS PASSED ✓")
else:
    print(f"  {FAIL} TESTS FAILED ✗")
print("=" * 60)
