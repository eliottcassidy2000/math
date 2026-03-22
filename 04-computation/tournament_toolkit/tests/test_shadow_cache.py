#!/usr/bin/env python3
"""
test_shadow_cache.py — Tests for shadow-based KV cache compression.

Tests:
1. Basic functionality (add, attend, shapes)
2. Compression ratio (must beat standard)
3. Attention quality (cosine similarity to standard)
4. Scaling behavior (quality vs compression at different seq lengths)
5. Importance tracking (recent and high-attention tokens preserved)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from tournament_toolkit.shadow_cache import ShadowCache, StandardCache

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
print("  SHADOW CACHE COMPRESSION TESTS")
print("=" * 60)
print()

np.random.seed(42)

# =================================================================
# TEST 1: BASIC FUNCTIONALITY
# =================================================================

print("TEST 1: Basic Functionality")
print("-" * 40)

d_model = 128
d_head = 64
n_heads = 2
max_seq = 1024

cache = ShadowCache(d_model, n_heads, max_seq)
std_cache = StandardCache(d_head, max_seq)

# Add some tokens
for t in range(100):
    k = np.random.randn(d_head)
    v = np.random.randn(d_head)
    cache.add_token(k, v)
    std_cache.add_token(k, v)

check("Cache length = 100", cache.length == 100)
check("Standard cache length = 100", std_cache.length == 100)

# Attend
query = np.random.randn(d_head)
out_shadow, attn_shadow = cache.attend(query)
out_std, attn_std = std_cache.attend(query)

check("Shadow output shape correct", out_shadow.shape == (d_head,))
check("Shadow attention shape correct", attn_shadow.shape == (100,))
check("Standard output shape correct", out_std.shape == (d_head,))
check("Attention sums to 1 (shadow)",
      abs(attn_shadow.sum() - 1.0) < 1e-6)
check("Attention sums to 1 (standard)",
      abs(attn_std.sum() - 1.0) < 1e-6)

print()

# =================================================================
# TEST 2: COMPRESSION RATIO
# =================================================================

print("TEST 2: Compression Ratio")
print("-" * 40)

stats = cache.stats()
print(f"  Shadow cache: {stats['memory_MB']:.2f} MB")
print(f"  Standard cache: {stats['standard_memory_MB']:.2f} MB")
print(f"  Compression ratio: {stats['compression_ratio']:.2f}")
print(f"  Savings: {stats['savings_percent']:.1f}%")

check("Shadow uses less memory than standard",
      stats['memory_MB'] < stats['standard_memory_MB'])
check("Savings > 30%", stats['savings_percent'] > 30)

print()

# =================================================================
# TEST 3: ATTENTION QUALITY
# =================================================================

print("TEST 3: Attention Quality")
print("-" * 40)

# Run many queries and measure cosine similarity of outputs
np.random.seed(123)
cache2 = ShadowCache(d_model, n_heads, max_seq, shadow_dim=d_head // 2)
std2 = StandardCache(d_head, max_seq)

n_tokens = 200
for t in range(n_tokens):
    k = np.random.randn(d_head) * 0.1
    v = np.random.randn(d_head) * 0.1
    cache2.add_token(k, v)
    std2.add_token(k, v)

cosine_sims = []
attn_correlations = []
n_queries = 50

for q in range(n_queries):
    query = np.random.randn(d_head) * 0.1
    out_s, attn_s = cache2.attend(query)
    out_r, attn_r = std2.attend(query)

    # Cosine similarity of outputs
    norm_s = np.linalg.norm(out_s)
    norm_r = np.linalg.norm(out_r)
    if norm_s > 1e-10 and norm_r > 1e-10:
        cos_sim = np.dot(out_s, out_r) / (norm_s * norm_r)
        cosine_sims.append(cos_sim)

    # Correlation of attention weights
    if len(attn_s) == len(attn_r):
        corr = np.corrcoef(attn_s, attn_r)[0, 1]
        if not np.isnan(corr):
            attn_correlations.append(corr)

mean_cos = np.mean(cosine_sims) if cosine_sims else 0
mean_corr = np.mean(attn_correlations) if attn_correlations else 0

print(f"  Mean output cosine similarity: {mean_cos:.4f}")
print(f"  Mean attention correlation: {mean_corr:.4f}")

check(f"Output cosine similarity > 0.5", mean_cos > 0.5)
check(f"Attention correlation > 0.5", mean_corr > 0.5)

print()

# =================================================================
# TEST 4: SCALING BEHAVIOR
# =================================================================

print("TEST 4: Scaling Behavior")
print("-" * 40)

print(f"  {'seq_len':>8s}  {'shadow_MB':>10s}  {'std_MB':>8s}  "
      f"{'ratio':>6s}  {'savings':>8s}  {'cos_sim':>8s}")
print(f"  {'─'*8}  {'─'*10}  {'─'*8}  {'─'*6}  {'─'*8}  {'─'*8}")

for seq_len in [50, 100, 200, 500, 1000]:
    np.random.seed(42)
    sc = ShadowCache(d_model, n_heads, max_seq=seq_len+10,
                     shadow_dim=d_head // 4, top_k_fraction=0.05)
    rc = StandardCache(d_head, max_seq=seq_len+10)

    for t in range(seq_len):
        k = np.random.randn(d_head) * 0.1
        v = np.random.randn(d_head) * 0.1
        sc.add_token(k, v)
        rc.add_token(k, v)

    # Measure quality
    sims = []
    for _ in range(20):
        q = np.random.randn(d_head) * 0.1
        o_s, _ = sc.attend(q)
        o_r, _ = rc.attend(q)
        n_s, n_r = np.linalg.norm(o_s), np.linalg.norm(o_r)
        if n_s > 1e-10 and n_r > 1e-10:
            sims.append(np.dot(o_s, o_r) / (n_s * n_r))

    st = sc.stats()
    mean_sim = np.mean(sims) if sims else 0
    print(f"  {seq_len:8d}  {st['memory_MB']:10.3f}  "
          f"{st['standard_memory_MB']:8.3f}  {st['compression_ratio']:6.2f}  "
          f"{st['savings_percent']:7.1f}%  {mean_sim:8.4f}")

print()

# =================================================================
# TEST 5: IMPORTANCE TRACKING
# =================================================================

print("TEST 5: Importance Tracking")
print("-" * 40)

np.random.seed(42)
cache3 = ShadowCache(d_model, n_heads, max_seq=500, top_k_fraction=0.1)

# Add tokens where some have distinctive keys
for t in range(200):
    if t == 50 or t == 100 or t == 150:
        # "Important" tokens: distinctive key direction
        k = np.ones(d_head) * 5.0  # very different from random
        v = np.ones(d_head) * 3.0
    else:
        k = np.random.randn(d_head) * 0.1
        v = np.random.randn(d_head) * 0.1
    cache3.add_token(k, v)

# Query aligned with distinctive tokens
query = np.ones(d_head)
_, attn = cache3.attend(query)

# The distinctive tokens should get high attention
check("Distinctive token 50 gets above-average attention",
      attn[50] > 1.5 / 200)
check("Distinctive token 100 gets above-average attention",
      attn[100] > 1.5 / 200)
check("Distinctive token 150 gets above-average attention",
      attn[150] > 1.5 / 200)

# Recent tokens should be in full_keys
check("Recent tokens in full precision",
      199 in cache3.full_keys)

print()

# =================================================================
# TEST 6: THE CODING THEORY PERSPECTIVE
# =================================================================

print("TEST 6: Coding Theory Metrics")
print("-" * 40)

np.random.seed(42)
d_head_ct = 32
cache_ct = ShadowCache(d_model=64, n_heads=2, max_seq=200, shadow_dim=16)
std_ct = StandardCache(d_head_ct, max_seq=200)

for t in range(100):
    k = np.random.randn(d_head_ct) * 0.1
    v = np.random.randn(d_head_ct) * 0.1
    cache_ct.add_token(k, v)
    std_ct.add_token(k, v)

# "Syndrome" = the shadow scores (compressed representation)
# "Codeword" = full KV vectors
# "Code rate" = shadow_dim / d_head

code_rate = cache_ct.shadow_dim / d_head_ct
print(f"  Code rate (shadow_dim / d_head): {code_rate:.3f}")
print(f"  = {cache_ct.shadow_dim} / {d_head_ct}")
print(f"  This means: each key is compressed to {code_rate:.0%} of original")

# "Decoding success" = how well we reconstruct the output
queries = [np.random.randn(d_head_ct) * 0.1 for _ in range(30)]
successes = 0
for q in queries:
    o_s, _ = cache_ct.attend(q)
    o_r, _ = std_ct.attend(q)
    n_s, n_r = np.linalg.norm(o_s), np.linalg.norm(o_r)
    if n_s > 1e-10 and n_r > 1e-10:
        cos = np.dot(o_s, o_r) / (n_s * n_r)
        if cos > 0.3:
            successes += 1

decoding_rate = successes / len(queries)
print(f"  Decoding success rate (cos > 0.3): {decoding_rate:.0%}")
print(f"  = fraction of queries with reasonable reconstruction")
print(f"  Note: with random data and random projection, quality is limited.")
print(f"  With LEARNED projections (trained), this would be much higher.")

check(f"Decoding rate > 50% (cos > 0.3 threshold)", decoding_rate > 0.5)

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
