#!/usr/bin/env python3
"""
gpt2_cartan_probe_s95.py — opus-2026-03-21-S95

THE CARTAN TOURNAMENT PROBE ON ACTUAL GPT-2 ATTENTION

Key experiment: Extract attention matrices from GPT-2, decompose them via
the Cartan decomposition gl(n) = so(n) + p + R*I, and measure:

1. How much energy is in the TOURNAMENT (antisymmetric) sector vs
   COOPERATION (symmetric traceless) sector vs SCALAR sector?
2. How does this compare to RANDOM attention?
3. Does the balance vary across layers? (Early vs middle vs late)
4. Does the balance vary across attention heads?
5. The COMMUTATOR ||[A_anti, A_sym]||: how entangled are the sectors?
6. Can we threshold attention to get actual tournaments and compute H(T)?

This is INV-180 step 1: the first real test of the Cartan bridge on LLM data.

Author: opus-2026-03-21-S95
"""

import sys
import numpy as np
import torch
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  GPT-2 CARTAN TOURNAMENT PROBE")
print("  opus-2026-03-21-S95")
print("=" * 72)

# ========================================================================
# PART 0: Load GPT-2 and extract attention matrices
# ========================================================================
print("\n--- Loading GPT-2 small ---")

from transformers import GPT2LMHeadModel, GPT2Tokenizer

tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
model = GPT2LMHeadModel.from_pretrained("gpt2", output_attentions=True)
model.eval()

# Test sentences with varying structure
sentences = [
    "The cat sat on the mat and looked at the dog.",
    "Mathematical proofs require careful reasoning and attention to detail.",
    "One plus one equals two, and two plus two equals four.",
    "The quick brown fox jumps over the lazy dog near the river.",
    "In tournament theory, every pair of vertices has exactly one directed edge.",
]

print(f"  Model: GPT-2 small (12 layers, 12 heads, d_model=768)")
print(f"  Testing on {len(sentences)} sentences\n")

def cartan_decompose(M):
    """Cartan decomposition of n x n matrix.
    Returns (anti_frac, sym_frac, scalar_frac, sigma, comm_norm2)
    """
    n = M.shape[0]
    anti = (M - M.T) / 2.0
    sym = (M + M.T) / 2.0
    sigma = np.trace(M) / n
    scalar = sigma * np.eye(n)
    sym_tl = sym - scalar

    norm2 = np.sum(M**2)
    if norm2 < 1e-15:
        return 0, 0, 0, sigma, 0

    anti_frac = np.sum(anti**2) / norm2
    sym_frac = np.sum(sym_tl**2) / norm2
    scalar_frac = np.sum(scalar**2) / norm2

    comm = anti @ sym_tl - sym_tl @ anti
    comm_norm2 = np.sum(comm**2)

    return anti_frac, sym_frac, scalar_frac, sigma, comm_norm2

def count_hp_small(T, n):
    """Count Hamiltonian paths for small n via DP."""
    if n > 12:
        return -1  # too slow
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if T[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

# ========================================================================
# PART 1: Extract and analyze attention patterns
# ========================================================================
print("=" * 72)
print("  PART 1: ATTENTION CARTAN DECOMPOSITION BY LAYER AND HEAD")
print("=" * 72)

all_results = []

for si, sentence in enumerate(sentences):
    inputs = tokenizer(sentence, return_tensors="pt")
    n_tokens = inputs['input_ids'].shape[1]
    print(f"\nSentence {si}: \"{sentence[:50]}...\" ({n_tokens} tokens)")

    with torch.no_grad():
        outputs = model(**inputs)
        attentions = outputs.attentions  # tuple of (batch, heads, seq, seq)

    # Analyze each layer and head
    for layer_idx in range(12):
        attn = attentions[layer_idx][0].numpy()  # (heads, seq, seq)
        for head_idx in range(12):
            A = attn[head_idx]  # n_tokens x n_tokens attention matrix
            n = A.shape[0]

            af, sf, scf, sigma, comm2 = cartan_decompose(A)

            # Tournament from attention (threshold)
            T = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i+1, n):
                    if A[i,j] > A[j,i]:
                        T[i][j] = 1
                    else:
                        T[j][i] = 1

            # Score sequence
            scores = T.sum(axis=1)
            S2 = np.sum((scores - (n-1)/2)**2)

            # H for small sequences
            H = count_hp_small(T, n) if n <= 12 else -1

            all_results.append({
                'sentence': si,
                'layer': layer_idx,
                'head': head_idx,
                'n': n,
                'anti_frac': af,
                'sym_frac': sf,
                'scalar_frac': scf,
                'sigma': sigma,
                'comm_norm2': comm2,
                'H': H,
                'S2': S2,
                'scores': tuple(sorted(scores)),
                'is_regular': (S2 < 0.01),
            })

# ========================================================================
# PART 2: Summary statistics
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: GLOBAL STATISTICS")
print("=" * 72)

anti_fracs = [r['anti_frac'] for r in all_results]
sym_fracs = [r['sym_frac'] for r in all_results]
scalar_fracs = [r['scalar_frac'] for r in all_results]
comms = [r['comm_norm2'] for r in all_results]
Hs = [r['H'] for r in all_results if r['H'] > 0]

print(f"\n  Total attention matrices analyzed: {len(all_results)}")
print(f"\n  Cartan fractions (mean ± std):")
print(f"    Tournament (antisymmetric): {np.mean(anti_fracs):.4f} ± {np.std(anti_fracs):.4f}")
print(f"    Cooperation (sym traceless): {np.mean(sym_fracs):.4f} ± {np.std(sym_fracs):.4f}")
print(f"    Scalar:                      {np.mean(scalar_fracs):.4f} ± {np.std(scalar_fracs):.4f}")
print(f"\n  Commutator ||[A,S]||² (mean): {np.mean(comms):.6f}")
print(f"  H (tournament HP count): mean={np.mean(Hs):.1f}, "
      f"min={min(Hs)}, max={max(Hs)}, median={np.median(Hs):.0f}")

# Compare with random baseline
print("\n  --- Random softmax attention baseline ---")
rng = np.random.RandomState(42)
n_test = all_results[0]['n']
rand_anti = []
for _ in range(1000):
    logits = rng.randn(n_test, n_test)
    exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
    A_rand = exp_l / exp_l.sum(axis=1, keepdims=True)
    af, sf, scf, _, _ = cartan_decompose(A_rand)
    rand_anti.append(af)

print(f"  Random anti_frac (n={n_test}): {np.mean(rand_anti):.4f} ± {np.std(rand_anti):.4f}")
print(f"  GPT-2  anti_frac:              {np.mean(anti_fracs):.4f} ± {np.std(anti_fracs):.4f}")

shift = np.mean(anti_fracs) - np.mean(rand_anti)
print(f"\n  SHIFT: GPT-2 - Random = {shift:+.4f}")
if shift > 0:
    print(f"  → Training INCREASES tournament (antisymmetric) energy!")
elif shift < 0:
    print(f"  → Training DECREASES tournament energy (more cooperation)")
else:
    print(f"  → No significant shift")

# ========================================================================
# PART 3: Layer-by-layer analysis
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: CARTAN FRACTIONS BY LAYER")
print("=" * 72)

print(f"\n  {'Layer':>5s} {'anti':>8s} {'sym':>8s} {'scalar':>8s} {'||[A,S]||²':>12s} {'H_mean':>8s}")

for layer in range(12):
    layer_results = [r for r in all_results if r['layer'] == layer]
    af = np.mean([r['anti_frac'] for r in layer_results])
    sf = np.mean([r['sym_frac'] for r in layer_results])
    scf = np.mean([r['scalar_frac'] for r in layer_results])
    cm = np.mean([r['comm_norm2'] for r in layer_results])
    Hs_layer = [r['H'] for r in layer_results if r['H'] > 0]
    H_mean = np.mean(Hs_layer) if Hs_layer else 0

    print(f"  {layer:5d} {af:8.4f} {sf:8.4f} {scf:8.4f} {cm:12.6f} {H_mean:8.1f}")

# ========================================================================
# PART 4: Head-by-head analysis
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: CARTAN FRACTIONS BY HEAD (averaged across layers)")
print("=" * 72)

print(f"\n  {'Head':>5s} {'anti':>8s} {'sym':>8s} {'scalar':>8s}")

for head in range(12):
    head_results = [r for r in all_results if r['head'] == head]
    af = np.mean([r['anti_frac'] for r in head_results])
    sf = np.mean([r['sym_frac'] for r in head_results])
    scf = np.mean([r['scalar_frac'] for r in head_results])
    print(f"  {head:5d} {af:8.4f} {sf:8.4f} {scf:8.4f}")

# ========================================================================
# PART 5: Tournament invariants of attention
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: TOURNAMENT INVARIANTS OF ATTENTION PATTERNS")
print("=" * 72)

print("\n  S₂ (Landau irregularity) by layer:")
print(f"  {'Layer':>5s} {'S₂_mean':>10s} {'S₂_std':>10s} {'reg_frac':>10s}")

for layer in range(12):
    layer_results = [r for r in all_results if r['layer'] == layer]
    S2s = [r['S2'] for r in layer_results]
    reg_frac = np.mean([r['is_regular'] for r in layer_results])
    print(f"  {layer:5d} {np.mean(S2s):10.2f} {np.std(S2s):10.2f} {reg_frac:10.3f}")

# ========================================================================
# PART 6: Score sequence analysis
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: HOW REGULAR ARE GPT-2 ATTENTION TOURNAMENTS?")
print("=" * 72)

# What fraction have near-regular score sequences?
n_test = all_results[0]['n']
all_S2 = [r['S2'] for r in all_results]
# Theoretical: transitive has S2 = n(n²-1)/12, regular has S2 = 0
S2_max = n_test * (n_test**2 - 1) / 12
print(f"\n  n = {n_test} tokens")
print(f"  S₂ range: [0, {S2_max:.1f}] (regular to transitive)")
print(f"  GPT-2 S₂ mean: {np.mean(all_S2):.2f}")
print(f"  GPT-2 S₂ median: {np.median(all_S2):.2f}")
print(f"  Random S₂ (uniform tournament): {n_test * (n_test - 1) / 4:.2f} (expected)")

# Random tournament S₂ for comparison
rand_S2 = []
rng = np.random.RandomState(42)
for _ in range(1000):
    logits = rng.randn(n_test, n_test)
    exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
    A_rand = exp_l / exp_l.sum(axis=1, keepdims=True)
    T_rand = np.zeros((n_test, n_test), dtype=int)
    for i in range(n_test):
        for j in range(i+1, n_test):
            if A_rand[i,j] > A_rand[j,i]:
                T_rand[i][j] = 1
            else:
                T_rand[j][i] = 1
    scores_rand = T_rand.sum(axis=1)
    rand_S2.append(np.sum((scores_rand - (n_test-1)/2)**2))

print(f"  Random S₂ (softmax attention): {np.mean(rand_S2):.2f} ± {np.std(rand_S2):.2f}")

# ========================================================================
# PART 7: H(T_attention) and Rédei check
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: H(T_attention) — HAMILTONIAN PATH COUNTS")
print("=" * 72)

Hs = [r['H'] for r in all_results if r['H'] > 0]
print(f"\n  Total attention tournaments with H computed: {len(Hs)}")
print(f"  H values: min={min(Hs)}, max={max(Hs)}, mean={np.mean(Hs):.1f}")
print(f"  H distribution:")

from collections import Counter
H_counts = Counter(Hs)
for H_val in sorted(H_counts.keys()):
    pct = 100 * H_counts[H_val] / len(Hs)
    print(f"    H={H_val:6d}: {H_counts[H_val]:4d} ({pct:5.1f}%)")

# Rédei check
all_odd = all(H % 2 == 1 for H in Hs)
print(f"\n  Rédei's theorem (H always odd): {'VERIFIED ✓' if all_odd else 'VIOLATED ✗'}")

# ========================================================================
# PART 8: Softmax temperature and tournament structure
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 8: ATTENTION SHARPNESS AND TOURNAMENT STRUCTURE")
print("=" * 72)
print("""
  The key question: does sharper attention (lower entropy) produce
  more tournament-like structure?

  We measure "sharpness" by the entropy of each attention row.
  Lower entropy = sharper = more tournament-like (one token dominates).
""")

for layer in range(12):
    layer_results = [r for r in all_results if r['layer'] == layer and r['sentence'] == 0]
    if not layer_results:
        continue

    # Get actual attention for this layer
    inputs = tokenizer(sentences[0], return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs)
        attentions = outputs.attentions

    attn = attentions[layer][0].numpy()  # (heads, seq, seq)

    # Compute entropy per head
    entropies = []
    anti_fracs_layer = []
    for head in range(12):
        A = attn[head]
        n = A.shape[0]
        # Row-wise entropy
        ent = -np.sum(A * np.log(A + 1e-12), axis=1)
        mean_ent = np.mean(ent)
        entropies.append(mean_ent)

        af, sf, scf, _, _ = cartan_decompose(A)
        anti_fracs_layer.append(af)

    if layer in [0, 5, 11]:
        corr = np.corrcoef(entropies, anti_fracs_layer)[0, 1]
        print(f"  Layer {layer:2d}: entropy↔anti_frac corr = {corr:+.4f}")

# ========================================================================
# PART 9: THE KEY FINDING — TRAINED vs RANDOM
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 9: THE KEY FINDING — TRAINED vs RANDOM CARTAN BALANCE")
print("=" * 72)

# Aggregate by layer
print(f"\n  {'Layer':>5s} {'GPT2_anti':>10s} {'Rand_anti':>10s} {'Shift':>8s} {'Direction':>15s}")

for layer in range(12):
    layer_results = [r for r in all_results if r['layer'] == layer]
    gpt2_anti = np.mean([r['anti_frac'] for r in layer_results])

    # Random baseline at same sequence length
    n_test = layer_results[0]['n']
    rng = np.random.RandomState(42 + layer)
    rand_anti_l = []
    for _ in range(200):
        logits = rng.randn(n_test, n_test)
        exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
        A_rand = exp_l / exp_l.sum(axis=1, keepdims=True)
        af, _, _, _, _ = cartan_decompose(A_rand)
        rand_anti_l.append(af)

    rand_anti = np.mean(rand_anti_l)
    shift = gpt2_anti - rand_anti
    direction = "MORE TOURN" if shift > 0.005 else ("MORE COOP" if shift < -0.005 else "~NEUTRAL")

    print(f"  {layer:5d} {gpt2_anti:10.4f} {rand_anti:10.4f} {shift:+8.4f} {direction:>15s}")

# ========================================================================
# PART 10: SYNTHESIS
# ========================================================================
print("\n\n" + "=" * 72)
print("  SYNTHESIS: THE CARTAN BRIDGE IN GPT-2")
print("=" * 72)

gpt2_mean_anti = np.mean(anti_fracs)
gpt2_mean_sym = np.mean(sym_fracs)
gpt2_mean_scalar = np.mean(scalar_fracs)

print(f"""
  GPT-2 ATTENTION CARTAN BALANCE:
    Tournament (antisymmetric): {gpt2_mean_anti:.4f}
    Cooperation (sym traceless): {gpt2_mean_sym:.4f}
    Scalar:                      {gpt2_mean_scalar:.4f}

  This tells us:
  1. How much of the attention computation is DIRECTED (tournament-like)
  2. How much is UNDIRECTED similarity (symmetric)
  3. How much is uniform baseline (scalar)

  The TOURNAMENT SECTOR is where our mathematical machinery applies:
  - OCF: H(T) = I(Omega(T), 2) for the thresholded attention tournament
  - Path homology: beta_k detects topological features of directed attention
  - Spectral theory: eigenvalues of B_T reveal rotation structure

  All of these are PARAMETER-FREE invariants that characterize LLM computation
  from pure mathematics alone.
""")

# Final summary: how does anti_frac vary across the model?
print("  LAYER GRADIENT:")
early_anti = np.mean([r['anti_frac'] for r in all_results if r['layer'] < 4])
mid_anti = np.mean([r['anti_frac'] for r in all_results if 4 <= r['layer'] < 8])
late_anti = np.mean([r['anti_frac'] for r in all_results if r['layer'] >= 8])

print(f"    Early (0-3):  {early_anti:.4f}")
print(f"    Middle (4-7): {mid_anti:.4f}")
print(f"    Late (8-11):  {late_anti:.4f}")

if early_anti < late_anti:
    print(f"    → GRADIENT: tournament energy INCREASES through layers")
elif early_anti > late_anti:
    print(f"    → GRADIENT: tournament energy DECREASES through layers")
else:
    print(f"    → No clear gradient")

print(f"\nDone. opus-2026-03-21-S95")
