#!/usr/bin/env python3
"""
gpt2_kurtosis_s97.py — opus-2026-03-21-S97

SPECTRAL KURTOSIS OF GPT-2 ATTENTION PATTERNS

Extends the Cartan probe (S95c) with spectral kurtosis analysis.

For each attention head, we compute:
1. The antisymmetric part A_anti = (A - A^T)/2
2. Its eigenvalues (purely imaginary: ±iθ_k)
3. The spectral kurtosis: how flat is the eigenvalue distribution?
4. Compare across layers — do deeper layers have flatter spectra?

HYPOTHESIS: Heads that act as "balanced reasoners" (considering multiple
tokens equally) will have LOW kurtosis (flat spectrum = Paley-like).
Heads that act as "focused spotlighters" (attending to one token) will
have HIGH kurtosis (peaked spectrum = transitive-like).

Author: opus-2026-03-21-S97
"""

import sys
import numpy as np
import torch
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  GPT-2 SPECTRAL KURTOSIS ANALYSIS")
print("  opus-2026-03-21-S97")
print("=" * 72)

import warnings
warnings.filterwarnings('ignore')
from transformers import GPT2LMHeadModel, GPT2Tokenizer

tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
model = GPT2LMHeadModel.from_pretrained("gpt2", output_attentions=True)
model.eval()

sentences = [
    "The cat sat on the mat and looked at the dog.",
    "Mathematical proofs require careful reasoning and attention to detail.",
    "One plus one equals two, and two plus two equals four.",
    "The quick brown fox jumps over the lazy dog near the river.",
    "In tournament theory, every pair of vertices has exactly one directed edge.",
]

print(f"  Model: GPT-2 small, {len(sentences)} sentences\n")

def spectral_kurtosis(M):
    """Compute spectral kurtosis of the antisymmetric part of M.
    Returns (kurtosis, theta_magnitudes, anti_frac, tr_A4, tr_A2).
    """
    n = M.shape[0]
    A_anti = (M - M.T) / 2.0

    # Eigenvalues of skew-symmetric matrix
    eigs = np.linalg.eigvals(A_anti)
    thetas = np.abs(eigs.imag)
    thetas = thetas[thetas > 1e-12]

    if len(thetas) == 0:
        return 0, [], 0, 0, 0

    # Spectral moments
    tr_A2 = -np.sum(thetas ** 2)  # = tr(A_anti^2)
    tr_A4 = np.sum(thetas ** 4)   # = tr(A_anti^4) [sign changes cancel]

    # Normalized kurtosis: tr(A^4)/tr(A^2)^2 * n
    kappa = tr_A4 / (tr_A2 ** 2) * n if tr_A2 != 0 else 0

    # Anti fraction
    norm2 = np.sum(M ** 2)
    anti_frac = np.sum(A_anti ** 2) / norm2 if norm2 > 0 else 0

    return kappa, thetas, anti_frac, tr_A4, tr_A2


# ========================================================================
# PART 1: Spectral kurtosis by layer and head
# ========================================================================
print("=" * 72)
print("  PART 1: SPECTRAL KURTOSIS BY LAYER AND HEAD")
print("=" * 72)

all_results = []

for si, sentence in enumerate(sentences):
    inputs = tokenizer(sentence, return_tensors="pt")
    n_tokens = inputs['input_ids'].shape[1]

    with torch.no_grad():
        outputs = model(**inputs)
        attentions = outputs.attentions

    for layer in range(12):
        attn = attentions[layer][0].numpy()
        for head in range(12):
            A = attn[head]
            kappa, thetas, af, tr4, tr2 = spectral_kurtosis(A)

            all_results.append({
                'sentence': si,
                'layer': layer,
                'head': head,
                'n': n_tokens,
                'kurtosis': kappa,
                'anti_frac': af,
                'tr_A4': tr4,
                'tr_A2': tr2,
                'n_thetas': len(thetas),
                'theta_max': max(thetas) if len(thetas) > 0 else 0,
                'theta_min': min(thetas) if len(thetas) > 0 else 0,
                'theta_spread': (max(thetas) - min(thetas)) if len(thetas) > 0 else 0,
            })

# ========================================================================
# PART 2: Layer-by-layer summary
# ========================================================================
print(f"\n{'Layer':>5s} {'κ_mean':>8s} {'κ_std':>8s} {'anti':>8s} "
      f"{'θ_spread':>10s} {'θ_max':>8s}")

for layer in range(12):
    lr = [r for r in all_results if r['layer'] == layer]
    km = np.mean([r['kurtosis'] for r in lr])
    ks = np.std([r['kurtosis'] for r in lr])
    af = np.mean([r['anti_frac'] for r in lr])
    ts = np.mean([r['theta_spread'] for r in lr])
    tm = np.mean([r['theta_max'] for r in lr])
    print(f"{layer:5d} {km:8.4f} {ks:8.4f} {af:8.4f} {ts:10.6f} {tm:8.6f}")

# ========================================================================
# PART 3: Head-by-head summary (averaged across layers/sentences)
# ========================================================================
print(f"\n{'Head':>5s} {'κ_mean':>8s} {'κ_std':>8s} {'anti':>8s}")

for head in range(12):
    hr = [r for r in all_results if r['head'] == head]
    km = np.mean([r['kurtosis'] for r in hr])
    ks = np.std([r['kurtosis'] for r in hr])
    af = np.mean([r['anti_frac'] for r in hr])
    print(f"{head:5d} {km:8.4f} {ks:8.4f} {af:8.4f}")

# ========================================================================
# PART 4: Compare trained vs random kurtosis
# ========================================================================
print(f"\n" + "=" * 72)
print("  PART 4: TRAINED vs RANDOM SPECTRAL KURTOSIS")
print("=" * 72)

# Random baseline
n_test = all_results[0]['n']
rng = np.random.RandomState(42)
rand_kurtoses = []
for _ in range(1000):
    logits = rng.randn(n_test, n_test)
    exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
    A = exp_l / exp_l.sum(axis=1, keepdims=True)
    k, _, _, _, _ = spectral_kurtosis(A)
    rand_kurtoses.append(k)

gpt2_kurtoses = [r['kurtosis'] for r in all_results]
rand_mean = np.mean(rand_kurtoses)
gpt2_mean = np.mean(gpt2_kurtoses)

print(f"\n  Random κ (n={n_test}): {rand_mean:.4f} ± {np.std(rand_kurtoses):.4f}")
print(f"  GPT-2 κ:              {gpt2_mean:.4f} ± {np.std(gpt2_kurtoses):.4f}")
print(f"  Shift: {gpt2_mean - rand_mean:+.4f}")

if gpt2_mean > rand_mean:
    print(f"  → GPT-2 has HIGHER kurtosis (more peaked spectrum)")
    print(f"  → Attention is MORE FOCUSED than random (concentrates on fewer tokens)")
elif gpt2_mean < rand_mean:
    print(f"  → GPT-2 has LOWER kurtosis (flatter spectrum)")
    print(f"  → Attention is MORE BALANCED than random (Paley-like)")
else:
    print(f"  → No significant shift")

# ========================================================================
# PART 5: Kurtosis gradient across layers
# ========================================================================
print(f"\n" + "=" * 72)
print("  PART 5: KURTOSIS GRADIENT ACROSS LAYERS")
print("=" * 72)

early = np.mean([r['kurtosis'] for r in all_results if r['layer'] < 4])
mid = np.mean([r['kurtosis'] for r in all_results if 4 <= r['layer'] < 8])
late = np.mean([r['kurtosis'] for r in all_results if r['layer'] >= 8])

print(f"  Early layers (0-3):  κ = {early:.4f}")
print(f"  Middle layers (4-7): κ = {mid:.4f}")
print(f"  Late layers (8-11):  κ = {late:.4f}")

if early > late:
    print(f"  → Kurtosis DECREASES through layers (early focused → late balanced)")
    print(f"  → Later layers develop more Paley-like attention structure")
elif early < late:
    print(f"  → Kurtosis INCREASES through layers (early balanced → late focused)")
    print(f"  → Later layers develop more focused attention")

# ========================================================================
# PART 6: Kurtosis × anti_frac interaction
# ========================================================================
print(f"\n" + "=" * 72)
print("  PART 6: KURTOSIS × ANTI_FRAC INTERACTION")
print("=" * 72)

# Correlation between kurtosis and anti_fraction
ks = [r['kurtosis'] for r in all_results]
afs = [r['anti_frac'] for r in all_results]
corr = np.corrcoef(ks, afs)[0, 1]

print(f"  corr(κ, anti_frac) = {corr:+.4f}")
if corr > 0.1:
    print(f"  → Higher kurtosis ↔ more tournament energy")
    print(f"  → Focused attention is MORE directional")
elif corr < -0.1:
    print(f"  → Higher kurtosis ↔ less tournament energy")
    print(f"  → Focused attention is LESS directional")
else:
    print(f"  → Weak correlation — kurtosis and directionality are independent")

# Per-layer correlation
print(f"\n  Per-layer corr(κ, anti_frac):")
for layer in range(12):
    lr = [r for r in all_results if r['layer'] == layer]
    ks_l = [r['kurtosis'] for r in lr]
    afs_l = [r['anti_frac'] for r in lr]
    if np.std(ks_l) > 1e-10 and np.std(afs_l) > 1e-10:
        c = np.corrcoef(ks_l, afs_l)[0, 1]
        print(f"    Layer {layer:2d}: {c:+.4f}")

# ========================================================================
# PART 7: THE RANKING QUALITY METRIC ON ATTENTION
# ========================================================================
print(f"\n" + "=" * 72)
print("  PART 7: RANKING QUALITY OF ATTENTION")
print("=" * 72)
print("""
  For each attention head, the thresholded tournament has a "ranking quality":
    Quality = 1 - κ(T)/κ_max

  Quality ≈ 1: Paley-like balanced attention (many alternative orderings)
  Quality ≈ 0: Transitive focused attention (one dominant ordering)
""")

# Compute quality for each layer
print(f"  {'Layer':>5s} {'Quality_mean':>12s} {'Quality_std':>12s}")

for layer in range(12):
    lr = [r for r in all_results if r['layer'] == layer]
    # Approximate κ_max for n tokens
    n = lr[0]['n']
    # For transitive tournament of size n: κ = n * tr(B^4) / tr(B^2)^2
    # tr(B^2) = n(n-1), tr(B^4) from computation...
    # Use the maximum observed kurtosis as κ_max
    all_ks = [r['kurtosis'] for r in all_results]
    k_max = max(all_ks) if max(all_ks) > 0 else 1

    qualities = [1 - r['kurtosis'] / k_max for r in lr]
    qm = np.mean(qualities)
    qs = np.std(qualities)
    print(f"  {layer:5d} {qm:12.4f} {qs:12.4f}")

# ========================================================================
# PART 8: SYNTHESIS
# ========================================================================
print(f"\n" + "=" * 72)
print("  SYNTHESIS: SPECTRAL KURTOSIS IN THE CARTAN BRIDGE")
print("=" * 72)

print(f"""
  GPT-2 SPECTRAL KURTOSIS PROFILE:

  1. OVERALL: GPT-2 κ = {gpt2_mean:.4f} vs Random κ = {rand_mean:.4f}
     → {'Higher' if gpt2_mean > rand_mean else 'Lower'} than random

  2. GRADIENT: Early κ = {early:.4f} → Late κ = {late:.4f}
     → {'Decreasing' if early > late else 'Increasing'} through layers

  3. Combined with S95c findings:
     - anti_frac increases through layers (0.34 → 0.43)
     - kurtosis {'decreases' if early > late else 'increases'} through layers

  INTERPRETATION:
  The spectral kurtosis complements the Cartan fraction.
  Together they give a 2D FINGERPRINT of each attention head:
    (anti_frac, kurtosis) = (how tournament-like, how focused)

  Four quadrants:
    HIGH anti_frac + LOW kurtosis = "Balanced Reasoner" (Paley-like)
    HIGH anti_frac + HIGH kurtosis = "Focused Decider" (decisive)
    LOW anti_frac + LOW kurtosis = "Symmetric Cooperator" (undirected)
    LOW anti_frac + HIGH kurtosis = "Degenerate" (neither directed nor balanced)
""")

print("Done. opus-2026-03-21-S97")
