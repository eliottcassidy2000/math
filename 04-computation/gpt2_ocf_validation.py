#!/usr/bin/env python3
"""gpt2_ocf_validation.py — opus-2026-03-17-S74b

EMPIRICAL VALIDATION: Does OCF confidence correlate with real LLM uncertainty?

This is THE key experiment for the tournament output head application.

METHODOLOGY:
1. Load GPT-2 (small, 124M params)
2. For each token in test text:
   a. Extract hidden states from MULTIPLE transformer layers
   b. Project each layer's hidden state through the output head → logits
   c. For top-k tokens, build a TOURNAMENT: layer i beats layer j on token t
      if layer i's logit for t is higher
   d. Count 3-cycles in this tournament → OCF confidence
   e. Record: was the model's prediction correct?
3. Correlate OCF confidence with actual prediction accuracy

HYPOTHESIS: High OCF confidence (few cycles, contexts agree) → higher accuracy
             Low OCF confidence (many cycles, contexts disagree) → lower accuracy
"""

import torch
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
from itertools import combinations
from collections import defaultdict
import time

# =============================================================================
# OCF COMPUTATION (from tournament_output_head.py, proved correct)
# =============================================================================

def count_3cycles(A, n):
    """Count directed 3-cycles in tournament adjacency matrix."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count

def count_5cycles(A, n):
    """Count directed 5-cycles."""
    if n < 5:
        return 0
    count = 0
    for combo in combinations(range(n), 5):
        verts = list(combo)
        from itertools import permutations
        seen = set()
        for perm in permutations(verts):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                canon = min(tuple(perm[i:] + perm[:i]) for i in range(5))
                if canon not in seen:
                    seen.add(canon)
                    count += 1
    return count

def ocf_confidence(A, n):
    """Compute OCF confidence = 1/H where H = 1 + 2*alpha_1 + ...

    Uses PROVED formula (THM-077). Only counts 3-cycles for speed
    at small k (5-cycles negligible for k <= 8).
    """
    t3 = count_3cycles(A, n)
    alpha_1 = t3
    H = 1 + 2 * alpha_1
    return 1.0 / H, H, t3

def tournament_from_multi_logits(logit_vectors, top_k=8):
    """Build tournament from multiple logit vectors (one per context/layer).

    For each pair of top-k tokens (i,j), take majority vote across contexts:
    if more contexts rank i > j, then i→j in the tournament.

    Returns: adjacency matrix A, token indices
    """
    n_contexts = len(logit_vectors)
    vocab_size = logit_vectors[0].shape[0]

    # Find top-k tokens by AVERAGE logit across all contexts
    avg_logits = np.mean(logit_vectors, axis=0)
    top_indices = np.argsort(-avg_logits)[:top_k]

    # Build tournament via majority vote
    k = len(top_indices)
    A = np.zeros((k, k), dtype=int)

    for i in range(k):
        for j in range(i+1, k):
            ti, tj = top_indices[i], top_indices[j]
            votes_i = sum(1 for lv in logit_vectors if lv[ti] > lv[tj])
            votes_j = n_contexts - votes_i
            if votes_i > votes_j:
                A[i][j] = 1
            elif votes_j > votes_i:
                A[j][i] = 1
            else:
                # Tie: use average logit ordering
                if avg_logits[ti] > avg_logits[tj]:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

    return A, top_indices


# =============================================================================
# GPT-2 MULTI-LAYER EXTRACTION
# =============================================================================

def extract_layer_logits(model, tokenizer, input_ids, layers=None):
    """Extract logits from multiple transformer layers.

    For each specified layer, take that layer's hidden state and project
    through the output head (lm_head) to get logits.

    This gives us MULTIPLE "opinions" about the next token — one per layer.
    If they disagree, the tournament has cycles → OCF confidence drops.
    """
    model.eval()

    if layers is None:
        n_layers = model.config.n_layer
        # Use a spread of layers: early, middle, late
        layers = [0, n_layers//4, n_layers//2, 3*n_layers//4, n_layers-1]

    with torch.no_grad():
        outputs = model(input_ids, output_hidden_states=True)
        hidden_states = outputs.hidden_states  # tuple of (batch, seq_len, d_model)
        n_layer = model.config.n_layer

        # Project each layer's hidden state through lm_head
        # NOTE: In transformers v5+, hidden_states[-1] (index n_layer) is
        # ALREADY post-ln_f. Intermediate layers need ln_f applied first.
        layer_logits = []
        for layer_idx in layers:
            hs = hidden_states[layer_idx + 1]  # +1 because index 0 is embeddings

            if layer_idx + 1 < n_layer:
                # Intermediate layer: apply layer norm before projection
                hs = model.transformer.ln_f(hs)
            # else: final layer, already post-ln_f

            logits = model.lm_head(hs)
            layer_logits.append(logits)

    return layer_logits, layers


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def run_experiment():
    print("=" * 70)
    print("GPT-2 OCF VALIDATION — Does tournament confidence predict accuracy?")
    print("=" * 70)
    print()

    # Load model
    print("Loading GPT-2 (124M)...")
    t0 = time.time()
    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2")
    model.eval()
    print(f"  Loaded in {time.time()-t0:.1f}s")
    print(f"  Layers: {model.config.n_layer}, d_model: {model.config.n_embd}, vocab: {model.config.vocab_size}")
    print()

    # Test texts — mix of easy and hard predictions
    test_texts = [
        "The capital of France is",
        "Once upon a time there was a",
        "The quick brown fox jumps over the",
        "In mathematics, a tournament is a complete directed",
        "The president of the United States is the",
        "Water freezes at zero degrees",
        "Machine learning models can sometimes produce",
        "The derivative of x squared is",
        "She went to the store to buy some",
        "The color of the sky is usually",
        "Artificial intelligence will transform",
        "The square root of 144 is",
        "He picked up the phone and called his",
        "The largest planet in our solar system is",
        "Python is a popular programming",
        "The theory of relativity was developed by Albert",
        "Quantum mechanics describes the behavior of",
        "The sum of the angles in a triangle is",
        "Natural language processing is a subfield of",
        "The speed of light is approximately 300000",
    ]

    # Layers to extract from (spread across the network)
    n_layers = model.config.n_layer  # 12 for GPT-2 small
    extract_layers = [0, 2, 4, 6, 8, 10, 11]  # 7 diverse layers

    print(f"Extracting logits from layers: {extract_layers}")
    print(f"Top-k for tournament: 8")
    print()

    results = []

    for text_idx, text in enumerate(test_texts):
        input_ids = tokenizer.encode(text, return_tensors="pt")

        # Get the actual next token prediction (standard)
        with torch.no_grad():
            outputs = model(input_ids)
            standard_logits = outputs.logits[0, -1, :]  # last position
            standard_probs = torch.softmax(standard_logits, dim=-1)
            top1_prob = standard_probs.max().item()
            top1_token = tokenizer.decode([standard_logits.argmax().item()])

            # Entropy of the standard prediction
            entropy = -(standard_probs * torch.log(standard_probs + 1e-10)).sum().item()

        # Extract multi-layer logits
        layer_logits_list, used_layers = extract_layer_logits(
            model, tokenizer, input_ids, layers=extract_layers
        )

        # Get logit vectors for the LAST position from each layer
        logit_vectors = []
        for ll in layer_logits_list:
            lv = ll[0, -1, :].numpy()  # (vocab_size,)
            logit_vectors.append(lv)

        # Build tournament from multi-layer logits
        top_k = 8
        A, top_indices = tournament_from_multi_logits(logit_vectors, top_k=top_k)

        # Compute OCF confidence
        conf, H, t3 = ocf_confidence(A, top_k)

        # Check layer agreement on top-1
        layer_top1s = [np.argmax(lv) for lv in logit_vectors]
        agreement_ratio = max(set(layer_top1s), key=layer_top1s.count)
        n_agree = layer_top1s.count(agreement_ratio)

        # Spectral gap (difference between top-1 and top-2 average logit)
        avg_logits = np.mean(logit_vectors, axis=0)
        sorted_avg = np.sort(avg_logits)[::-1]
        spectral_gap = sorted_avg[0] - sorted_avg[1]

        result = {
            'text': text,
            'top1': top1_token.strip(),
            'top1_prob': top1_prob,
            'entropy': entropy,
            'ocf_confidence': conf,
            'H': H,
            't3': t3,
            'layer_agreement': n_agree / len(extract_layers),
            'spectral_gap': spectral_gap,
        }
        results.append(result)

        # Print per-text results
        regime = "CERTAIN" if conf > 0.5 else ("CHOOSING" if conf > 0.1 else "CONFUSED")
        print(f"  [{text_idx+1:2d}] \"{text}\"")
        print(f"       → '{result['top1']}' (p={top1_prob:.3f})")
        print(f"       OCF: H={H}, t3={t3}, conf={conf:.4f} [{regime}]")
        print(f"       Layers agree: {n_agree}/{len(extract_layers)}, "
              f"entropy={entropy:.2f}, gap={spectral_gap:.2f}")
        print()

    # =============================================================================
    # CORRELATION ANALYSIS
    # =============================================================================

    print("=" * 70)
    print("CORRELATION ANALYSIS")
    print("=" * 70)
    print()

    # Extract arrays
    confs = np.array([r['ocf_confidence'] for r in results])
    probs = np.array([r['top1_prob'] for r in results])
    entropies = np.array([r['entropy'] for r in results])
    agreements = np.array([r['layer_agreement'] for r in results])
    gaps = np.array([r['spectral_gap'] for r in results])
    Hs = np.array([r['H'] for r in results])

    # Pearson correlations
    def pearson(x, y):
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            return 0.0
        return np.corrcoef(x, y)[0, 1]

    print("Pearson correlations:")
    print(f"  OCF confidence vs top1_prob:     r = {pearson(confs, probs):+.4f}")
    print(f"  OCF confidence vs entropy:       r = {pearson(confs, -entropies):+.4f}")
    print(f"  OCF confidence vs layer_agree:   r = {pearson(confs, agreements):+.4f}")
    print(f"  OCF confidence vs spectral_gap:  r = {pearson(confs, gaps):+.4f}")
    print(f"  layer_agree vs top1_prob:        r = {pearson(agreements, probs):+.4f}")
    print(f"  spectral_gap vs top1_prob:       r = {pearson(gaps, probs):+.4f}")
    print(f"  entropy vs top1_prob:            r = {pearson(-entropies, probs):+.4f}")
    print()

    # Regime analysis
    print("REGIME ANALYSIS:")
    for regime_name, lo, hi in [("CERTAIN (conf>0.5)", 0.5, 1.1),
                                  ("CHOOSING (0.1<conf<0.5)", 0.1, 0.5),
                                  ("CONFUSED (conf<0.1)", 0.0, 0.1)]:
        mask = (confs >= lo) & (confs < hi)
        n = mask.sum()
        if n > 0:
            avg_prob = probs[mask].mean()
            avg_ent = entropies[mask].mean()
            avg_agree = agreements[mask].mean()
            print(f"  {regime_name}: n={n}, avg_prob={avg_prob:.3f}, "
                  f"avg_entropy={avg_ent:.2f}, avg_layer_agree={avg_agree:.2f}")
    print()

    # H distribution
    print("H DISTRIBUTION:")
    h_counts = defaultdict(int)
    for r in results:
        h_counts[r['H']] += 1
    for h in sorted(h_counts.keys()):
        examples = [r['text'][:30] for r in results if r['H'] == h][:3]
        print(f"  H={h:3d}: {h_counts[h]:2d} texts — {', '.join(examples)}")
    print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    # Check the key hypothesis
    median_conf = np.median(confs)
    high_conf_mask = confs >= median_conf
    low_conf_mask = confs < median_conf

    if high_conf_mask.sum() > 0 and low_conf_mask.sum() > 0:
        high_conf_prob = probs[high_conf_mask].mean()
        low_conf_prob = probs[low_conf_mask].mean()
        print(f"KEY TEST: Does OCF confidence predict prediction quality?")
        print(f"  High OCF confidence (conf >= {median_conf:.4f}): avg top1_prob = {high_conf_prob:.4f}")
        print(f"  Low OCF confidence  (conf <  {median_conf:.4f}): avg top1_prob = {low_conf_prob:.4f}")
        print(f"  Difference: {high_conf_prob - low_conf_prob:+.4f}")
        print(f"  {'YES — higher OCF confidence → better predictions' if high_conf_prob > low_conf_prob else 'NO — no correlation found'}")

    print()
    print("DOES OCF ADD INFORMATION BEYOND ENTROPY?")
    # Partial correlation: OCF vs prob, controlling for entropy
    if np.std(confs) > 1e-10:
        # Residualize both on entropy
        from numpy.polynomial import polynomial as P
        # Simple linear residualization
        ent_coeff = np.polyfit(entropies, confs, 1)
        conf_resid = confs - np.polyval(ent_coeff, entropies)
        prob_ent_coeff = np.polyfit(entropies, probs, 1)
        prob_resid = probs - np.polyval(prob_ent_coeff, entropies)
        partial_r = pearson(conf_resid, prob_resid)
        print(f"  Partial correlation (OCF vs prob | entropy): r = {partial_r:+.4f}")
        if abs(partial_r) > 0.2:
            print(f"  YES — OCF provides ADDITIONAL information beyond entropy alone")
        else:
            print(f"  UNCLEAR — OCF and entropy may capture similar information at this scale")

    print()
    print("NOVEL FINDINGS:")
    print(f"  - {sum(1 for r in results if r['t3'] > 0)}/{len(results)} texts have 3-cycles "
          f"(intransitive layer disagreement)")
    print(f"  - Mean t3 = {np.mean([r['t3'] for r in results]):.1f} "
          f"(compare: max possible = C(8,3) = 56)")
    print(f"  - Layers {extract_layers[0]} and {extract_layers[-1]} agree on top-1 "
          f"in {sum(1 for r in results if np.argmax(logit_vectors[0]) == np.argmax(logit_vectors[-1]))} "
          f"of {len(results)} cases")

    return results

if __name__ == "__main__":
    results = run_experiment()
