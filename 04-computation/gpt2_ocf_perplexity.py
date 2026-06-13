#!/usr/bin/env python3
"""gpt2_ocf_perplexity.py — opus-2026-03-17-S74b

REFINED EXPERIMENT: OCF confidence vs actual next-token prediction accuracy.

KEY INSIGHT FROM FIRST EXPERIMENT:
  OCF measures ORDERING CONSENSUS across layers — whether different parts
  of the network agree on the relative ranking of top candidates.
  This is DIFFERENT from entropy (how peaked the distribution is).

  The right question is NOT "does OCF correlate with top1 probability?"
  but "does OCF predict whether the model will be CORRECT?"

METHODOLOGY:
1. Run GPT-2 on a text where we know the ground truth (the text itself)
2. For each position, predict next token
3. Check if prediction matches actual next token
4. Compare accuracy in high-OCF vs low-OCF regimes
5. Also: does OCF add information beyond what entropy already provides?

If OCF detects cases where the model is internally inconsistent,
those should be positions where the model is more likely to be WRONG.
"""

import torch
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
from itertools import combinations
from collections import defaultdict
import time

# =============================================================================
# OCF (fast version — 3-cycles only, sufficient for k=8)
# =============================================================================

def count_3cycles_fast(A, n):
    """Count directed 3-cycles using matrix trace: t3 = (Tr(A^3) - 3*sum(s_i(s_i-1)))/6."""
    # Matrix method — much faster for larger n
    A2 = A @ A
    A3 = A2 @ A
    trace = np.trace(A3)
    # Score sequence adjustment
    scores = A.sum(axis=1)
    score_correction = np.sum(scores * (scores - 1))
    # Wrong — the trace method gives t3 = Tr(A^3)/3 for a simple directed graph
    # For tournaments: Tr(A^3) counts each directed 3-cycle 3 times (once per starting vertex)
    # But we also get non-cycle closed walks of length 3
    # Correct: t3 = (Tr(A^3) - sum(s_i * (n-1-s_i) where the walk revisits)
    # Actually for tournaments: Tr(A^3) = 3*t3 exactly! (no 2-walks possible since
    # tournaments have exactly one edge between each pair)
    # Wait — A^3[i][i] = number of walks i→j→k→i of length 3
    # In a tournament, this counts directed 3-cycles through i
    # Each 3-cycle {i,j,k} is counted 3 times in Tr(A^3) (once as i→j→k→i, etc.)
    # Actually no — each cycle i→j→k→i is counted once in A^3[i][i]
    # And there are two orientations: i→j→k→i and i→k→j→i
    # Only one of them exists in a tournament
    # So Tr(A^3) = sum over all directed 3-cycles of 3 (counted at each vertex)
    # Wait — a directed 3-cycle i→j→k→i contributes:
    #   to A^3[i][i]: the walk i→j, j→k, k→i — yes, contributes 1
    #   to A^3[j][j]: the walk j→k, k→i, i→j — yes, contributes 1
    #   to A^3[k][k]: the walk k→i, i→j, j→k — yes, contributes 1
    # So Tr(A^3) = 3 * (number of directed 3-cycles) exactly
    # But in an undirected 3-cycle on vertices {i,j,k}, there are TWO directed cycles
    # In a tournament, exactly ONE directed cycle exists per triple
    # SO: Tr(A^3) = 3 * t3 where t3 = number of cyclic triples
    t3 = trace // 3
    return int(t3)

def tournament_from_multi_logits(logit_vectors, top_k=8):
    """Build tournament from multiple logit vectors via majority vote."""
    avg_logits = np.mean(logit_vectors, axis=0)
    top_indices = np.argsort(-avg_logits)[:top_k]
    n_contexts = len(logit_vectors)
    k = len(top_indices)
    A = np.zeros((k, k), dtype=int)

    for i in range(k):
        for j in range(i+1, k):
            ti, tj = top_indices[i], top_indices[j]
            votes_i = sum(1 for lv in logit_vectors if lv[ti] > lv[tj])
            if votes_i > n_contexts // 2:
                A[i][j] = 1
            elif votes_i < n_contexts // 2:
                A[j][i] = 1
            else:
                if avg_logits[ti] > avg_logits[tj]:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

    return A, top_indices


def extract_all_layer_logits(model, input_ids, layers):
    """Extract logits from specified layers for ALL positions at once.

    NOTE: In transformers v5+, hidden_states[-1] (= hidden_states[n_layer])
    is ALREADY post-ln_f (same as last_hidden_state). Intermediate layers
    (hidden_states[1] through hidden_states[n_layer-1]) are pre-ln_f.
    So we apply ln_f to intermediates but NOT to the final layer.
    """
    model.eval()
    n_layer = model.config.n_layer  # 12 for GPT-2

    with torch.no_grad():
        outputs = model(input_ids, output_hidden_states=True)
        hidden_states = outputs.hidden_states

        layer_logits = []
        for layer_idx in layers:
            hs = hidden_states[layer_idx + 1]
            # hidden_states[n_layer] = last entry = post-ln_f
            # hidden_states[1..n_layer-1] = pre-ln_f (need normalization)
            if layer_idx + 1 < n_layer:
                # Intermediate layer: apply ln_f before projection
                hs = model.transformer.ln_f(hs)
            # else: final layer, already post-ln_f
            logits = model.lm_head(hs)
            layer_logits.append(logits[0].numpy())  # (seq_len, vocab)

    return layer_logits


def run_perplexity_experiment():
    print("=" * 70)
    print("GPT-2 OCF PERPLEXITY — Does intransitivity predict errors?")
    print("=" * 70)
    print()

    # Load
    print("Loading GPT-2...")
    t0 = time.time()
    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2")
    model.eval()
    print(f"  Loaded in {time.time()-t0:.1f}s")
    print()

    # Test text — use diverse content
    test_texts = [
        "The Pythagorean theorem states that in a right-angled triangle, "
        "the square of the hypotenuse is equal to the sum of the squares "
        "of the other two sides. This fundamental result in geometry has "
        "been proved in hundreds of different ways throughout history.",

        "Once upon a time in a small village by the sea, there lived an "
        "old fisherman who had not caught a single fish in eighty-four days. "
        "The boy who had been with him had gone to another boat after forty days.",

        "In computer science, a hash table is a data structure that implements "
        "an associative array abstract data type, a structure that can map keys "
        "to values. A hash table uses a hash function to compute an index.",

        "The quick brown fox jumped over the lazy sleeping dog and ran across "
        "the green meadow towards the distant hills. The sun was setting and "
        "the sky was painted in shades of orange and purple.",
    ]

    layers = [0, 2, 4, 6, 8, 10, 11]
    top_k = 8

    all_results = []

    for text_idx, text in enumerate(test_texts):
        input_ids = tokenizer.encode(text, return_tensors="pt")
        seq_len = input_ids.shape[1]
        tokens = [tokenizer.decode([input_ids[0, i].item()]) for i in range(seq_len)]

        print(f"Text {text_idx+1}: \"{text[:60]}...\" ({seq_len} tokens)")

        # Extract all layer logits at once
        layer_logits = extract_all_layer_logits(model, input_ids, layers)
        # layer_logits[l][pos] = logit vector for layer l at position pos

        # For each position (except last), check prediction
        for pos in range(seq_len - 1):
            actual_next = input_ids[0, pos + 1].item()

            # Get logit vectors from each layer at this position
            logit_vecs = [ll[pos] for ll in layer_logits]

            # Standard prediction (from last layer)
            standard_logits = logit_vecs[-1]  # layer 11
            top1_pred = np.argmax(standard_logits)
            correct = (top1_pred == actual_next)

            # Softmax probability of actual next token
            max_logit = np.max(standard_logits)
            exp_logits = np.exp(standard_logits - max_logit)
            probs = exp_logits / np.sum(exp_logits)
            actual_prob = probs[actual_next]
            top1_prob = probs[top1_pred]

            # Entropy
            entropy = -np.sum(probs * np.log(probs + 1e-10))

            # Build tournament from multi-layer logits
            A, top_indices = tournament_from_multi_logits(logit_vecs, top_k=top_k)
            t3 = count_3cycles_fast(A, top_k)
            H = 1 + 2 * t3
            conf = 1.0 / H

            # Layer agreement on top-1
            layer_top1s = [np.argmax(lv) for lv in logit_vecs]
            n_agree = max(set(layer_top1s), key=layer_top1s.count)
            n_agree = layer_top1s.count(n_agree)

            # Rank of actual token in standard prediction
            rank = np.sum(standard_logits > standard_logits[actual_next]) + 1

            all_results.append({
                'text_idx': text_idx,
                'pos': pos,
                'token': tokens[pos],
                'next_token': tokens[pos+1] if pos+1 < len(tokens) else '?',
                'correct': correct,
                'rank': rank,
                'top1_prob': top1_prob,
                'actual_prob': actual_prob,
                'entropy': entropy,
                'H': H,
                't3': t3,
                'ocf_conf': conf,
                'layer_agree': n_agree / len(layers),
            })

        print(f"  Processed {seq_len-1} positions")

    print()
    n = len(all_results)
    print(f"Total: {n} predictions")
    print()

    # =============================================================================
    # ANALYSIS
    # =============================================================================

    correct = np.array([r['correct'] for r in all_results])
    confs = np.array([r['ocf_conf'] for r in all_results])
    entropies = np.array([r['entropy'] for r in all_results])
    top1_probs = np.array([r['top1_prob'] for r in all_results])
    actual_probs = np.array([r['actual_prob'] for r in all_results])
    Hs = np.array([r['H'] for r in all_results])
    t3s = np.array([r['t3'] for r in all_results])
    ranks = np.array([r['rank'] for r in all_results])
    agreements = np.array([r['layer_agree'] for r in all_results])

    print("=" * 70)
    print("1. OVERALL STATISTICS")
    print("=" * 70)
    print(f"  Accuracy (top-1): {correct.mean():.4f}")
    print(f"  Mean entropy: {entropies.mean():.2f}")
    print(f"  Mean top1_prob: {top1_probs.mean():.4f}")
    print(f"  Mean OCF confidence: {confs.mean():.4f}")
    print(f"  Mean t3: {t3s.mean():.2f}")
    print(f"  Fraction with t3>0: {(t3s > 0).mean():.4f}")
    print(f"  Mean layer agreement: {agreements.mean():.4f}")
    print()

    print("=" * 70)
    print("2. OCF REGIME ANALYSIS — Does intransitivity predict errors?")
    print("=" * 70)
    print()

    regime_bins = [
        ("TRANSITIVE (H=1)", Hs == 1),
        ("MILD CYCLING (H=3-5)", (Hs >= 3) & (Hs <= 5)),
        ("MODERATE CYCLING (H=7-15)", (Hs >= 7) & (Hs <= 15)),
        ("HEAVY CYCLING (H>15)", Hs > 15),
    ]

    print(f"  {'Regime':<30s} | {'N':>4} | {'Acc':>6} | {'AvgProb':>7} | {'AvgEnt':>6} | {'AvgRank':>7} | {'LayAgr':>6}")
    print(f"  {'-'*80}")
    for name, mask in regime_bins:
        if mask.sum() > 0:
            print(f"  {name:<30s} | {mask.sum():>4d} | {correct[mask].mean():>6.3f} | "
                  f"{actual_probs[mask].mean():>7.4f} | {entropies[mask].mean():>6.2f} | "
                  f"{ranks[mask].mean():>7.1f} | {agreements[mask].mean():>6.3f}")
    print()

    # Also try: transitive vs ANY cycling
    trans_mask = Hs == 1
    cycle_mask = Hs > 1
    if trans_mask.sum() > 0 and cycle_mask.sum() > 0:
        print("  BINARY SPLIT:")
        print(f"    Transitive (H=1):  n={trans_mask.sum():4d}, accuracy={correct[trans_mask].mean():.4f}, "
              f"avg_rank={ranks[trans_mask].mean():.1f}")
        print(f"    Has cycles (H>1):  n={cycle_mask.sum():4d}, accuracy={correct[cycle_mask].mean():.4f}, "
              f"avg_rank={ranks[cycle_mask].mean():.1f}")
        acc_diff = correct[trans_mask].mean() - correct[cycle_mask].mean()
        print(f"    Difference: {acc_diff:+.4f} {'(transitive = more accurate ✓)' if acc_diff > 0 else '(unexpected!)'}")
    print()

    print("=" * 70)
    print("3. H-VALUE DISTRIBUTION")
    print("=" * 70)
    h_counts = defaultdict(int)
    h_accuracy = defaultdict(list)
    for r in all_results:
        h_counts[r['H']] += 1
        h_accuracy[r['H']].append(r['correct'])
    print(f"  {'H':>5} | {'Count':>5} | {'Accuracy':>8} | {'Fraction':>8}")
    print(f"  {'-'*35}")
    for h in sorted(h_counts.keys()):
        acc = np.mean(h_accuracy[h])
        frac = h_counts[h] / n
        print(f"  {h:>5d} | {h_counts[h]:>5d} | {acc:>8.4f} | {frac:>8.4f}")
    print()

    print("=" * 70)
    print("4. CORRELATIONS")
    print("=" * 70)
    def pearson(x, y):
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            return 0.0
        return np.corrcoef(x, y)[0, 1]

    # Point-biserial correlation (correct is binary)
    print(f"  OCF conf vs correct:       r = {pearson(confs, correct.astype(float)):+.4f}")
    print(f"  t3 vs correct:             r = {pearson(t3s, correct.astype(float)):+.4f}")
    print(f"  entropy vs correct:        r = {pearson(entropies, -correct.astype(float)):+.4f}")
    print(f"  layer_agree vs correct:    r = {pearson(agreements, correct.astype(float)):+.4f}")
    print(f"  top1_prob vs correct:      r = {pearson(top1_probs, correct.astype(float)):+.4f}")
    print()
    print(f"  OCF conf vs actual_prob:   r = {pearson(confs, actual_probs):+.4f}")
    print(f"  OCF conf vs entropy:       r = {pearson(confs, entropies):+.4f}")
    print(f"  OCF conf vs layer_agree:   r = {pearson(confs, agreements):+.4f}")
    print(f"  t3 vs entropy:             r = {pearson(t3s, entropies):+.4f}")
    print()

    # Does OCF add info beyond entropy?
    print("=" * 70)
    print("5. DOES OCF ADD INFORMATION BEYOND ENTROPY?")
    print("=" * 70)

    # Bin by entropy terciles, then check OCF within each bin
    ent_sorted = np.argsort(entropies)
    n_bins = 3
    bin_size = n // n_bins

    print(f"\n  Within each entropy bin, does OCF still predict accuracy?")
    print(f"  {'Entropy bin':<20s} | {'N':>4} | {'Trans Acc':>9} | {'Cycle Acc':>9} | {'Diff':>8}")
    print(f"  {'-'*60}")

    for b in range(n_bins):
        lo = b * bin_size
        hi = min((b + 1) * bin_size, n) if b < n_bins - 1 else n
        idx = ent_sorted[lo:hi]

        ent_lo = entropies[idx].min()
        ent_hi = entropies[idx].max()
        label = f"[{ent_lo:.1f}, {ent_hi:.1f}]"

        trans_in_bin = np.array([i for i in idx if Hs[i] == 1])
        cycle_in_bin = np.array([i for i in idx if Hs[i] > 1])

        if len(trans_in_bin) > 0 and len(cycle_in_bin) > 0:
            t_acc = correct[trans_in_bin].mean()
            c_acc = correct[cycle_in_bin].mean()
            diff = t_acc - c_acc
            print(f"  {label:<20s} | {len(idx):>4d} | {t_acc:>9.4f} | {c_acc:>9.4f} | {diff:>+8.4f}")
        else:
            t_acc = correct[trans_in_bin].mean() if len(trans_in_bin) > 0 else float('nan')
            c_acc = correct[cycle_in_bin].mean() if len(cycle_in_bin) > 0 else float('nan')
            print(f"  {label:<20s} | {len(idx):>4d} | {t_acc:>9.4f} | {'N/A':>9s} | {'N/A':>8s}")

    print()

    # =============================================================================
    # EXAMPLE INTERESTING CASES
    # =============================================================================

    print("=" * 70)
    print("6. EXAMPLE CASES — High intransitivity (H ≥ 7)")
    print("=" * 70)
    high_h = [r for r in all_results if r['H'] >= 7]
    for r in high_h[:15]:
        status = "✓" if r['correct'] else "✗"
        print(f"  {status} \"{r['token']}\" → \"{r['next_token']}\" "
              f"H={r['H']:3d} t3={r['t3']:2d} rank={r['rank']:3d} "
              f"p={r['actual_prob']:.4f} ent={r['entropy']:.2f} "
              f"agree={r['layer_agree']:.2f}")
    print()

    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()

    # The key finding
    if trans_mask.sum() > 0 and cycle_mask.sum() > 0:
        trans_acc = correct[trans_mask].mean()
        cycle_acc = correct[cycle_mask].mean()
        trans_prob = actual_probs[trans_mask].mean()
        cycle_prob = actual_probs[cycle_mask].mean()

        if trans_acc > cycle_acc:
            print("  ★ POSITIVE RESULT: Transitive tournaments (all layers agree on ordering)")
            print(f"    predict BETTER next tokens: {trans_acc:.4f} vs {cycle_acc:.4f}")
            print(f"    OCF confidence IS a useful signal of model certainty.")
        else:
            print("  ⊘ NULL RESULT: OCF confidence does NOT predict accuracy in this setup.")
            print("    Possible reasons:")
            print("    1. Early layers may disagree for DIFFERENT reasons than late-layer uncertainty")
            print("    2. Layer 0 (embeddings) may carry fundamentally different information")
            print("    3. The tournament among 8 tokens may be too coarse")
            print("    4. Need context-diverse prompts (paraphrases), not layer-diverse logits")

        print()
        print("  QUANTITATIVE RESULT:")
        print(f"    Transitive: {trans_mask.sum()} positions, accuracy={trans_acc:.4f}, avg_prob={trans_prob:.4f}")
        print(f"    Intransitive: {cycle_mask.sum()} positions, accuracy={cycle_acc:.4f}, avg_prob={cycle_prob:.4f}")

    return all_results


if __name__ == "__main__":
    results = run_perplexity_experiment()
