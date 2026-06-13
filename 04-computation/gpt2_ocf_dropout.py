#!/usr/bin/env python3
"""gpt2_ocf_dropout.py — opus-2026-03-17-S74b

MC DROPOUT OCF: Using dropout passes as tournament contexts.

INSIGHT FROM PREVIOUS EXPERIMENTS:
  Transformer layers are BAD contexts for OCF because they measure
  processing depth (positive signal), not uncertainty.

  MC Dropout gives GENUINE uncertainty estimates: each forward pass
  with dropout enabled samples a different sub-network. If these
  sub-networks disagree on the top-k ranking, that's real uncertainty.

METHODOLOGY:
1. Enable dropout in GPT-2 (set model.train() but disable batch norm)
2. Run K forward passes with dropout
3. Build tournament from K logit vectors via majority vote
4. Compute OCF confidence
5. Check correlation with prediction accuracy
"""

import torch
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
from collections import defaultdict
import time


def count_3cycles_trace(A):
    """Count 3-cycles via Tr(A^3)/3."""
    A3 = A @ A @ A
    return int(np.trace(A3)) // 3


def tournament_from_multi_logits(logit_vectors, top_k=8):
    """Build tournament from multiple logit vectors via majority vote."""
    avg_logits = np.mean(logit_vectors, axis=0)
    top_indices = np.argsort(-avg_logits)[:top_k]
    n_contexts = len(logit_vectors)
    k = len(top_indices)
    A = np.zeros((k, k), dtype=np.int32)

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


def enable_dropout(model):
    """Enable dropout for MC sampling while keeping everything else in eval mode."""
    for module in model.modules():
        if isinstance(module, torch.nn.Dropout):
            module.train()


def run_mc_dropout_experiment():
    print("=" * 70)
    print("GPT-2 OCF MC DROPOUT — Genuine uncertainty via sub-network sampling")
    print("=" * 70)
    print()

    # Load
    print("Loading GPT-2...")
    t0 = time.time()
    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2")
    print(f"  Loaded in {time.time()-t0:.1f}s")
    print(f"  Dropout rate: {model.config.resid_pdrop}")
    print()

    # Enable MC dropout
    model.eval()  # Base state: eval
    enable_dropout(model)  # But enable dropout layers

    # Test texts
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

    n_mc_passes = 11  # Odd number for clean majority votes
    top_k = 8

    all_results = []

    for text_idx, text in enumerate(test_texts):
        input_ids = tokenizer.encode(text, return_tensors="pt")
        seq_len = input_ids.shape[1]
        tokens = [tokenizer.decode([input_ids[0, i].item()]) for i in range(seq_len)]

        print(f"Text {text_idx+1}: \"{text[:60]}...\" ({seq_len} tokens)")

        # Get deterministic prediction (eval mode, no dropout)
        model.eval()
        with torch.no_grad():
            det_outputs = model(input_ids)
            det_logits = det_outputs.logits[0].numpy()  # (seq_len, vocab)

        # MC dropout passes
        enable_dropout(model)
        mc_logits = []
        for mc_pass in range(n_mc_passes):
            with torch.no_grad():
                mc_out = model(input_ids)
                mc_logits.append(mc_out.logits[0].numpy())  # (seq_len, vocab)

        # For each position, compute OCF and check accuracy
        for pos in range(seq_len - 1):
            actual_next = input_ids[0, pos + 1].item()

            # Deterministic prediction
            std_logits = det_logits[pos]
            top1_pred = np.argmax(std_logits)
            correct = (top1_pred == actual_next)

            # Softmax
            max_l = np.max(std_logits)
            probs = np.exp(std_logits - max_l)
            probs /= probs.sum()
            top1_prob = probs[top1_pred]
            actual_prob = probs[actual_next]
            entropy = -np.sum(probs * np.log(probs + 1e-10))

            # MC dropout logit vectors at this position
            logit_vecs = [mc[pos] for mc in mc_logits]

            # Build tournament
            A, top_indices = tournament_from_multi_logits(logit_vecs, top_k=top_k)
            t3 = count_3cycles_trace(A)
            H = 1 + 2 * t3
            conf = 1.0 / H

            # MC agreement: how many passes agree on top-1?
            mc_top1s = [np.argmax(lv) for lv in logit_vecs]
            mode_token = max(set(mc_top1s), key=mc_top1s.count)
            mc_agreement = mc_top1s.count(mode_token) / n_mc_passes

            # MC variance of top-1 logit (uncertainty measure)
            top1_logits_mc = [lv[top1_pred] for lv in logit_vecs]
            mc_variance = np.var(top1_logits_mc)

            rank = np.sum(std_logits > std_logits[actual_next]) + 1

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
                'mc_agreement': mc_agreement,
                'mc_variance': mc_variance,
            })

        print(f"  Processed {seq_len-1} positions, {n_mc_passes} MC passes each")

    print()
    n = len(all_results)
    print(f"Total: {n} predictions × {n_mc_passes} MC passes")
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
    mc_agreements = np.array([r['mc_agreement'] for r in all_results])
    mc_variances = np.array([r['mc_variance'] for r in all_results])

    def pearson(x, y):
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            return 0.0
        return np.corrcoef(x, y)[0, 1]

    print("=" * 70)
    print("1. OVERALL STATISTICS")
    print("=" * 70)
    print(f"  Accuracy (top-1): {correct.mean():.4f}")
    print(f"  Mean OCF confidence: {confs.mean():.4f}")
    print(f"  Mean t3: {t3s.mean():.2f}")
    print(f"  Fraction with t3>0: {(t3s > 0).mean():.4f}")
    print(f"  Mean MC agreement: {mc_agreements.mean():.4f}")
    print(f"  Mean MC variance: {mc_variances.mean():.4f}")
    print()

    print("=" * 70)
    print("2. OCF REGIME ANALYSIS")
    print("=" * 70)
    print()

    regime_bins = [
        ("TRANSITIVE (H=1)", Hs == 1),
        ("MILD CYCLING (H=3-5)", (Hs >= 3) & (Hs <= 5)),
        ("MODERATE+ (H≥7)", Hs >= 7),
    ]

    print(f"  {'Regime':<25s} | {'N':>4} | {'Acc':>6} | {'AvgProb':>7} | {'MC_Agr':>6} | {'MC_Var':>8}")
    print(f"  {'-'*70}")
    for name, mask in regime_bins:
        if mask.sum() > 0:
            print(f"  {name:<25s} | {mask.sum():>4d} | {correct[mask].mean():>6.3f} | "
                  f"{actual_probs[mask].mean():>7.4f} | {mc_agreements[mask].mean():>6.3f} | "
                  f"{mc_variances[mask].mean():>8.4f}")
    print()

    # Binary split
    trans_mask = Hs == 1
    cycle_mask = Hs > 1
    if trans_mask.sum() > 0 and cycle_mask.sum() > 0:
        print("  BINARY SPLIT:")
        print(f"    Transitive (H=1):  n={trans_mask.sum():4d}, accuracy={correct[trans_mask].mean():.4f}, "
              f"mc_agree={mc_agreements[trans_mask].mean():.3f}")
        print(f"    Has cycles (H>1):  n={cycle_mask.sum():4d}, accuracy={correct[cycle_mask].mean():.4f}, "
              f"mc_agree={mc_agreements[cycle_mask].mean():.3f}")
        acc_diff = correct[trans_mask].mean() - correct[cycle_mask].mean()
        print(f"    Accuracy diff: {acc_diff:+.4f} "
              f"{'(transitive = more accurate ✓)' if acc_diff > 0 else '(null result)'}")
    print()

    print("=" * 70)
    print("3. CORRELATIONS")
    print("=" * 70)
    print(f"  OCF conf vs correct:       r = {pearson(confs, correct.astype(float)):+.4f}")
    print(f"  t3 vs correct:             r = {pearson(t3s, correct.astype(float)):+.4f}")
    print(f"  MC agreement vs correct:   r = {pearson(mc_agreements, correct.astype(float)):+.4f}")
    print(f"  MC variance vs correct:    r = {pearson(mc_variances, correct.astype(float)):+.4f}")
    print(f"  entropy vs correct:        r = {pearson(entropies, -correct.astype(float)):+.4f}")
    print(f"  top1_prob vs correct:      r = {pearson(top1_probs, correct.astype(float)):+.4f}")
    print()
    print(f"  OCF conf vs MC agreement:  r = {pearson(confs, mc_agreements):+.4f}")
    print(f"  OCF conf vs MC variance:   r = {pearson(confs, mc_variances):+.4f}")
    print(f"  OCF conf vs entropy:       r = {pearson(confs, entropies):+.4f}")
    print(f"  MC agreement vs entropy:   r = {pearson(mc_agreements, entropies):+.4f}")
    print()

    print("=" * 70)
    print("4. H DISTRIBUTION")
    print("=" * 70)
    h_counts = defaultdict(int)
    h_accuracy = defaultdict(list)
    h_mc_agree = defaultdict(list)
    for r in all_results:
        h_counts[r['H']] += 1
        h_accuracy[r['H']].append(r['correct'])
        h_mc_agree[r['H']].append(r['mc_agreement'])
    print(f"  {'H':>5} | {'Count':>5} | {'Accuracy':>8} | {'MC_Agree':>8}")
    print(f"  {'-'*35}")
    for h in sorted(h_counts.keys()):
        acc = np.mean(h_accuracy[h])
        agr = np.mean(h_mc_agree[h])
        print(f"  {h:>5d} | {h_counts[h]:>5d} | {acc:>8.4f} | {agr:>8.4f}")
    print()

    # Does OCF add info beyond entropy? (within entropy bins)
    print("=" * 70)
    print("5. DOES MC-DROPOUT OCF ADD INFORMATION BEYOND ENTROPY?")
    print("=" * 70)

    ent_sorted = np.argsort(entropies)
    n_bins = 3
    bin_size = n // n_bins

    print(f"\n  {'Entropy bin':<20s} | {'N':>4} | {'Trans Acc':>9} | {'Cycle Acc':>9} | {'Diff':>8}")
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
            t_n = len(trans_in_bin) if len(trans_in_bin) > 0 else 0
            c_n = len(cycle_in_bin) if len(cycle_in_bin) > 0 else 0
            t_acc = correct[trans_in_bin].mean() if t_n > 0 else float('nan')
            c_acc = correct[cycle_in_bin].mean() if c_n > 0 else float('nan')
            print(f"  {label:<20s} | {len(idx):>4d} | {str(f'{t_acc:.4f}') if t_n > 0 else 'N/A':>9s} (n={t_n}) | "
                  f"{str(f'{c_acc:.4f}') if c_n > 0 else 'N/A':>9s} (n={c_n})")
    print()

    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()

    if trans_mask.sum() > 0 and cycle_mask.sum() > 0:
        trans_acc = correct[trans_mask].mean()
        cycle_acc = correct[cycle_mask].mean()

        if trans_acc > cycle_acc + 0.02:
            print("  ★ POSITIVE RESULT: MC dropout-based OCF predicts accuracy!")
            print(f"    Transitive: {trans_acc:.4f}, Intransitive: {cycle_acc:.4f}")
            print(f"    OCF adds genuine uncertainty information from sub-network sampling.")
        elif cycle_acc > trans_acc + 0.02:
            print("  ⊖ REVERSE RESULT: Intransitive MC dropout = MORE accurate")
            print(f"    Transitive: {trans_acc:.4f}, Intransitive: {cycle_acc:.4f}")
            print(f"    MC dropout may be too mild at GPT-2's default 0.1 dropout rate.")
        else:
            print("  ⊘ NULL RESULT: No clear relationship between MC dropout OCF and accuracy.")
            print(f"    Transitive: {trans_acc:.4f}, Intransitive: {cycle_acc:.4f}")
            print(f"    GPT-2's 0.1 dropout rate may be too low to create meaningful")
            print(f"    sub-network diversity. Try higher dropout or ensemble approach.")

    # Additional analysis: is MC agreement a better predictor?
    print()
    print("  COMPARISON OF PREDICTORS:")
    print(f"    top1_prob vs correct:    r = {pearson(top1_probs, correct.astype(float)):+.4f} (best baseline)")
    print(f"    entropy vs correct:      r = {pearson(entropies, -correct.astype(float)):+.4f}")
    print(f"    MC agreement vs correct: r = {pearson(mc_agreements, correct.astype(float)):+.4f}")
    print(f"    OCF conf vs correct:     r = {pearson(confs, correct.astype(float)):+.4f}")
    print(f"    MC variance vs correct:  r = {pearson(mc_variances, correct.astype(float)):+.4f}")

    return all_results


if __name__ == "__main__":
    results = run_mc_dropout_experiment()
