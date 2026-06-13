#!/usr/bin/env python3
"""gpt2_ocf_dropout_fast.py — opus-2026-03-17-S74b

FAST MC Dropout sweep — targeted configs with just 2 texts.

Focuses on the KEY question: does higher dropout create enough
intransitivity for OCF to become useful?
"""

import torch
import torch.nn as nn
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
from collections import defaultdict
import time


def count_3cycles_trace(A):
    A3 = A @ A @ A
    return int(np.trace(A3)) // 3


def tournament_from_multi_logits(logit_vectors, top_k=8):
    avg_logits = np.mean(logit_vectors, axis=0)
    top_indices = np.argsort(-avg_logits)[:top_k]
    n_ctx = len(logit_vectors)
    k = len(top_indices)
    A = np.zeros((k, k), dtype=np.int32)
    for i in range(k):
        for j in range(i+1, k):
            ti, tj = top_indices[i], top_indices[j]
            vi = sum(1 for lv in logit_vectors if lv[ti] > lv[tj])
            if vi > n_ctx // 2:
                A[i][j] = 1
            elif vi < n_ctx // 2:
                A[j][i] = 1
            else:
                A[i][j] = 1 if avg_logits[ti] > avg_logits[tj] else 0
                A[j][i] = 1 - A[i][j]
    return A, top_indices


def set_dropout_rate(model, rate):
    for m in model.modules():
        if isinstance(m, nn.Dropout):
            m.p = rate
            m.train()


def pearson(x, y):
    if np.std(x) < 1e-10 or np.std(y) < 1e-10:
        return 0.0
    return np.corrcoef(x, y)[0, 1]


def run():
    print("=" * 70)
    print("FAST MC DROPOUT SWEEP — Key question: does high dropout help?")
    print("=" * 70)
    print()

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2")

    # Just 2 texts to keep it fast
    texts = [
        "The Pythagorean theorem states that in a right-angled triangle, "
        "the square of the hypotenuse is equal to the sum of the squares "
        "of the other two sides. This fundamental result in geometry has "
        "been proved in hundreds of different ways throughout history.",

        "Once upon a time in a small village by the sea, there lived an "
        "old fisherman who had not caught a single fish in eighty-four days. "
        "The boy who had been with him had gone to another boat after forty days.",
    ]

    # Deterministic predictions
    model.eval()
    det_data = []
    for text in texts:
        ids = tokenizer.encode(text, return_tensors="pt")
        with torch.no_grad():
            out = model(ids)
        det_data.append((ids, out.logits[0].numpy()))

    total_pos = sum(ids.shape[1] - 1 for ids, _ in det_data)
    print(f"Total positions: {total_pos}")
    print()

    configs = [
        (0.1, 11, 8),
        (0.2, 11, 8),
        (0.3, 11, 8),
        (0.5, 11, 8),
        (0.3, 11, 16),
        (0.5, 11, 16),
        (0.5, 21, 16),
    ]

    print(f"{'Config':<25s} | {'%t3>0':>6} | {'r(OCF,c)':>8} | {'r(MC,c)':>8} | "
          f"{'TrAcc':>6} | {'CyAcc':>6} | {'n_cy':>5} | {'time':>5}")
    print("-" * 85)

    for dr, nmc, tk in configs:
        t0 = time.time()
        set_dropout_rate(model, dr)

        all_correct = []
        all_conf = []
        all_t3 = []
        all_mc_agree = []
        all_H = []

        for ids, det_logits in det_data:
            seq_len = ids.shape[1]
            mc_logits = []
            for _ in range(nmc):
                with torch.no_grad():
                    mc_logits.append(model(ids).logits[0].numpy())

            for pos in range(seq_len - 1):
                actual = ids[0, pos + 1].item()
                top1 = np.argmax(det_logits[pos])
                correct = int(top1 == actual)

                lvs = [mc[pos] for mc in mc_logits]
                A, ti = tournament_from_multi_logits(lvs, top_k=tk)
                t3 = count_3cycles_trace(A)
                H = 1 + 2 * t3

                mc_top1s = [np.argmax(lv) for lv in lvs]
                mode = max(set(mc_top1s), key=mc_top1s.count)
                mc_agree = mc_top1s.count(mode) / nmc

                all_correct.append(correct)
                all_conf.append(1.0 / H)
                all_t3.append(t3)
                all_mc_agree.append(mc_agree)
                all_H.append(H)

        ca = np.array(all_correct, dtype=float)
        cfa = np.array(all_conf)
        t3a = np.array(all_t3)
        mca = np.array(all_mc_agree)
        Ha = np.array(all_H)

        frac = (t3a > 0).mean()
        r_ocf = pearson(cfa, ca)
        r_mc = pearson(mca, ca)

        tm = Ha == 1
        cm = Ha > 1
        ta = ca[tm].mean() if tm.sum() > 0 else float('nan')
        cyc_a = ca[cm].mean() if cm.sum() > 0 else float('nan')
        elapsed = time.time() - t0

        label = f"d={dr:.1f} mc={nmc} k={tk}"
        print(f"  {label:<23s} | {frac:>5.1%} | {r_ocf:>+8.4f} | {r_mc:>+8.4f} | "
              f"{ta:>6.3f} | {cyc_a if cm.sum() > 0 else 'N/A':>6} | {cm.sum():>5d} | {elapsed:>5.1f}s")

        # If this config has meaningful cycling, print H distribution
        if cm.sum() >= 5:
            h_counts = defaultdict(int)
            h_acc = defaultdict(list)
            for i in range(len(all_H)):
                h_counts[all_H[i]] += 1
                h_acc[all_H[i]].append(all_correct[i])
            print(f"    H distribution:")
            for h in sorted(h_counts.keys()):
                print(f"      H={h:5d}: n={h_counts[h]:3d}, acc={np.mean(h_acc[h]):.3f}")

    print()
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print("""
  KEY QUESTION: Does higher dropout create enough intransitivity?

  If %t3>0 stays under ~10% even at dropout=0.5 with k=16,
  then MC dropout fundamentally cannot produce enough tournament
  diversity for OCF to work.

  If %t3>0 exceeds ~20% and r(OCF,correct) becomes meaningfully
  positive (+0.1 or more), then OCF adds real value at high dropout.

  EITHER WAY: MC agreement (scalar) is always useful (r ≈ 0.2-0.3).
  The question is whether the TOURNAMENT STRUCTURE adds anything
  beyond the simple scalar agreement measure.
""")


if __name__ == "__main__":
    run()
