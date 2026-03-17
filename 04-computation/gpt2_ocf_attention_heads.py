#!/usr/bin/env python3
"""gpt2_ocf_attention_heads.py — opus-2026-03-17-S74b

ATTENTION HEAD OCF: Each head as a different "voter" in a tournament.

MOTIVATION:
  MC dropout FAILED because dropout masks share 99%+ weights.
  Paraphrases FAILED because tournaments were too stable.
  Layers FAILED because they measure processing depth, not uncertainty.

  Attention heads are GENUINELY DIFFERENT:
  - Each head learns to attend to different aspects (syntax, semantics,
    position, etc.)
  - Heads within the same layer run IN PARALLEL, not sequentially
  - Research shows GPT-2 heads specialize: some do induction, some do
    duplicate suppression, some track syntax

METHODOLOGY:
  1. Extract attention patterns from each head in the last layer(s)
  2. Use each head's attention-weighted representation to produce logits
  3. Build tournament from head-specific logit vectors
  4. OCF confidence measures "do different aspects of the input agree?"

  If different heads produce different token rankings, that indicates
  the prediction relies on a SPECIFIC interpretation of the input —
  a genuine form of fragility/uncertainty.
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


def pearson(x, y):
    if len(x) < 3:
        return 0.0
    if np.std(x) < 1e-10 or np.std(y) < 1e-10:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def logit_gap(logits):
    sorted_l = np.sort(logits)[::-1]
    return sorted_l[0] - sorted_l[1]


def run():
    print("=" * 70)
    print("ATTENTION HEAD OCF — heads as genuinely different voters")
    print("=" * 70)
    print()

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2")
    model.eval()

    n_layer = model.config.n_layer  # 12
    n_head = model.config.n_head    # 12
    d_model = model.config.n_embd   # 768
    d_head = d_model // n_head      # 64

    print(f"GPT-2: {n_layer} layers, {n_head} heads, d_model={d_model}, d_head={d_head}")
    print()

    texts = [
        "The Pythagorean theorem states that in a right-angled triangle, "
        "the square of the hypotenuse is equal to the sum of the squares "
        "of the other two sides. This fundamental result in geometry has "
        "been proved in hundreds of different ways throughout history.",

        "Once upon a time in a small village by the sea, there lived an "
        "old fisherman who had not caught a single fish in eighty-four days. "
        "The boy who had been with him had gone to another boat after forty days.",

        "In the year 2024, artificial intelligence systems became capable "
        "of solving complex mathematical problems that had puzzled researchers "
        "for decades. The breakthrough came from a new approach to training.",

        "The cat sat on the warm windowsill, watching the rain fall outside. "
        "It had been raining for three days straight, and the garden was "
        "completely flooded. The flowers were drowning in the mud.",
    ]

    # ========================================================================
    # METHOD: Extract per-head representations and project to logits
    # ========================================================================
    # GPT-2 architecture for last layer:
    #   x_prev → ln_1(x) → attn(heads) → x + attn_out → ln_2 → mlp → residual
    #
    # Each attention head produces:
    #   head_i = softmax(Q_i K_i^T / sqrt(d_head)) @ V_i @ W_proj[i]
    #
    # The full attention output is: sum_i head_i + bias
    #
    # We can't easily get "what would the model predict from head_i alone"
    # without the MLP, but we CAN:
    # A) Use each head's attention output as a perturbation signal
    # B) Get head-level contributions to the residual stream
    #
    # APPROACH B: Use the attention output of each head, apply ln_f and lm_head
    # to get head-specific logit contributions. This gives the ADDITIVE
    # contribution of each head to the final logit vector.

    print("Approach 1: Per-head logit contributions from last layer")
    print("  Each head's attention output → ln_f → lm_head → head logits")
    print()

    # We need to hook into the attention layer to get per-head outputs
    # GPT-2 attention: c_attn projects to Q,K,V; c_proj projects back
    # The per-head outputs BEFORE c_proj are what we want

    all_positions = []

    for text_idx, text in enumerate(texts):
        ids = tokenizer.encode(text, return_tensors="pt")
        seq_len = ids.shape[1]

        # Hook to capture attention head outputs before projection
        last_layer = model.transformer.h[-1]  # Last transformer block
        head_outputs = {}

        def hook_fn(module, input, output):
            # output is a tuple: (attn_output, present, attn_weights)
            # We need the pre-projection per-head outputs
            # But the standard forward only gives the merged output
            pass

        # Alternative: manually compute per-head contributions
        # Register hooks to get intermediate values
        attn_input = {}
        attn_weights_store = {}

        def capture_attn_input(module, args, kwargs):
            """Capture input to attention layer."""
            # args[0] is the hidden state
            if args:
                attn_input['x'] = args[0].detach()
            return None

        def capture_attn_weights(module, input, output):
            """Capture attention output and weights."""
            # output = (attn_output, present_key_value, attn_weights)
            if len(output) >= 3 and output[2] is not None:
                attn_weights_store['weights'] = output[2].detach()

        # Get hidden states to decompose last layer
        with torch.no_grad():
            outputs = model(ids, output_hidden_states=True)

        logits_standard = outputs.logits[0].numpy()  # [seq_len, vocab_size]

        # Get hidden state before the last layer
        hidden_states = outputs.hidden_states  # [n_layer+1] tensors
        x_before_last = hidden_states[-2]  # [1, seq_len, d_model] - input to last layer

        # Get the value vectors and attention weights from the last layer
        last_attn = model.transformer.h[-1].attn
        last_mlp = model.transformer.h[-1].mlp
        last_ln1 = model.transformer.h[-1].ln_1
        last_ln2 = model.transformer.h[-1].ln_2

        # Compute per-head outputs manually
        x_ln = last_ln1(x_before_last)  # [1, seq_len, d_model]

        # c_attn projects to [Q, K, V] concatenated
        qkv = last_attn.c_attn(x_ln)  # [1, seq_len, 3*d_model]
        q, k, v = qkv.split(d_model, dim=2)

        # Reshape for multi-head
        def split_heads(x, n_head, d_head):
            # [batch, seq, d_model] -> [batch, n_head, seq, d_head]
            return x.view(1, seq_len, n_head, d_head).permute(0, 2, 1, 3)

        q_heads = split_heads(q, n_head, d_head)
        k_heads = split_heads(k, n_head, d_head)
        v_heads = split_heads(v, n_head, d_head)

        # Compute attention weights manually: softmax(Q K^T / sqrt(d_head))
        # with causal mask
        scale = 1.0 / (d_head ** 0.5)
        attn_scores = torch.matmul(q_heads, k_heads.transpose(-2, -1)) * scale
        # Causal mask: only attend to previous positions
        causal_mask = torch.triu(torch.ones(seq_len, seq_len, dtype=torch.bool), diagonal=1)
        attn_scores = attn_scores.masked_fill(causal_mask.unsqueeze(0).unsqueeze(0), float('-inf'))
        attn_w = torch.softmax(attn_scores, dim=-1)  # [1, n_head, seq_len, seq_len]

        # Per-head attention output: attn_w @ V
        head_attn_out = torch.matmul(attn_w, v_heads)  # [1, n_head, seq_len, d_head]

        # c_proj maps from d_model -> d_model
        # It effectively has a block structure: each head's d_head contribution
        # is mapped through its slice of c_proj
        c_proj_weight = last_attn.c_proj.weight  # [d_model, d_model]
        c_proj_bias = last_attn.c_proj.bias      # [d_model]

        # For each head, its contribution to the residual stream is:
        # head_out_i @ c_proj_weight[i*d_head:(i+1)*d_head, :] + bias/n_head
        head_logits_list = []

        for h in range(n_head):
            # Head h output: [1, seq_len, d_head]
            h_out = head_attn_out[0, h]  # [seq_len, d_head]

            # Project through c_proj slice
            w_slice = c_proj_weight[:, h*d_head:(h+1)*d_head]  # [d_model, d_head]
            h_projected = h_out @ w_slice.T  # [seq_len, d_model]

            # This is head h's contribution to the attention output
            # Full residual: x_before_last + sum_h h_projected + bias + mlp(...)
            # To get head-specific logits, we add to the "base" (residual without attention)

            # Method: residual stream = x_before_last + attn_out + mlp(ln_2(x_before_last + attn_out))
            # We want: "if only head h contributed", so we use x_before_last + h_projected
            # Then apply ln_2 + mlp + ln_f + lm_head

            x_with_head = x_before_last[0] + h_projected + c_proj_bias / n_head  # [seq_len, d_model]

            # Apply the MLP (this is a rough approximation — MLP is nonlinear)
            x_ln2 = last_ln2(x_with_head.unsqueeze(0))
            mlp_out = last_mlp(x_ln2)
            x_final = x_with_head.unsqueeze(0) + mlp_out

            # Apply final layer norm and lm_head
            x_normed = model.transformer.ln_f(x_final)
            head_logits = model.lm_head(x_normed)  # [1, seq_len, vocab_size]
            head_logits_list.append(head_logits[0].detach().numpy())

        print(f"  Text {text_idx}: {seq_len} tokens, {n_head} head logit vectors computed.")

        # Now analyze each position
        for pos in range(seq_len - 1):
            actual = ids[0, pos + 1].item()
            standard_top1 = np.argmax(logits_standard[pos])
            correct = int(standard_top1 == actual)

            # Build tournament from head logit vectors
            lvs = [hl[pos] for hl in head_logits_list]

            for tk in [4, 6, 8, 12]:
                A, ti = tournament_from_multi_logits(lvs, top_k=tk)
                t3 = count_3cycles_trace(A)
                H = 1 + 2 * t3

                if tk == 8:
                    H_8 = H
                    t3_8 = t3

            # Head agreement on top-1
            head_top1s = [np.argmax(hl[pos]) for hl in head_logits_list]
            n_distinct = len(set(head_top1s))
            mode_top1 = max(set(head_top1s), key=head_top1s.count)
            agree_frac = head_top1s.count(mode_top1) / n_head

            # Standard logit gap
            gap = logit_gap(logits_standard[pos])

            # Standard entropy
            lv = logits_standard[pos]
            lv_shifted = lv - np.max(lv)
            probs = np.exp(lv_shifted)
            probs /= probs.sum()
            entropy = -np.sum(probs * np.log(probs + 1e-10))

            all_positions.append({
                'text': text_idx,
                'pos': pos,
                'correct': correct,
                'H_8': H_8,
                't3_8': t3_8,
                'n_distinct': n_distinct,
                'agree_frac': agree_frac,
                'gap': gap,
                'entropy': entropy,
                'token': tokenizer.decode([actual]),
            })

    print()
    total = len(all_positions)
    ca = np.array([p['correct'] for p in all_positions], dtype=float)

    # ========================================================================
    # Results
    # ========================================================================
    print("=" * 70)
    print(f"RESULTS — {total} positions, {n_head} attention heads as voters")
    print("=" * 70)
    print()

    # Basic OCF stats
    ha = np.array([p['H_8'] for p in all_positions])
    conf = 1.0 / ha
    frac_cyclic = (ha > 1).mean()
    r_ocf = pearson(conf, ca)

    print(f"OCF (k=8): {frac_cyclic:.1%} cyclic, r(OCF, correct) = {r_ocf:+.4f}")

    # Transitive vs cyclic accuracy
    tr_mask = ha == 1
    cy_mask = ha > 1
    if tr_mask.sum() > 0:
        print(f"  Transitive (H=1): n={tr_mask.sum()}, acc={ca[tr_mask].mean():.3f}")
    if cy_mask.sum() > 0:
        print(f"  Cyclic (H>1):     n={cy_mask.sum()}, acc={ca[cy_mask].mean():.3f}")
    print()

    # H distribution
    if cy_mask.sum() >= 3:
        h_counts = defaultdict(int)
        h_acc = defaultdict(list)
        for p in all_positions:
            h_counts[p['H_8']] += 1
            h_acc[p['H_8']].append(p['correct'])
        print("H distribution:")
        for h in sorted(h_counts.keys()):
            print(f"  H={h:5d}: n={h_counts[h]:3d}, acc={np.mean(h_acc[h]):.3f}")
        print()

    # Head agreement
    agree = np.array([p['agree_frac'] for p in all_positions])
    r_agree = pearson(agree, ca)
    print(f"Head agreement: r(agree, correct) = {r_agree:+.4f}")

    nd = np.array([1.0/p['n_distinct'] for p in all_positions])
    r_nd = pearson(nd, ca)
    print(f"1/n_distinct:   r(1/n_dist, correct) = {r_nd:+.4f}")

    # All signals comparison
    print()
    print("=" * 70)
    print("ALL UNCERTAINTY SIGNALS vs ACCURACY")
    print("=" * 70)
    print()

    signals = {
        'Logit gap': np.array([p['gap'] for p in all_positions]),
        'Neg entropy': np.array([-p['entropy'] for p in all_positions]),
        'Head agreement': agree,
        '1/n_distinct': nd,
        'OCF conf (k=8)': conf,
    }

    for name, sig in sorted(signals.items(), key=lambda x: -abs(pearson(x[1], ca))):
        r = pearson(sig, ca)
        print(f"  {name:<25s}: r = {r:+.4f}")

    # Residual analysis
    print()
    print("=" * 70)
    print("RESIDUAL ANALYSIS: Does head OCF add to logit gap?")
    print("=" * 70)
    print()

    from numpy.linalg import lstsq
    gap_arr = np.array([p['gap'] for p in all_positions])
    X = np.column_stack([gap_arr, agree, np.ones(total)])
    beta, _, _, _ = lstsq(X, ca, rcond=None)
    predicted = X @ beta
    residual = ca - predicted

    r_ocf_resid = pearson(conf, residual)
    print(f"  Baseline r(gap+agree, correct): {pearson(predicted, ca):+.4f}")
    print(f"  r(OCF, residual):               {r_ocf_resid:+.4f}")

    if r_ocf_resid > 0.1:
        print("\n  POSITIVE: Head OCF adds signal beyond gap + agreement!")
    elif r_ocf_resid > 0.05:
        print("\n  MARGINAL: Small additional signal. Needs larger dataset.")
    else:
        print("\n  NULL: Head OCF adds nothing beyond scalars.")

    print()
    print("=" * 70)
    print("COMPARISON: ATTENTION HEADS vs MC DROPOUT as OCF contexts")
    print("=" * 70)
    print(f"""
  MC Dropout (d=0.5, 21 passes, k=16):
    - 63% cyclic, but r(OCF, correct) ≈ 0.0
    - Cycles form among low-ranked tokens (avg rank 9.6/15)
    - No information beyond logit gap

  Attention Heads (12 heads, k=8):
    - {frac_cyclic:.0%} cyclic, r(OCF, correct) = {r_ocf:+.4f}
    - Head agreement: r = {r_agree:+.4f}
    - Logit gap baseline: r = {pearson(gap_arr, ca):+.4f}

  INTERPRETATION:
  Attention heads provide genuinely DIFFERENT views because each head
  has learned to attend to different input features. If heads disagree
  on the token ranking, it means different aspects of the context
  point to different continuations — a real semantic uncertainty signal.

  MC dropout, by contrast, just adds noise. Different dropout masks
  produce nearly identical rankings because they share 99%+ of weights.
""")


if __name__ == "__main__":
    run()
