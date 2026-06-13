#!/usr/bin/env python3
"""ocf_llm_synthesis.py — opus-2026-03-17-S74b

SYNTHESIS: What we learned about OCF confidence for LLMs.

This script consolidates all experimental findings into a clear picture
of where OCF works, where it doesn't, and why.
"""

print("""
===========================================================================
OCF FOR LLMs: EXPERIMENTAL SYNTHESIS
opus-2026-03-17-S74b
===========================================================================

EXPERIMENTS CONDUCTED:
1. Layer-based OCF (7 GPT-2 layers as contexts)
2. MC Dropout OCF (11 passes at 0.1 dropout)
3. MC Dropout Deep (7 configs: dropout 0.1-0.5, MC 11-21, k 8-16)
4. Paraphrase OCF (7 rephrasings per question, 9 test cases)
5. MC Dropout Deep Analysis (4 texts, d=0.5, mc=21, k=16, residual analysis)
6. Attention Head OCF (12 heads from last layer as voters)

===========================================================================
FINDING 1: LAYER DISAGREEMENT ≠ UNCERTAINTY
===========================================================================

Transformer layers are BAD contexts for OCF because they measure
PROCESSING DEPTH, not uncertainty.

  Bug discovery: transformers v5+ returns hidden_states[-1] as
  post-ln_f (= last_hidden_state). Applying ln_f again gives garbage.
  Initial "positive" result (H=1 → 15% acc vs H>1 → 5% acc) was
  an artifact of this double-application.

  Corrected result: H=1 → 39.0% accuracy, H>1 → 50.0% accuracy.
  Layer disagreement actually HELPS — it means more refinement.

  LESSON: Different transformer layers are not "voters" in a tournament.
  They are sequential processing stages. Using them as tournament
  contexts conflates "how much processing happened" with "how uncertain
  the model is."

===========================================================================
FINDING 2: MC DROPOUT CREATES INSUFFICIENT DIVERSITY
===========================================================================

GPT-2's default dropout (0.1) produces almost no tournament intransitivity.

  Initial experiment: Only 8/173 positions (4.6%) had t3 > 0.
  The 11 MC passes produce nearly identical logit vectors.

  Key metrics:
  - MC agreement (scalar): r = +0.31 with correct (useful!)
  - OCF confidence (tournament): r = -0.04 with correct (useless)
  - Entropy: r = +0.47 with correct (strong baseline)
  - Top1 prob: r = +0.56 with correct (best baseline)

  WHY: Dropout creates SMALL perturbations, not genuinely different
  viewpoints. The top-k ordering is robust to 10% random neuron
  dropping. OCF needs QUALITATIVE differences, not quantitative noise.

  DEEP SWEEP RESULTS (7 configs, 2+4 texts):
  - d=0.1 k=8:  4.2% cyclic, r(OCF,correct) = -0.02
  - d=0.2 k=8:  10.5% cyclic, r = +0.06
  - d=0.3 k=8:  6.3% cyclic, r = +0.05
  - d=0.5 k=8:  11.6% cyclic, r = -0.03
  - d=0.3 k=16: 47.4% cyclic, r = -0.05
  - d=0.5 k=16: 48.4% cyclic, r = -0.10
  - d=0.5 mc=21 k=16: 63.2% cyclic, r = +0.08

  CRITICAL: Even at 63% intransitivity, OCF does NOT predict accuracy.
  Cycles form among LOW-RANKED tokens (avg rank 9.6/15 in k=16 tournament).
  Only 1% of cycles involve top-3 tokens. They are tied-loser noise.

  RESIDUAL ANALYSIS: After regressing out logit gap + MC agreement,
  OCF adds ZERO incremental information (r = -0.04 to -0.12).

===========================================================================
FINDING 3: PARAPHRASE TOP-1 VARIES, BUT TOURNAMENT IS STILL TRANSITIVE
===========================================================================

Even when different paraphrases produce different top-1 predictions
(e.g., 4 distinct top-1 for "Capital of France"), the tournament
among the top-8 candidates is H=1 for ALL 9 test cases.

  This is because:
  1. With 7 prompts, you need 4+ to disagree on a SPECIFIC PAIR
     to create a cycle via majority vote
  2. Even when top-1 changes, the pairwise ordering among tokens 2-8
     is stable across paraphrases
  3. The signal is in WHETHER top-1 changes, not in the full ordering

  USEFUL SIGNAL: "distinct top-1 count" across paraphrases
    Factual (known answer): avg 2.1 distinct top-1s
    Ambiguous (no single answer): avg 3.0 distinct top-1s
  This captures framing sensitivity without needing OCF.

===========================================================================
FINDING 4: WHAT ACTUALLY WORKS FROM THE PROJECT
===========================================================================

A. STAGED EVALUATION (hot-256 early exit)
   STATUS: WORKS. Drop-in replacement for LLM output head.
   Prediction match: 7/8 with default cache, 8/8 with warmed cache.
   Theoretical savings: 500x FLOPS on confident tokens (V=128K).
   This does NOT depend on OCF — it's pure engineering.

B. OCF FOR RANKING PROBLEMS (non-LLM)
   STATUS: MATHEMATICALLY PROVEN. The OCF formula H = 1 + 2α₁ + 4α₂ + ...
   gives EXACT Hamiltonian path counts. For Chatbot Arena, A/B testing,
   and other ranking applications, OCF provides the mathematically exact
   measure of ranking ambiguity. No empirical validation needed.

C. ARCTANH ATTENTION (formal group)
   STATUS: MATHEMATICALLY CORRECT. arctanh linearizes the formal group
   F(x,y) = (x+y)/(1+xy). This enables additive evidence accumulation.
   Needs empirical validation on a real training run.

D. WALSH SPARSITY
   STATUS: PROVED. H has O(n²) nonzero Walsh coefficients out of 2^{C(n,2)}.
   Compression factor: 341x at n=5, ~100,000x at n=7.
   Application to LLM output compression needs investigation.

E. OCF AS LLM UNCERTAINTY SIGNAL
   STATUS: NULL across ALL tested approaches.
   - Layer-based: NULL (measures processing depth, not uncertainty)
   - MC dropout (default 0.1): NULL (insufficient diversity, 4% cyclic)
   - MC dropout (extreme 0.5, k=16): NULL (63% cyclic, but noise in low ranks)
   - Paraphrase: NULL (tournaments too stable, all H=1)
   - Attention heads: NULL (0% cyclic — MLP dominates head contributions)
   CONCLUSION: Single-model OCF for LLM uncertainty is a dead end.
   Remaining possibility: multi-MODEL ensembles (genuinely different architectures).

===========================================================================
THE FUNDAMENTAL ISSUE
===========================================================================

OCF detects QUALITATIVE disagreement: does the tournament have CYCLES?
This is a binary/coarse signal: either all contexts agree on the
ordering (H=1, transitive) or they don't (H>1, intransitive).

For LLM uncertainty, the signal is typically QUANTITATIVE: how MUCH
do the logits differ? Scalar measures like MC agreement, entropy,
and top-1 probability capture this more directly.

OCF shines when the contexts are GENUINELY DIFFERENT VIEWPOINTS:
- Different human judges (Chatbot Arena)
- Different experiments (A/B testing)
- Different models in an ensemble
- Different prompt framings (but needs MORE paraphrases, not 7)

OCF does NOT add value when contexts are PERTURBATIONS OF THE SAME VIEW:
- Different dropout masks (small noise, even at 50% dropout)
- Different transformer layers (sequential stages, not independent)
- Few paraphrases of the same prompt (insufficient diversity)
- Different attention heads (MLP dominates, rankings are identical)

===========================================================================
RECOMMENDATION
===========================================================================

FOR THE REPO'S ENGINEERING MANDATE:

1. SHIP the staged evaluation (TournamentHead) — it's a real win
2. SHIP the OCF for ranking applications (Chatbot Arena toolkit)
3. KEEP the arctanh attention for future training experiments
4. DO NOT claim OCF detects hallucination in single-model LLMs — REFUTED
5. Multi-MODEL ensemble OCF is the remaining viable approach
6. For single-model uncertainty: use entropy (r=+0.47) or logit gap (r=+0.40)

FOR KIND-PASTEUR:
- The "5 polynomial coefficients" framework (A_0-A_4) uses only 3-cycles
  instead of all odd cycles → 53.1% error rate at n=5, 98.4% at n=7
- The arbitrary weights (9, 4, 7/3, 3/2) are not theorem-backed
- The "regime 30/36/42" thresholds are creative but arbitrary
- The arctanh linearization IS correct and genuinely novel
- The staged evaluation IS a real engineering win
""")
